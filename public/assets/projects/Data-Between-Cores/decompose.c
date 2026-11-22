/**
 * Toy domain decomposition + renumbering benchmark (no CFD, runs in ~seconds).
 *
 * Mesh: 2D L-shaped domain with a circular obstacle removed, cells numbered in
 * random order (like a mesh generator that does not care about locality).
 * The mesh is only a graph: cells + face neighbours (CSR).
 *
 * Decompositions into P parts (simplified versions of the real algorithms):
 *   random  : every cell goes to a random part
 *   simple  : geometric blocks px x py, like OpenFOAM's "simple" method
 *   scotch* : recursive graph bisection (BFS from a pseudo-peripheral cell,
 *             split the BFS order in half, recurse) + boundary refinement
 *   metis*  : direct k-way graph growing from far-apart seeds + boundary
 *             refinement (no multilevel coarsening, which real METIS has)
 *
 * Renumbering: sweep a Jacobi-like kernel over the graph with the original
 * (random) numbering vs a BFS (Cuthill-McKee-like) numbering.
 *
 * Build: gcc -O3 -march=native decompose.c -o decompose
 * Run:   ./decompose [P] [outdir-for-ppm-images]
 */
#define _POSIX_C_SOURCE 200112L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#ifndef NX
#define NX 1600
#define NY 800
#endif
#define CX (NX / 4.0)
#define CY (NY / 2.0)
#define R  (NX * 0.06)
/* L-shaped domain: the whole bottom half plus the left half of the top */
#define INSIDE(i, j) ((j) < NY / 2 || (i) < NX / 2)

static int n;                 /* number of cells */
static int *xadj, *adj;       /* CSR adjacency */
static int *cx, *cy;          /* cell coordinates (for "simple" and images) */
static int cell_at[NY][NX];   /* -1 = obstacle */

static double now(void) {
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}
static unsigned long long rng = 88172645463325252ULL;
static unsigned rnd(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return (unsigned)rng; }

/* ---------------- mesh ---------------- */
static void build_mesh(void) {
    int count = 0;
    for (int j = 0; j < NY; j++)
        for (int i = 0; i < NX; i++) {
            double dx = i + 0.5 - CX, dy = j + 0.5 - CY;
            cell_at[j][i] = (!INSIDE(i, j) || dx * dx + dy * dy < R * R) ? -1 : count++;
        }
    n = count;
    /* random numbering: shuffle ids */
    int *perm = malloc(n * sizeof(int));
    for (int k = 0; k < n; k++) perm[k] = k;
    for (int k = n - 1; k > 0; k--) { int r = rnd() % (k + 1), t = perm[k]; perm[k] = perm[r]; perm[r] = t; }
    cx = malloc(n * sizeof(int)); cy = malloc(n * sizeof(int));
    for (int j = 0; j < NY; j++)
        for (int i = 0; i < NX; i++)
            if (cell_at[j][i] >= 0) {
                int id = perm[cell_at[j][i]];
                cell_at[j][i] = id; cx[id] = i; cy[id] = j;
            }
    free(perm);
    xadj = malloc((n + 1) * sizeof(int)); adj = malloc(4 * (size_t)n * sizeof(int));
    int e = 0;
    for (int c = 0; c < n; c++) {
        xadj[c] = e;
        int i = cx[c], j = cy[c];
        const int di[4] = {1, -1, 0, 0}, dj[4] = {0, 0, 1, -1};
        for (int d = 0; d < 4; d++) {
            int a = i + di[d], b = j + dj[d];
            if (a >= 0 && a < NX && b >= 0 && b < NY && cell_at[b][a] >= 0) adj[e++] = cell_at[b][a];
        }
    }
    xadj[n] = e;
}

/* BFS over cells with mask[c]==tag, from src. Writes visit order; returns count.
   dist (optional) gets BFS distance. */
static int bfs(int src, const int *mask, int tag, int *order, int *dist, int *seen, int stamp) {
    int head = 0, tail = 0;
    order[tail++] = src; seen[src] = stamp; if (dist) dist[src] = 0;
    while (head < tail) {
        int c = order[head++];
        for (int e = xadj[c]; e < xadj[c + 1]; e++) {
            int nb = adj[e];
            if (seen[nb] != stamp && (!mask || mask[nb] == tag)) {
                seen[nb] = stamp; if (dist) dist[nb] = dist[c] + 1; order[tail++] = nb;
            }
        }
    }
    return tail;
}

static int *seen, stamp_ctr = 0;

/* pseudo-peripheral node: BFS twice, take the last cell reached */
static int peripheral(int start, const int *mask, int tag, int *buf) {
    int cnt = bfs(start, mask, tag, buf, NULL, seen, ++stamp_ctr);
    int far = buf[cnt - 1];
    cnt = bfs(far, mask, tag, buf, NULL, seen, ++stamp_ctr);
    return buf[cnt - 1];
}

/* ---------------- decompositions ---------------- */
static void dec_random(int *part, int P) { for (int c = 0; c < n; c++) part[c] = rnd() % P; }

static void dec_simple(int *part, int P) {
    int px = P, py = 1;                      /* pick px*py = P closest to the aspect ratio */
    for (int a = 1; a <= P; a++) if (P % a == 0) {
        int b = P / a;
        if (fabs((double)a / b - (double)NX / NY) < fabs((double)px / py - (double)NX / NY)) { px = a; py = b; }
    }
    for (int c = 0; c < n; c++) part[c] = (cx[c] * px / NX) + px * (cy[c] * py / NY);
}

/* recursive bisection on the graph (Scotch-like) */
static int *rb_buf;
static void rb(int *cells, int m, int *part, int p0, int np) {
    if (np == 1) { for (int k = 0; k < m; k++) part[cells[k]] = p0; return; }
    /* mark subset with a unique tag */
    static int tagctr = 0;
    int tag = -2 - (++tagctr);
    for (int k = 0; k < m; k++) part[cells[k]] = tag;
    int src = peripheral(cells[0], part, tag, rb_buf);
    /* BFS order of the subset (restart if disconnected) */
    int stamp = ++stamp_ctr, got = bfs(src, part, tag, rb_buf, NULL, seen, stamp);
    for (int k = 0; k < m && got < m; k++)
        if (seen[cells[k]] != stamp) {
            int *tmp = rb_buf + got; got += bfs(cells[k], part, tag, tmp, NULL, seen, stamp);
        }
    memcpy(cells, rb_buf, m * sizeof(int));
    int nl = np / 2, ml = (int)((long)m * nl / np);
    rb(cells, ml, part, p0, nl);
    rb(cells + ml, m - ml, part, p0 + nl, np - nl);
}
static void refine(int *part, int P);
static void dec_scotch(int *part, int P) {
    int *cells = malloc(n * sizeof(int));
    for (int c = 0; c < n; c++) cells[c] = c;
    rb(cells, n, part, 0, P);
    free(cells);
    refine(part, P);
}

/* boundary refinement (greedy, FM-like), shared by scotch* and metis*:
   move a boundary cell to the neighbouring part that holds most of its
   neighbours if that cuts fewer faces (or equal faces but better balance),
   never letting a part exceed 3% above average */
static void refine(int *part, int P) {
    int *size = calloc(P, sizeof(int));
    for (int c = 0; c < n; c++) size[part[c]]++;
    int maxsize = (int)(1.03 * n / P);
    for (int pass = 0; pass < 30; pass++) {
        int moves = 0;
        for (int c = 0; c < n; c++) {
            int a = part[c], cnt[64] = {0}, b = a;
            for (int e = xadj[c]; e < xadj[c + 1]; e++) cnt[part[adj[e]]]++;
            for (int e = xadj[c]; e < xadj[c + 1]; e++) {
                int q = part[adj[e]];
                if (cnt[q] > cnt[b] || (cnt[q] == cnt[b] && size[q] < size[b])) b = q;
            }
            if (b == a || cnt[a] == 0) continue;
            int gain = cnt[b] - cnt[a];
            if ((gain > 0 && size[b] < maxsize) || (gain == 0 && size[b] + 1 < size[a])) {
                part[c] = b; size[a]--; size[b]++; moves++;
            }
        }
        if (!moves) break;
    }
    free(size);
}

/* direct k-way graph growing (METIS-like): P seeds far apart, then the part
   that is currently the smallest always grabs the next cell on its front */
static void dec_metis(int *part, int P) {
    int *dist = malloc(n * sizeof(int)), *queue = malloc(n * sizeof(int)), *seeds = malloc(P * sizeof(int));
    seeds[0] = peripheral(0, NULL, 0, queue);
    for (int s = 1; s < P; s++) {           /* next seed = cell farthest from all seeds */
        for (int c = 0; c < n; c++) dist[c] = -1;
        int head = 0, tail = 0;
        for (int k = 0; k < s; k++) { dist[seeds[k]] = 0; queue[tail++] = seeds[k]; }
        while (head < tail) {
            int c = queue[head++];
            for (int e = xadj[c]; e < xadj[c + 1]; e++)
                if (dist[adj[e]] < 0) { dist[adj[e]] = dist[c] + 1; queue[tail++] = adj[e]; }
        }
        seeds[s] = queue[tail - 1];
    }
    int *size = calloc(P, sizeof(int));
    int **q = malloc(P * sizeof(int *)), *qh = calloc(P, sizeof(int)), *qt = calloc(P, sizeof(int)),
        *qcap = malloc(P * sizeof(int));
    for (int c = 0; c < n; c++) part[c] = -1;
    for (int p = 0; p < P; p++) { qcap[p] = 1 << 16; q[p] = malloc(qcap[p] * sizeof(int)); q[p][qt[p]++] = seeds[p]; }
    for (int assigned = 0; assigned < n;) {
        int p = -1;                            /* smallest part that can still grow */
        for (int k = 0; k < P; k++) if (qh[k] < qt[k] && (p < 0 || size[k] < size[p])) p = k;
        if (p < 0) break;
        int c = q[p][qh[p]++];
        if (part[c] >= 0) continue;
        part[c] = p; size[p]++; assigned++;
        for (int e = xadj[c]; e < xadj[c + 1]; e++)
            if (part[adj[e]] < 0) {
                if (qt[p] == qcap[p]) { qcap[p] *= 2; q[p] = realloc(q[p], qcap[p] * sizeof(int)); }
                q[p][qt[p]++] = adj[e];
            }
    }
    for (int p = 0; p < P; p++) free(q[p]);
    free(q); free(qh); free(qt); free(qcap); free(size); free(dist); free(queue); free(seeds);
    refine(part, P);
}

/* ---------------- metrics ---------------- */
static void report(const char *name, const int *part, int P, double ms) {
    long cut = 0; int *size = calloc(P, sizeof(int));
    char *nbr = calloc((size_t)P * P, 1);
    for (int c = 0; c < n; c++) {
        size[part[c]]++;
        for (int e = xadj[c]; e < xadj[c + 1]; e++)
            if (part[adj[e]] != part[c]) { if (adj[e] > c) cut++; nbr[part[c] * P + part[adj[e]]] = 1; }
    }
    int maxs = 0, maxn = 0, mins = n, idle = 0;
    for (int p = 0; p < P; p++) {
        if (size[p] > maxs) maxs = size[p];
        if (size[p] < mins) mins = size[p];
        if (size[p] == 0) idle++;
        int k = 0; for (int q = 0; q < P; q++) k += nbr[p * P + q];
        if (k > maxn) maxn = k;
    }
    long faces = xadj[n] / 2;
    /* halo per iteration: each cut face -> one double sent each way */
    printf("%-9s %10ld %8.2f%% %10.1f %9.1f%% %9.1f%% %5d %8d %8.0f\n", name, cut, 100.0 * cut / faces,
           cut * 2 * 8 / 1024.0, 100.0 * ((double)maxs * P / n - 1),
           100.0 * ((double)mins * P / n - 1), idle, maxn, ms);
    free(size); free(nbr);
}

static void write_ppm(const char *dir, const char *name, const int *part) {
    if (!dir) return;
    static const unsigned char pal[8][3] = {   /* validated categorical palette */
        {0x39,0x87,0xe5},{0xd9,0x59,0x26},{0x19,0x9e,0x70},{0xc9,0x85,0x00},
        {0xd5,0x51,0x81},{0x0c,0xa3,0x0c},{0x90,0x85,0xe9},{0xe6,0x67,0x67}};
    char path[512]; snprintf(path, sizeof path, "%s/%s.ppm", dir, name);
    FILE *f = fopen(path, "wb"); if (!f) { perror(path); return; }
    fprintf(f, "P6 %d %d 255\n", NX, NY);
    for (int j = NY - 1; j >= 0; j--)
        for (int i = 0; i < NX; i++) {
            int c = cell_at[j][i];
            unsigned char px[3] = {0x08, 0x09, 0x0a};
            if (c >= 0) memcpy(px, pal[part[c] % 8], 3);
            fwrite(px, 1, 3, f);
        }
    fclose(f);
}

/* ---------------- renumbering + sweep ---------------- */
static double sweep_ns(const int *old_of_new) {
    int *new_of_old = malloc(n * sizeof(int));
    for (int k = 0; k < n; k++) new_of_old[old_of_new[k]] = k;
    int *xa = malloc((n + 1) * sizeof(int)), *ad = malloc(xadj[n] * sizeof(int));
    int e2 = 0;
    for (int k = 0; k < n; k++) {
        int c = old_of_new[k]; xa[k] = e2;
        for (int e = xadj[c]; e < xadj[c + 1]; e++) ad[e2++] = new_of_old[adj[e]];
    }
    xa[n] = e2;
    double *x = malloc(n * sizeof(double)), *y = malloc(n * sizeof(double));
    for (int c = 0; c < n; c++) x[new_of_old[c]] = 1.0 + (c & 7);   /* same data, any numbering */
    const int SWEEPS = 50;
    double t0 = now();
    for (int it = 0; it < SWEEPS; it++) {
        for (int k = 0; k < n; k++) {
            double s = 0.0;
            for (int e = xa[k]; e < xa[k + 1]; e++) s += x[ad[e]];
            y[k] = 0.25 * s + 0.01;
        }
        double *t = x; x = y; y = t;
    }
    double ns = (now() - t0) * 1e9 / ((double)SWEEPS * n);
    double chk = 0; for (int k = 0; k < n; k++) chk += x[new_of_old[k]] * (k % 3);
    fprintf(stderr, "  (checksum %.6e)\n", chk);
    free(new_of_old); free(xa); free(ad); free(x); free(y);
    return ns;
}

#ifndef NO_MAIN
int main(int argc, char **argv) {
    int P = (argc > 1) ? atoi(argv[1]) : 8;
    const char *dir = (argc > 2) ? argv[2] : NULL;
    if (P < 2 || P > 64) { fprintf(stderr, "P must be in [2, 64]\n"); return 1; }
    build_mesh();
    seen = calloc(n, sizeof(int)); rb_buf = malloc(n * sizeof(int));
    printf("Mesh: L-shaped domain in a %d x %d box, with cylinder: %d cells, %d faces, P = %d parts\n\n",
           NX, NY, n, xadj[n] / 2, P);

    printf("%-9s %10s %9s %10s %10s %10s %5s %8s %8s\n", "method", "cut faces", "cut %", "halo KB/it",
           "largest", "smallest", "idle", "max nbrs", "time ms");
    int *part = malloc(n * sizeof(int));
    struct { const char *name; void (*fn)(int *, int); } m[] = {
        {"random", dec_random}, {"simple", dec_simple}, {"scotch*", dec_scotch}, {"metis*", dec_metis}};
    for (int k = 0; k < 4; k++) {
        double t0 = now(); m[k].fn(part, P); double ms = (now() - t0) * 1e3;
        report(m[k].name, part, P, ms);
        write_ppm(dir, m[k].name[strlen(m[k].name) - 1] == '*' ?
                  (k == 2 ? "scotch" : "metis") : m[k].name, part);
    }

    /* renumbering: original random order vs BFS (Cuthill-McKee-like) order */
    printf("\nSweep over the mesh (sum of face neighbours, %d cells):\n", n);
    int *order = malloc(n * sizeof(int));
    for (int k = 0; k < n; k++) order[k] = k;
    double t_rand = sweep_ns(order);
    int src = peripheral(0, NULL, 0, order);
    bfs(src, NULL, 0, order, NULL, seen, ++stamp_ctr);
    double t_bfs = sweep_ns(order);
    printf("  random numbering (as generated) : %6.2f ns per cell\n", t_rand);
    printf("  BFS renumbering (renumberMesh)  : %6.2f ns per cell  -> %.1fx faster\n", t_bfs, t_rand / t_bfs);
    return 0;
}
#endif /* NO_MAIN */
