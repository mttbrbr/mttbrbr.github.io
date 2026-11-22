/**
 * Mesh skewness and the skew correction.
 *
 * Periodic 2D mesh of N x N parallelograms on [0,1)^2: the mesh nodes of row j
 * are shifted horizontally by s_j = a h (-1)^j, so cells lean alternately left
 * and right. Vertical faces stay aligned with the line joining the cell
 * centres; horizontal faces stay orthogonal to it (no non-orthogonality), but
 * their centre is displaced by a h from the point where that line crosses the
 * face: pure skewness, of size a (as a fraction of the cell-centre distance).
 *
 * For phi = sin(2 pi x) sin(2 pi y) the program measures the error of
 *   - the Gauss gradient  grad(phi)_P = 1/V sum phi_f S_f
 *   - the convective term div(U phi) with U = (1, 0.5)
 * with face values from linear interpolation, or linear + skew correction
 * phi_f = phi_lin + (x_f - x_i) . grad(phi)_f   (OpenFOAM: skewCorrected linear).
 *
 * Build: gcc -O3 -march=native skew.c -o skew -lm     Run: ./skew outdir
 */
#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static int N; static double h, a;
#define C(i, j) ((((j) + N) % N) * N + (((i) + N) % N))
static double sh(int j) { return a * h * (((j % 2) + 2) % 2 ? -1 : 1); }   /* node shift of row j */
static double phi(double x, double y) { return sin(2 * M_PI * x) * sin(2 * M_PI * y); }
static double phx(double x, double y) { return 2 * M_PI * cos(2 * M_PI * x) * sin(2 * M_PI * y); }
static double phy(double x, double y) { return 2 * M_PI * sin(2 * M_PI * x) * cos(2 * M_PI * y); }

/* cell centres are at ((i+.5) h, (j+.5) h) for every row (the shifts of the two
   edges cancel); face centres: east/west at x = i h, north at x = (i+.5) h + s_{j+1} */
static void gradient(const double *p, const double *gxo, const double *gyo, int corr, double *gx, double *gy) {
    double V = h * h;
    for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
        int c = C(i, j);
        double sj = sh(j), sn = sh(j + 1), sx = 0, sy = 0;
        /* east face: S = (h, 2 s_j)... edge from node (i+1, j) to (i+1, j+1) */
        double Sex = h, Sey = sn - sj;          /* edge (sn - sj, h) -> normal (h, -(sn - sj)) */
        Sey = -(sn - sj);
        double fe = 0.5 * (p[c] + p[C(i + 1, j)]);
        double fw = 0.5 * (p[c] + p[C(i - 1, j)]);
        double fn = 0.5 * (p[c] + p[C(i, j + 1)]);
        double fs = 0.5 * (p[c] + p[C(i, j - 1)]);
        if (corr) {   /* skew correction on north/south faces: displacement (s, 0) */
            fn += sn * 0.5 * (gxo[c] + gxo[C(i, j + 1)]);
            fs += sj * 0.5 * (gxo[c] + gxo[C(i, j - 1)]);
        }
        sx += fe * Sex - fw * Sex; sy += fe * Sey - fw * Sey;
        sy += fn * h - fs * h;
        gx[c] = sx / V; gy[c] = sy / V;
    }
}
static double conv_err(const double *p, const double *gx, const double *gy, int corr) {
    double Ux = 1, Uy = 0.5, V = h * h, e2 = 0, n2 = 0;
    for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
        int c = C(i, j);
        double sj = sh(j), sn = sh(j + 1);
        double Sex = h, Sey = -(sn - sj);
        double fe = 0.5 * (p[c] + p[C(i + 1, j)]), fw = 0.5 * (p[c] + p[C(i - 1, j)]);
        double fn = 0.5 * (p[c] + p[C(i, j + 1)]), fs = 0.5 * (p[c] + p[C(i, j - 1)]);
        if (corr) { fn += sn * 0.5 * (gx[c] + gx[C(i, j + 1)]); fs += sj * 0.5 * (gx[c] + gx[C(i, j - 1)]); }
        double div = ((Ux * Sex + Uy * Sey) * (fe - fw) + Uy * h * (fn - fs)) / V;
        double x = (i + 0.5) * h, y = (j + 0.5) * h, ex = Ux * phx(x, y) + Uy * phy(x, y);
        e2 += (div - ex) * (div - ex); n2 += ex * ex;
    }
    return sqrt(e2 / n2);
}

int main(int argc, char **argv) {
    const char *dir = argc > 1 ? argv[1] : ".";
    char path[512]; snprintf(path, sizeof path, "%s/skewness.csv", dir);
    FILE *f = fopen(path, "w");
    fprintf(f, "skew,N,grad_linear,grad_corrected,conv_linear,conv_corrected\n");
    double as[] = {0, 0.05, 0.1, 0.2, 0.3, 0.4};
    printf("relative L2 error of the Gauss gradient / convective term (linear | skew-corrected)\n");
    for (int q = 0; q < 6; q++) {
        a = as[q];
        printf("skew %.2f:", a);
        for (N = 16; N <= 512; N *= 2) {
            h = 1.0 / N;
            double *p = malloc(N * N * 8), *gx = calloc(N * N, 8), *gy = calloc(N * N, 8), *gx2 = calloc(N * N, 8), *gy2 = calloc(N * N, 8);
            for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) p[C(i, j)] = phi((i + 0.5) * h, (j + 0.5) * h);
            double e[2], cv[2];
            for (int corr = 0; corr <= 1; corr++) {
                gradient(p, gx, gy, 0, gx2, gy2);                       /* linear gradient */
                for (int it = 0; corr && it < 3; it++) {                  /* corrected: iterate using the latest gradient */
                    gradient(p, gx2, gy2, 1, gx, gy);
                    for (int k = 0; k < N * N; k++) { gx2[k] = gx[k]; gy2[k] = gy[k]; }
                }
                double e2 = 0, n2 = 0;
                for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
                    double x = (i + 0.5) * h, y = (j + 0.5) * h, dx = gx2[C(i, j)] - phx(x, y), dy = gy2[C(i, j)] - phy(x, y);
                    e2 += dx * dx + dy * dy; n2 += phx(x, y) * phx(x, y) + phy(x, y) * phy(x, y);
                }
                e[corr] = sqrt(e2 / n2);
                cv[corr] = conv_err(p, gx2, gy2, corr);
            }
            fprintf(f, "%g,%d,%.6e,%.6e,%.6e,%.6e\n", a, N, e[0], e[1], cv[0], cv[1]);
            if (N == 32 || N == 512) printf("  N=%d grad %.1e | %.1e  conv %.1e | %.1e", N, e[0], e[1], cv[0], cv[1]);
            free(p); free(gx); free(gy); free(gx2); free(gy2);
        }
        printf("\n");
    }
    fclose(f);
    return 0;
}
