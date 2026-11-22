/**
 * Halo exchange on a decomposed unstructured mesh, with real MPI messages.
 *
 * Every rank builds the same mesh and the same decomposition (decompose.c,
 * deterministic), keeps only its own cells plus a layer of "halo" copies of
 * the neighbouring cells owned by other ranks, and then iterates:
 *
 *   1. pack the values the neighbours need into send buffers
 *   2. exchange them (MPI_Isend / MPI_Irecv / MPI_Waitall)
 *   3. unpack the received values into the halo cells
 *   4. sweep: every owned cell = 0.25 * sum of its face neighbours + 0.01
 *   5. global reduction (MPI_Allreduce of a residual), like a convergence check
 *
 * Steps 1-3 are timed as "exchange", step 4 as "compute", step 5 as "sync":
 * since the reduction needs every rank, the sync time is mostly time spent
 * waiting for the slowest rank (load imbalance).
 * The global result is checked against a serial run (must be identical).
 *
 * Build: mpicc -O3 -march=native halo.c -o halo -lm
 * Run:   mpirun -np 4 --bind-to core ./halo
 */
#define NO_MAIN
#include "decompose.c"
#include <mpi.h>

#ifndef ITERS
#define ITERS 200
#endif

/* serial reference: same arithmetic, same neighbour order -> bitwise equal */
static double serial(const int *gorder, double *t_iter) {
    double *x = malloc(n * sizeof(double)), *y = malloc(n * sizeof(double));
    for (int c = 0; c < n; c++) x[c] = 1.0 + (c & 7);
    /* renumber in BFS order for a fair (cache friendly) serial baseline */
    int *g2l = malloc(n * sizeof(int));
    for (int k = 0; k < n; k++) g2l[gorder[k]] = k;
    int *xa = malloc((n + 1) * sizeof(int)), *ad = malloc(xadj[n] * sizeof(int));
    double *xl = malloc(n * sizeof(double)), *yl = malloc(n * sizeof(double));
    int e2 = 0;
    for (int k = 0; k < n; k++) {
        int c = gorder[k]; xa[k] = e2; xl[k] = x[c];
        for (int e = xadj[c]; e < xadj[c + 1]; e++) ad[e2++] = g2l[adj[e]];
    }
    xa[n] = e2;
    double t0 = MPI_Wtime();
    for (int it = 0; it < ITERS; it++) {
        for (int k = 0; k < n; k++) {
            double s = 0.0;
            for (int e = xa[k]; e < xa[k + 1]; e++) s += xl[ad[e]];
            yl[k] = 0.25 * s + 0.01;
        }
        double *t = xl; xl = yl; yl = t;
    }
    *t_iter = (MPI_Wtime() - t0) / ITERS;
    double chk = 0;
    for (int k = 0; k < n; k++) chk += xl[k] * (gorder[k] % 3);
    free(x); free(y); free(g2l); free(xa); free(ad); free(xl); free(yl);
    return chk;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, P;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    build_mesh();
    seen = calloc(n, sizeof(int)); rb_buf = malloc(n * sizeof(int));

    /* global BFS order: used to number cells inside every subdomain */
    int *gorder = malloc(n * sizeof(int));
    int src = peripheral(0, NULL, 0, gorder);
    bfs(src, NULL, 0, gorder, NULL, seen, ++stamp_ctr);

    double t_serial = 0, chk_serial = 0;
    if (!rank) chk_serial = serial(gorder, &t_serial);
    MPI_Bcast(&t_serial, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (!rank) {
        printf("Mesh: %d cells, %d faces, %d ranks, %d iterations\n", n, xadj[n] / 2, P, ITERS);
        printf("Serial (1 core): %.3f ms per iteration\n\n", t_serial * 1e3);
        printf("method,cells_max,halo_cells_total,msgs_per_iter,MB_per_iter,"
               "compute_ms_max,exchange_ms_avg,exchange_ms_max,sync_ms_avg,iter_ms,speedup,check\n");
    }

    int *part = malloc(n * sizeof(int));
    struct { const char *name; void (*fn)(int *, int); } m[] = {
        {"random", dec_random}, {"simple", dec_simple}, {"scotch*", dec_scotch}, {"metis*", dec_metis}};

    for (int k = 0; k < 4; k++) {
        m[k].fn(part, P);

        /* ---- local numbering: owned cells first (BFS order), then halo grouped by owner ---- */
        int *g2l = malloc(n * sizeof(int));
        for (int c = 0; c < n; c++) g2l[c] = -1;
        int nown = 0;
        for (int i = 0; i < n; i++) if (part[gorder[i]] == rank) g2l[gorder[i]] = nown++;
        char *need = calloc(n, 1);                 /* halo cells this rank needs */
        for (int c = 0; c < n; c++) if (part[c] == rank)
            for (int e = xadj[c]; e < xadj[c + 1]; e++) if (part[adj[e]] != rank) need[adj[e]] = 1;
        int *roff = calloc(P + 1, sizeof(int)), nloc = nown;
        for (int q = 0; q < P; q++) {
            roff[q] = nloc - nown;
            for (int i = 0; i < n; i++) { int c = gorder[i]; if (need[c] && part[c] == q) g2l[c] = nloc++; }
        }
        roff[P] = nloc - nown;
        /* send lists: my cells adjacent to rank q, in the same (BFS) order q expects */
        int *soff = calloc(P + 1, sizeof(int)), *slist = malloc(4 * (size_t)n * sizeof(int)), ns = 0;
        for (int q = 0; q < P; q++) {
            soff[q] = ns;
            if (q == rank) continue;
            for (int i = 0; i < n; i++) {
                int c = gorder[i];
                if (part[c] != rank) continue;
                for (int e = xadj[c]; e < xadj[c + 1]; e++)
                    if (part[adj[e]] == q) { slist[ns++] = g2l[c]; break; }
            }
        }
        soff[P] = ns;
        /* local CSR, same neighbour order as the global mesh */
        int *xa = malloc((nown + 1) * sizeof(int)), *ad = malloc(4 * (size_t)nown * sizeof(int)), e2 = 0;
        int *gid = malloc(nloc * sizeof(int));
        for (int c = 0; c < n; c++) if (g2l[c] >= 0) gid[g2l[c]] = c;
        for (int l = 0; l < nown; l++) {
            int c = gid[l]; xa[l] = e2;
            for (int e = xadj[c]; e < xadj[c + 1]; e++) ad[e2++] = g2l[adj[e]];
        }
        xa[nown] = e2;
        double *x = malloc(nloc * sizeof(double)), *y = malloc(nloc * sizeof(double));
        for (int l = 0; l < nloc; l++) x[l] = 1.0 + (gid[l] & 7);
        double *sbuf = malloc((ns + 1) * sizeof(double)), *rbuf = malloc((roff[P] + 1) * sizeof(double));
        MPI_Request *req = malloc(2 * P * sizeof(MPI_Request));

        /* ---- iterate ---- */
        double t_exch = 0, t_comp = 0, t_sync = 0;
        MPI_Barrier(MPI_COMM_WORLD);
        double t_start = MPI_Wtime();
        for (int it = 0; it < ITERS; it++) {
            double t0 = MPI_Wtime();
            int nr = 0;
            for (int q = 0; q < P; q++) if (roff[q + 1] > roff[q])
                MPI_Irecv(rbuf + roff[q], roff[q + 1] - roff[q], MPI_DOUBLE, q, it, MPI_COMM_WORLD, &req[nr++]);
            for (int s = 0; s < ns; s++) sbuf[s] = x[slist[s]];          /* pack */
            for (int q = 0; q < P; q++) if (soff[q + 1] > soff[q])
                MPI_Isend(sbuf + soff[q], soff[q + 1] - soff[q], MPI_DOUBLE, q, it, MPI_COMM_WORLD, &req[nr++]);
            MPI_Waitall(nr, req, MPI_STATUSES_IGNORE);
            memcpy(x + nown, rbuf, roff[P] * sizeof(double));             /* unpack */
            double t1 = MPI_Wtime();
            for (int l = 0; l < nown; l++) {
                double s = 0.0;
                for (int e = xa[l]; e < xa[l + 1]; e++) s += x[ad[e]];
                y[l] = 0.25 * s + 0.01;
            }
            double t2 = MPI_Wtime();
            double res = 0, res_all;
            for (int l = 0; l < nown; l += 64) res += fabs(y[l] - x[l]);
            MPI_Allreduce(&res, &res_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            double t3 = MPI_Wtime();
            t_exch += t1 - t0; t_comp += t2 - t1; t_sync += t3 - t2;
            double *t = x; x = y; y = t;
        }
        double t_total = MPI_Wtime() - t_start;

        /* ---- statistics ---- */
        double chk = 0;
        for (int l = 0; l < nown; l++) chk += x[l] * (gid[l] % 3);
        /* sum in BFS order on rank 0 would be needed for bitwise equality of the
           checksum itself; compare with a relative tolerance instead */
        double chk_all, comp_max, exch_max, exch_sum, tot_max, sync_sum;
        long halo = roff[P], halo_all, msgs = 0, msgs_all, cells_max, own = nown;
        for (int q = 0; q < P; q++) msgs += (soff[q + 1] > soff[q]);
        MPI_Reduce(&chk, &chk_all, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_comp, &comp_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_exch, &exch_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_exch, &exch_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_sync, &sync_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&t_total, &tot_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&halo, &halo_all, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&msgs, &msgs_all, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&own, &cells_max, 1, MPI_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
        if (!rank) {
            double it_ms = tot_max / ITERS * 1e3;
            printf("%s,%ld,%ld,%ld,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.2f,%s\n", m[k].name, cells_max, halo_all,
                   msgs_all, halo_all * 8.0 / 1e6, comp_max / ITERS * 1e3, exch_sum / P / ITERS * 1e3,
                   exch_max / ITERS * 1e3, sync_sum / P / ITERS * 1e3, it_ms, t_serial * 1e3 / it_ms,
                   fabs(chk_all - chk_serial) <= 1e-9 * fabs(chk_serial) ? "ok" : "MISMATCH");
        }
        free(g2l); free(need); free(roff); free(soff); free(slist); free(xa); free(ad); free(gid);
        free(x); free(y); free(sbuf); free(rbuf); free(req);
    }
    MPI_Finalize();
    return 0;
}
