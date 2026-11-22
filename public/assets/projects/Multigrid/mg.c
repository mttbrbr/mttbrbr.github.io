/**
 * Why multigrid: the same Poisson equation (the pressure equation of an
 * incompressible solver) solved with four methods on grids of growing size.
 *
 *   -lap(u) = f on the unit square, u = 0 on the boundary,
 *   exact u = x(1-x) y(1-y^2)  ->  f = 2 y (1-y^2) + 6 x y (1-x)
 *   (not an eigenvector of the matrix, so no method gets lucky)
 *   5-point finite differences on N x N interior points, h = 1/(N+1), N = 2^k - 1.
 *
 *   GS     : Gauss-Seidel (lexicographic) iterations
 *   PCG    : conjugate gradient, Jacobi (diagonal) preconditioner
 *   MG     : geometric multigrid V-cycles (2+2 Gauss-Seidel sweeps,
 *            full-weighting restriction, bilinear prolongation)
 *   MG-PCG : conjugate gradient preconditioned by one symmetric V-cycle
 *            (the role GAMG plays as a preconditioner in OpenFOAM)
 *
 * All stop when the residual has dropped by 1e8.
 * Build: gcc -O3 -march=native mg.c -o mg -lm      Run: ./mg outdir
 */
#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define TOL 1e-8
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static FILE *hist;
static char tag[64];

/* grids include a ghost ring of zeros: arrays are (N+2)^2, index (i,j) -> i + (N+2) j */
#define ID(i, j, n) ((i) + ((n) + 2) * (j))

static double *alloc(int n) { return calloc((size_t)(n + 2) * (n + 2), sizeof(double)); }

static void residual(int n, const double *u, const double *f, double *r) {   /* r = f - A u, A = -lap * h^2 */
    for (int j = 1; j <= n; j++) for (int i = 1; i <= n; i++)
        r[ID(i, j, n)] = f[ID(i, j, n)] - (4 * u[ID(i, j, n)] - u[ID(i - 1, j, n)] - u[ID(i + 1, j, n)]
                                          - u[ID(i, j - 1, n)] - u[ID(i, j + 1, n)]);
}
static double norm(int n, const double *r) {
    double s = 0; for (int j = 1; j <= n; j++) for (int i = 1; i <= n; i++) s += r[ID(i, j, n)] * r[ID(i, j, n)];
    return sqrt(s);
}
static void gs_sweep(int n, double *u, const double *f, int reverse) {
    if (!reverse) {
        for (int j = 1; j <= n; j++) for (int i = 1; i <= n; i++)
            u[ID(i, j, n)] = 0.25 * (f[ID(i, j, n)] + u[ID(i - 1, j, n)] + u[ID(i + 1, j, n)] + u[ID(i, j - 1, n)] + u[ID(i, j + 1, n)]);
    } else {
        for (int j = n; j >= 1; j--) for (int i = n; i >= 1; i--)
            u[ID(i, j, n)] = 0.25 * (f[ID(i, j, n)] + u[ID(i - 1, j, n)] + u[ID(i + 1, j, n)] + u[ID(i, j - 1, n)] + u[ID(i, j + 1, n)]);
    }
}

/* one V-cycle on A u = f (scaled by h^2: the coarse operator is the same stencil
   with f scaled by 4) */
static void vcycle(int n, double *u, const double *f, int symmetric, double weight) {
    if (n <= 3) { for (int s = 0; s < 50; s++) gs_sweep(n, u, f, 0); return; }
    for (int s = 0; s < 2; s++) gs_sweep(n, u, f, 0);
    int nc = (n - 1) / 2;
    double *r = alloc(n), *fc = alloc(nc), *uc = alloc(nc);
    residual(n, u, f, r);
    for (int J = 1; J <= nc; J++) for (int I = 1; I <= nc; I++) {         /* full weighting */
        int i = 2 * I, j = 2 * J;
        double v = 4 * r[ID(i, j, n)] + 2 * (r[ID(i - 1, j, n)] + r[ID(i + 1, j, n)] + r[ID(i, j - 1, n)] + r[ID(i, j + 1, n)])
                 + r[ID(i - 1, j - 1, n)] + r[ID(i + 1, j - 1, n)] + r[ID(i - 1, j + 1, n)] + r[ID(i + 1, j + 1, n)];
        fc[ID(I, J, nc)] = 4 * v / 16;      /* x4: coarse h^2 is 4 times the fine one */
    }
    vcycle(nc, uc, fc, symmetric, weight / 4);
    for (int J = 0; J <= nc; J++) for (int I = 0; I <= nc; I++) {         /* bilinear prolongation */
        int i = 2 * I, j = 2 * J;
        double a = uc[ID(I, J, nc)], b = uc[ID(I + 1, J, nc)], c = uc[ID(I, J + 1, nc)], d = uc[ID(I + 1, J + 1, nc)];
        if (i >= 1 && j >= 1 && i <= n && j <= n) u[ID(i, j, n)] += a;
        if (i + 1 <= n && j >= 1 && j <= n) u[ID(i + 1, j, n)] += 0.5 * (a + b);
        if (j + 1 <= n && i >= 1 && i <= n) u[ID(i, j + 1, n)] += 0.5 * (a + c);
        if (i + 1 <= n && j + 1 <= n) u[ID(i + 1, j + 1, n)] += 0.25 * (a + b + c + d);
    }
    for (int s = 0; s < 2; s++) gs_sweep(n, u, f, symmetric);              /* reverse order: symmetric cycle */
    free(r); free(fc); free(uc);
}

typedef struct { long iters; double seconds, err; int ok; } Res;

static double error_inf(int n, const double *u) {
    double h = 1.0 / (n + 1), e = 0;
    for (int j = 1; j <= n; j++) for (int i = 1; i <= n; i++) {
        double x = i * h, y = j * h, ex = x * (1 - x) * y * (1 - y * y), d = fabs(u[ID(i, j, n)] - ex);
        if (d > e) e = d;
    }
    return e;
}
static void setup_f(int n, double *f) {
    double h = 1.0 / (n + 1);
    for (int j = 1; j <= n; j++) for (int i = 1; i <= n; i++)
        { double x = i * h, y = j * h; f[ID(i, j, n)] = h * h * (2 * y * (1 - y * y) + 6 * x * y * (1 - x)); }
}

static Res run_gs(int n, long maxit) {
    double *u = alloc(n), *f = alloc(n), *r = alloc(n); setup_f(n, f);
    residual(n, u, f, r); double r0 = norm(n, r), t0 = now(); Res R = {0};
    for (long it = 1; it <= maxit; it++) {
        gs_sweep(n, u, f, 0);
        if (it % 10 == 0 || it < 10) {
            residual(n, u, f, r); double rr = norm(n, r) / r0;
            if (hist) fprintf(hist, "%s,%ld,%.6e\n", tag, it, rr);
            if (rr < TOL) { R.iters = it; R.ok = 1; break; }
        }
        R.iters = it;
    }
    R.seconds = now() - t0; R.err = error_inf(n, u);
    free(u); free(f); free(r); return R;
}
static Res run_mg(int n) {
    double *u = alloc(n), *f = alloc(n), *r = alloc(n); setup_f(n, f);
    residual(n, u, f, r); double r0 = norm(n, r), t0 = now(); Res R = {0};
    for (int it = 1; it <= 200; it++) {
        vcycle(n, u, f, 0, 1);
        residual(n, u, f, r); double rr = norm(n, r) / r0;
        if (hist) fprintf(hist, "%s,%d,%.6e\n", tag, it, rr);
        R.iters = it;
        if (rr < TOL) { R.ok = 1; break; }
    }
    R.seconds = now() - t0; R.err = error_inf(n, u);
    free(u); free(f); free(r); return R;
}
/* PCG with preconditioner: 0 = Jacobi, 1 = one symmetric V-cycle */
static Res run_pcg(int n, int mgprec) {
    double *x = alloc(n), *b = alloc(n), *r = alloc(n), *z = alloc(n), *p = alloc(n), *q = alloc(n);
    setup_f(n, b);
    memcpy(r, b, sizeof(double) * (n + 2) * (n + 2));
    double r0 = norm(n, r), t0 = now(); Res R = {0};
    #define PREC(zz, rr) do { if (mgprec) { memset(zz, 0, sizeof(double) * (n + 2) * (n + 2)); vcycle(n, zz, rr, 1, 1); } \
                              else for (int k = 0; k < (n + 2) * (n + 2); k++) zz[k] = rr[k] / 4; } while (0)
    PREC(z, r); memcpy(p, z, sizeof(double) * (n + 2) * (n + 2));
    double rz = 0; for (int j = 1; j <= n; j++) for (int i = 1; i <= n; i++) rz += r[ID(i, j, n)] * z[ID(i, j, n)];
    for (int it = 1; it <= 100000; it++) {
        for (int j = 1; j <= n; j++) for (int i = 1; i <= n; i++)
            q[ID(i, j, n)] = 4 * p[ID(i, j, n)] - p[ID(i - 1, j, n)] - p[ID(i + 1, j, n)] - p[ID(i, j - 1, n)] - p[ID(i, j + 1, n)];
        double pq = 0; for (int j = 1; j <= n; j++) for (int i = 1; i <= n; i++) pq += p[ID(i, j, n)] * q[ID(i, j, n)];
        double al = rz / pq;
        for (int j = 1; j <= n; j++) for (int i = 1; i <= n; i++) { x[ID(i, j, n)] += al * p[ID(i, j, n)]; r[ID(i, j, n)] -= al * q[ID(i, j, n)]; }
        double rr = norm(n, r) / r0;
        if (hist) fprintf(hist, "%s,%d,%.6e\n", tag, it, rr);
        R.iters = it;
        if (rr < TOL) { R.ok = 1; break; }
        PREC(z, r);
        double rz2 = 0; for (int j = 1; j <= n; j++) for (int i = 1; i <= n; i++) rz2 += r[ID(i, j, n)] * z[ID(i, j, n)];
        double be = rz2 / rz; rz = rz2;
        for (int j = 1; j <= n; j++) for (int i = 1; i <= n; i++) p[ID(i, j, n)] = z[ID(i, j, n)] + be * p[ID(i, j, n)];
    }
    R.seconds = now() - t0; R.err = error_inf(n, x);
    free(x); free(b); free(r); free(z); free(p); free(q); return R;
}

int main(int argc, char **argv) {
    const char *dir = argc > 1 ? argv[1] : ".";
    char path[512];
    snprintf(path, sizeof path, "%s/history.csv", dir); hist = fopen(path, "w"); fprintf(hist, "case,iter,residual\n");
    snprintf(path, sizeof path, "%s/scaling.csv", dir); FILE *sc = fopen(path, "w");
    fprintf(sc, "N,unknowns,method,iters,seconds,max_error\n");
    printf("%6s %9s | %-26s | %-26s | %-26s | %-26s\n", "N", "unknowns", "Gauss-Seidel", "PCG (Jacobi)", "Multigrid V-cycle", "PCG + multigrid");
    printf("%6s %9s |", "", ""); for (int k = 0; k < 4; k++) printf(" %8s %9s %6s |", "iters", "time (s)", ""); printf("\n");
    int Ns[] = {15, 31, 63, 127, 255, 511, 1023};
    for (int k = 0; k < 7; k++) {
        int n = Ns[k]; FILE *h0 = hist;
        Res g = {0}, c, m, mc;
        if (n != 127) hist = NULL;
        snprintf(tag, sizeof tag, "gs"); if (n <= 255) g = run_gs(n, 2000000);
        snprintf(tag, sizeof tag, "pcg"); c = run_pcg(n, 0);
        snprintf(tag, sizeof tag, "mg"); m = run_mg(n);
        snprintf(tag, sizeof tag, "mgpcg"); mc = run_pcg(n, 1);
        hist = h0;
        printf("%6d %9d |", n, n * n);
        if (n <= 255) printf(" %8ld %9.3f %6s |", g.iters, g.seconds, ""); else printf(" %8s %9s %6s |", "-", "-", "");
        printf(" %8ld %9.3f %6s | %8ld %9.3f %6s | %8ld %9.3f %6s |\n", c.iters, c.seconds, "", m.iters, m.seconds, "", mc.iters, mc.seconds, "");
        if (n <= 255) fprintf(sc, "%d,%d,gs,%ld,%.5f,%.3e\n", n, n * n, g.iters, g.seconds, g.err);
        fprintf(sc, "%d,%d,pcg,%ld,%.5f,%.3e\n%d,%d,mg,%ld,%.5f,%.3e\n%d,%d,mgpcg,%ld,%.5f,%.3e\n",
                n, n * n, c.iters, c.seconds, c.err, n, n * n, m.iters, m.seconds, m.err, n, n * n, mc.iters, mc.seconds, mc.err);
        fflush(stdout);
    }
    fclose(sc); fclose(hist);
    return 0;
}
