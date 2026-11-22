/**
 * Same simulation, different number of "processes", different result.
 *
 * 2D incompressible Navier-Stokes, vorticity-streamfunction, periodic box
 * [0, 2pi)^2, Kolmogorov forcing F sin(k y) in x-momentum, linear drag.
 *   w_t + J(psi, w) = nu lap(w) - alpha w - F k cos(k y),   lap(psi) = -w
 * Arakawa Jacobian, 5-point Laplacians, SSP-RK3 in time. The Poisson equation
 * is solved at every stage with conjugate gradient, like the pressure in a
 * CFD code. The ONLY difference between the runs: every dot product of CG is
 * computed as P partial sums (one per "process", contiguous blocks of rows)
 * added together at the end, which is what MPI_Allreduce does.
 * Mathematically identical; in floating point the last bit can differ.
 *
 * Build: gcc -O3 -march=native kolmogorov.c -o kolmogorov -lm
 * Run:   ./kolmogorov outdir nu tend F
 */
#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 64
#define NR 5                      /* runs: P = 1, 2, 4, 8, and P = 1 with a 1-ulp perturbation */
static const int Ps[NR] = {1, 2, 4, 8, 1};
static double dx, nu = 0.005, alpha = 0.1, F = 0.2; static int kf = 4;

#define I(i, j) ((((j) + N) % N) * N + (((i) + N) % N))

static int P = 1;   /* partitions used by the dot products of the current run */
static double dot(const double *a, const double *b) {
    double part[64] = {0};
    int rows = N / P;
    for (int p = 0; p < P; p++) {               /* each "process" sums its own block of rows */
        double s = 0;
        for (int k = p * rows * N; k < (p + 1) * rows * N; k++) s += a[k] * b[k];
        part[p] = s;
    }
    double s = 0; for (int p = 0; p < P; p++) s += part[p];   /* the reduction */
    return s;
}
static void lap(const double *x, double *y) {
    for (int j = 0; j < N; j++) for (int i = 0; i < N; i++)
        y[I(i, j)] = (x[I(i + 1, j)] + x[I(i - 1, j)] + x[I(i, j + 1)] + x[I(i, j - 1)] - 4 * x[I(i, j)]) / (dx * dx);
}
/* solve -lap(psi) = w with CG (warm start from psi), zero mean */
static int poisson(double *psi, const double *w) {
    static double r[N * N], p[N * N], q[N * N];
    lap(psi, q);
    for (int k = 0; k < N * N; k++) { r[k] = w[k] + q[k]; p[k] = r[k]; }
    double rr = dot(r, r), bb = dot(w, w); int it;
    for (it = 0; it < 2000 && rr > 1e-24 * bb; it++) {
        lap(p, q); for (int k = 0; k < N * N; k++) q[k] = -q[k];
        double al = rr / dot(p, q);
        for (int k = 0; k < N * N; k++) { psi[k] += al * p[k]; r[k] -= al * q[k]; }
        double rn = dot(r, r), be = rn / rr; rr = rn;
        for (int k = 0; k < N * N; k++) p[k] = r[k] + be * p[k];
    }
    return it;
}
/* Arakawa Jacobian J(psi, w) */
static void jacobian(const double *ps, const double *w, double *J) {
    double c = 1.0 / (12 * dx * dx);
    for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
        double j1 = (ps[I(i+1,j)] - ps[I(i-1,j)]) * (w[I(i,j+1)] - w[I(i,j-1)]) - (ps[I(i,j+1)] - ps[I(i,j-1)]) * (w[I(i+1,j)] - w[I(i-1,j)]);
        double j2 = ps[I(i+1,j)] * (w[I(i+1,j+1)] - w[I(i+1,j-1)]) - ps[I(i-1,j)] * (w[I(i-1,j+1)] - w[I(i-1,j-1)])
                  - ps[I(i,j+1)] * (w[I(i+1,j+1)] - w[I(i-1,j+1)]) + ps[I(i,j-1)] * (w[I(i+1,j-1)] - w[I(i-1,j-1)]);
        double j3 = w[I(i,j+1)] * (ps[I(i+1,j+1)] - ps[I(i-1,j+1)]) - w[I(i,j-1)] * (ps[I(i+1,j-1)] - ps[I(i-1,j-1)])
                  - w[I(i+1,j)] * (ps[I(i+1,j+1)] - ps[I(i+1,j-1)]) + w[I(i-1,j)] * (ps[I(i-1,j+1)] - ps[I(i-1,j-1)]);
        J[I(i, j)] = (j1 + j2 + j3) * c;
    }
}
static double *forcing;
static void rhs(double *w, double *psi, double *out) {
    static double J[N * N], L[N * N];
    poisson(psi, w);
    jacobian(psi, w, J); lap(w, L);
    for (int k = 0; k < N * N; k++) out[k] = -J[k] + nu * L[k] - alpha * w[k] + forcing[k];
}

int main(int argc, char **argv) {
    const char *dir = argc > 1 ? argv[1] : ".";
    nu = argc > 2 ? atof(argv[2]) : 0.005;
    F = argc > 4 ? atof(argv[4]) : 0.2;
    double tend = argc > 3 ? atof(argv[3]) : 200, dt = 0.01;
    dx = 2 * M_PI / N;
    forcing = malloc(N * N * 8);
    for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) forcing[I(i, j)] = -F * kf * cos(kf * j * dx);
    double *w[NR], *psi[NR], *w1 = malloc(N * N * 8), *w2 = malloc(N * N * 8), *k1 = malloc(N * N * 8);
    /* identical initial condition: laminar profile + a small deterministic perturbation */
    for (int r = 0; r < NR; r++) {
        w[r] = malloc(N * N * 8); psi[r] = calloc(N * N, 8);
        for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
            double x = i * dx, y = j * dx;
            w[r][I(i, j)] = -0.5 * kf * cos(kf * y) + 0.05 * sin(x + 0.3) * cos(2 * y + 0.7) + 0.03 * cos(3 * x + 1.1) * sin(y);
        }
    }
    w[4][I(10, 10)] = nextafter(w[4][I(10, 10)], 1e9);      /* run 5: one cell changed by 1 ulp */

    char path[512]; snprintf(path, sizeof path, "%s/divergence_nu%g_F%g.csv", dir, nu, F);
    FILE *f = fopen(path, "w");
    fprintf(f, "t,diff_P2,diff_P4,diff_P8,diff_ulp,E_P1,E_P2,E_P4,E_P8,Z_P1,Z_P2,Z_P4,Z_P8\n");
    long steps = (long)(tend / dt + 0.5);
    double Esum[NR] = {0}, Zsum[NR] = {0}; long nstat = 0;
    for (long s = 0; s <= steps; s++) {
        if (s % 25 == 0) {
            double n1 = sqrt(dot(w[0], w[0])), E[NR], Z[NR];
            for (int r = 0; r < NR; r++) {
                static double uu[N * N];
                /* energy = 1/2 <|grad psi|^2> = 1/2 <psi w>, enstrophy = 1/2 <w^2> */
                P = 1; poisson(psi[r], w[r]);
                double e = 0, z = 0; for (int k = 0; k < N * N; k++) { e += psi[r][k] * w[r][k]; z += w[r][k] * w[r][k]; }
                E[r] = 0.5 * e / (N * N); Z[r] = 0.5 * z / (N * N); (void)uu;
            }
            double d[NR];
            for (int r = 1; r < NR; r++) { double s2 = 0; for (int k = 0; k < N * N; k++) s2 += (w[r][k] - w[0][k]) * (w[r][k] - w[0][k]); d[r] = sqrt(s2) / n1; }
            fprintf(f, "%.3f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n", s * dt, d[1], d[2], d[3], d[4],
                    E[0], E[1], E[2], E[3], Z[0], Z[1], Z[2], Z[3]);
            if (s * dt >= tend / 2) { for (int r = 0; r < NR; r++) { Esum[r] += E[r]; Zsum[r] += Z[r]; } nstat++; }
        }
        if (s == steps) break;
        for (int r = 0; r < NR; r++) {                       /* SSP-RK3 */
            P = Ps[r];
            rhs(w[r], psi[r], k1);
            for (int k = 0; k < N * N; k++) w1[k] = w[r][k] + dt * k1[k];
            rhs(w1, psi[r], k1);
            for (int k = 0; k < N * N; k++) w2[k] = 0.75 * w[r][k] + 0.25 * (w1[k] + dt * k1[k]);
            rhs(w2, psi[r], k1);
            for (int k = 0; k < N * N; k++) w[r][k] = w[r][k] / 3 + 2.0 / 3 * (w2[k] + dt * k1[k]);
        }
    }
    fclose(f);
    for (int r = 0; r < NR; r++) {                   /* final vorticity fields */
        snprintf(path, sizeof path, "%s/vorticity_nu%g_F%g_run%d.csv", dir, nu, F, r);
        FILE *g = fopen(path, "w");
        for (int j = 0; j < N; j++) { for (int i = 0; i < N; i++) fprintf(g, "%s%.6e", i ? "," : "", w[r][I(i, j)]); fprintf(g, "\n"); }
        fclose(g);
    }
    printf("nu = %g, F = %g: time-averaged statistics over t in [%g, %g]\n", nu, F, tend / 2, tend);
    const char *lab[NR] = {"P = 1", "P = 2", "P = 4", "P = 8", "P = 1, 1 ulp"};
    for (int r = 0; r < NR; r++) printf("  %-14s energy %.6f  enstrophy %.6f\n", lab[r], Esum[r] / nstat, Zsum[r] / nstat);
    return 0;
}
