/**
 * Why a density-based solver struggles at low Mach number (1D Euler equations).
 *
 *  A. Stiffness and dissipation. A density blob is carried by a uniform flow
 *     (u, p constant) once across a periodic domain, with an explicit
 *     density-based finite-volume scheme (Rusanov / local Lax-Friedrichs flux).
 *     The time step is set by the fastest wave, |u| + c; the blob moves at u.
 *     Reference: what an incompressible (pressure-based) solver reduces to for
 *     this flow: the blob advected at u with upwind, time step set by u only.
 *  B. Measured numerical viscosity of the density-based scheme vs Mach number,
 *     from the spreading of the blob over a short time.
 *  C. Round-off: recover the pressure fluctuation (~ M^2 p0) from the conserved
 *     variables (rho, rho u, E) in single and double precision.
 *  D. Eigenvalues (wave speeds) of the 1D Euler equations with and without a
 *     Turkel-type low-Mach preconditioner, and their ratio (stiffness).
 *
 * Build: gcc -O3 -march=native lowmach.c -o lowmach -lm
 * Run:   ./lowmach outdir
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define GAMMA 1.4
#define N 200
#define CFL 0.9
#define AMP 0.1          /* density blob amplitude */
#define SIG 0.05         /* blob width */

static const char *dir = ".";
static FILE *out(const char *name) {
    char p[512]; snprintf(p, sizeof p, "%s/%s", dir, name);
    FILE *f = fopen(p, "w"); if (!f) { perror(p); exit(1); } return f;
}
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

/* sound speed c0 = 1 at rho = 1: p0 = 1/gamma. Velocity u = M. */
static double blob(double x) { return 1.0 + AMP * exp(-0.5 * pow((x - 0.5) / SIG, 2)); }

/* ---- density-based: explicit Rusanov on (rho, rho u, E) ---- */
typedef struct { double r, m, E; } Q;
static void flux(Q q, double *f, double *s) {
    double u = q.m / q.r, p = (GAMMA - 1) * (q.E - 0.5 * q.m * u);
    f[0] = q.m; f[1] = q.m * u + p; f[2] = (q.E + p) * u;
    *s = fabs(u) + sqrt(GAMMA * p / q.r);
}
/* advance until time T; returns number of steps; rho written to out */
static long density_based(double M, double T, double *rho, double *wall) {
    double dx = 1.0 / N, p0 = 1.0 / GAMMA;
    Q *q = malloc(N * sizeof(Q)), *qn = malloc(N * sizeof(Q));
    double (*F)[3] = malloc(N * sizeof *F);
    for (int i = 0; i < N; i++) {
        double r = blob((i + 0.5) * dx);
        q[i] = (Q){r, r * M, p0 / (GAMMA - 1) + 0.5 * r * M * M};
    }
    double t = 0; long steps = 0, t0c = 0; (void)t0c;
    double w0 = now();
    while (t < T - 1e-12) {
        double smax = 0;
        for (int i = 0; i < N; i++) { double f[3], s; flux(q[i], f, &s); if (s > smax) smax = s; }
        double dt = CFL * dx / smax; if (t + dt > T) dt = T - t;
        for (int i = 0; i < N; i++) {                 /* face i+1/2 */
            Q L = q[i], R = q[(i + 1) % N];
            double fl[3], fr[3], sl, sr;
            flux(L, fl, &sl); flux(R, fr, &sr);
            double s = sl > sr ? sl : sr;
            F[i][0] = 0.5 * (fl[0] + fr[0]) - 0.5 * s * (R.r - L.r);
            F[i][1] = 0.5 * (fl[1] + fr[1]) - 0.5 * s * (R.m - L.m);
            F[i][2] = 0.5 * (fl[2] + fr[2]) - 0.5 * s * (R.E - L.E);
        }
        for (int i = 0; i < N; i++) {
            int im = (i - 1 + N) % N;
            qn[i].r = q[i].r - dt / dx * (F[i][0] - F[im][0]);
            qn[i].m = q[i].m - dt / dx * (F[i][1] - F[im][1]);
            qn[i].E = q[i].E - dt / dx * (F[i][2] - F[im][2]);
        }
        memcpy(q, qn, N * sizeof(Q)); t += dt; steps++;
    }
    *wall = now() - w0;
    for (int i = 0; i < N; i++) rho[i] = q[i].r;
    free(q); free(qn); free(F);
    return steps;
}

/* ---- incompressible reference: blob advected at u, upwind, dt = CFL dx / u ---- */
static long incompressible(double M, double T, double *rho, double *wall) {
    double dx = 1.0 / N, *r = malloc(N * sizeof(double)), *rn = malloc(N * sizeof(double));
    for (int i = 0; i < N; i++) r[i] = blob((i + 0.5) * dx);
    double t = 0; long steps = 0, w0 = 0; double ws = now(); (void)w0;
    while (t < T - 1e-12) {
        double dt = CFL * dx / M; if (t + dt > T) dt = T - t;
        double C = M * dt / dx;
        for (int i = 0; i < N; i++) rn[i] = r[i] - C * (r[i] - r[(i - 1 + N) % N]);
        memcpy(r, rn, N * sizeof(double)); t += dt; steps++;
    }
    *wall = now() - ws;
    memcpy(rho, r, N * sizeof(double)); free(r); free(rn);
    return steps;
}

static double variance(const double *rho) {
    double m0 = 0, m1 = 0, m2 = 0, dx = 1.0 / N;
    for (int i = 0; i < N; i++) { double w = rho[i] - 1.0, x = (i + 0.5) * dx; m0 += w; m1 += w * x; m2 += w * x * x; }
    double mu = m1 / m0; return m2 / m0 - mu * mu;
}
/* centre the periodic profile on the original blob position (shift by u T) */
static double peak(const double *rho) { double mx = 0; for (int i = 0; i < N; i++) if (rho[i] - 1 > mx) mx = rho[i] - 1; return mx / AMP; }

int main(int argc, char **argv) {
    if (argc > 1) dir = argv[1];
    double Ms[] = {0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001};
    int nM = sizeof Ms / sizeof *Ms;
    double *rho = malloc(N * sizeof(double));

    /* A. one crossing of the domain: T = 1 / u */
    FILE *f = out("crossing.csv");
    fprintf(f, "M,solver,steps,wall_ms,peak_left\n");
    printf("A. blob carried once across the domain (T = L/u), %d cells, CFL %.1f\n", N, CFL);
    printf("   %7s | %-28s | %-28s\n", "", "density-based (Rusanov)", "incompressible (upwind)");
    printf("   %7s | %9s %8s %9s | %9s %8s %9s\n", "Mach", "steps", "ms", "peak", "steps", "ms", "peak");
    FILE *pf = out("profiles.csv");
    fprintf(pf, "M,x,initial,density_based,incompressible\n");
    for (int k = 0; k < nM; k++) {
        double M = Ms[k], T = 1.0 / M, w1, w2;
        long s1 = density_based(M, T, rho, &w1); double p1 = peak(rho);
        double *r1 = malloc(N * sizeof(double)); memcpy(r1, rho, N * sizeof(double));
        long s2 = incompressible(M, T, rho, &w2); double p2 = peak(rho);
        printf("   %7.3f | %9ld %8.2f %8.1f%% | %9ld %8.2f %8.1f%%\n", M, s1, w1 * 1e3, 100 * p1, s2, w2 * 1e3, 100 * p2);
        fprintf(f, "%g,density,%ld,%.3f,%.4f\n%g,incompressible,%ld,%.3f,%.4f\n", M, s1, w1 * 1e3, p1, M, s2, w2 * 1e3, p2);
        if (M == 0.5 || M == 0.1 || M == 0.01)
            for (int i = 0; i < N; i++) {
                double x = (i + 0.5) / N;
                fprintf(pf, "%g,%.4f,%.6f,%.6f,%.6f\n", M, x, blob(x), r1[i], rho[i]);
            }
        free(r1);
    }
    fclose(f); fclose(pf);

    /* B. numerical viscosity from short-time spreading */
    f = out("viscosity.csv");
    fprintf(f, "M,nu_measured,nu_rusanov_pred,nu_upwind_incomp\n");
    printf("\nB. numerical viscosity of the density-based scheme (dx = 1/%d)\n", N);
    printf("   %7s %14s %14s %14s %10s\n", "Mach", "measured", "(|u|+c)dx/2", "u dx/2", "ratio");
    for (int k = 0; k < nM; k++) {
        double M = Ms[k], dx = 1.0 / N, w;
        double s = M + 1.0, nu_pred = 0.5 * dx * (s - CFL * M * M / s);
        double v0, *r0 = malloc(N * sizeof(double));
        for (int i = 0; i < N; i++) r0[i] = blob((i + 0.5) * dx);
        v0 = variance(r0); free(r0);
        double T = 0.5 * SIG * SIG / (2 * nu_pred);        /* variance grows by ~50% */
        density_based(M, T, rho, &w);
        double nu = (variance(rho) - v0) / (2 * T);
        printf("   %7.3f %14.4e %14.4e %14.4e %9.0fx\n", M, nu, nu_pred, 0.5 * M * dx, nu / (0.5 * M * dx));
        fprintf(f, "%g,%.6e,%.6e,%.6e\n", M, nu, nu_pred, 0.5 * M * dx);
    }
    fclose(f);

    /* C. round-off in the pressure */
    f = out("roundoff.csv");
    fprintf(f, "M,rel_err_float,rel_err_double\n");
    printf("\nC. pressure fluctuation dp = M^2 p0 recovered from (rho, rho u, E)\n");
    printf("   %8s %16s %16s\n", "Mach", "error float", "error double");
    for (double M = 1e-1; M > 1e-9; M /= sqrt(10.0)) {
        double p0 = 1.0 / GAMMA, dp = M * M * p0, r = 1.0, u = M;
        /* exact conserved variables, rounded to each precision */
        float rf = (float)r, mf = (float)(r * u), Ef = (float)((p0 + dp) / (GAMMA - 1) + 0.5 * r * u * u);
        float pf_ = (float)(GAMMA - 1) * (Ef - 0.5f * mf * mf / rf);
        double rd = r, md = r * u, Ed = (p0 + dp) / (GAMMA - 1) + 0.5 * r * u * u;
        double pd = (GAMMA - 1) * (Ed - 0.5 * md * md / rd);
        double ef = fabs(((double)pf_ - p0) - dp) / dp, ed = fabs((pd - p0) - dp) / dp;
        printf("   %8.1e %15.2e%% %15.2e%%\n", M, 100 * ef, 100 * ed);
        fprintf(f, "%.3e,%.6e,%.6e\n", M, ef, ed);
    }
    fclose(f);

    /* D. wave speeds with and without preconditioning */
    f = out("eigen.csv");
    fprintf(f, "M,lam_min,lam_max,kappa,lamp_min,lamp_max,kappa_precond\n");
    printf("\nD. wave speeds of the 1D Euler equations (c = 1, u = M)\n");
    printf("   %7s | %9s %9s %10s | %9s %9s %8s\n", "Mach", "slowest", "fastest", "ratio", "slowest", "fastest", "ratio");
    printf("   %7s | %-30s | %-28s\n", "", "plain", "Turkel preconditioned");
    for (int k = 0; k < nM; k++) {
        double M = Ms[k], u = M, c = 1.0;
        /* plain: u, u + c, u - c */
        double l[3] = {u, u + c, u - c}, mn = 1e30, mx = 0;
        for (int i = 0; i < 3; i++) { double a = fabs(l[i]); if (a < mn) mn = a; if (a > mx) mx = a; }
        /* preconditioned (p,u) subsystem: P = diag(eps, 1), A = [[u, rho c^2], [1/rho, u]]
           eigenvalues of P A: 1/2 [ (1+eps) u +- sqrt((1-eps)^2 u^2 + 4 eps c^2) ], eps = M^2 */
        double eps = M * M < 1 ? M * M : 1;
        double sq = sqrt((1 - eps) * (1 - eps) * u * u + 4 * eps * c * c);
        double lp[3] = {u, 0.5 * ((1 + eps) * u + sq), 0.5 * ((1 + eps) * u - sq)}, pmn = 1e30, pmx = 0;
        for (int i = 0; i < 3; i++) { double a = fabs(lp[i]); if (a < pmn) pmn = a; if (a > pmx) pmx = a; }
        printf("   %7.3f | %9.4f %9.4f %10.1f | %9.4f %9.4f %8.2f\n", M, mn, mx, mx / mn, pmn, pmx, pmx / pmn);
        fprintf(f, "%g,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e\n", M, mn, mx, mx / mn, pmn, pmx, pmx / pmn);
    }
    fclose(f);
    free(rho);
    return 0;
}
