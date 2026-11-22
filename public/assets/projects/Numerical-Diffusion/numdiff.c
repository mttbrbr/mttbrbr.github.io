/**
 * Numerical diffusion, measured. Three experiments, all in well under a second.
 *
 *  A. 1D linear advection u_t + a u_x = 0 on a periodic domain, one revolution,
 *     with three schemes written as one flux-limited family (Sweby):
 *        F_{i+1/2} = a u_i + a/2 (1 - C) phi(r_i) (u_{i+1} - u_i)
 *        phi = 0         -> first-order upwind
 *        phi = 1         -> Lax-Wendroff (second order)
 *        phi = van Leer  -> TVD, second order where smooth, no new extrema
 *     Writes profiles (profiles.csv) and L1 errors vs grid (convergence.csv).
 *
 *  B. Upwind on a Gaussian: measure how fast its variance grows and turn it into
 *     a viscosity, nu = d(sigma^2)/dt / 2. Compare with the modified equation,
 *     nu_num = a dx (1 - C) / 2  (viscosity.csv).
 *
 *  C. 2D steady advection of a step at angle theta on an N x N grid with upwind
 *     (false diffusion). Writes the field (field2d_<theta>.csv) and the width of
 *     the smeared step vs distance (falsediff.csv).
 *
 * Build: gcc -O3 -march=native numdiff.c -o numdiff -lm
 * Run:   ./numdiff outdir
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static const char *dir = ".";
static FILE *out(const char *name) {
    char p[512]; snprintf(p, sizeof p, "%s/%s", dir, name);
    FILE *f = fopen(p, "w"); if (!f) { perror(p); exit(1); } return f;
}

/* ---------------- 1D ---------------- */
enum { UPWIND, LAXWENDROFF, VANLEER };
static const char *sname[] = {"upwind", "lax-wendroff", "van-leer"};

static double limiter(int s, double r) {
    if (s == UPWIND) return 0.0;
    if (s == LAXWENDROFF) return 1.0;
    return (r + fabs(r)) / (1.0 + fabs(r));          /* van Leer */
}

/* initial conditions on [0,1): 1 = square, 2 = gaussian, 3 = both */
#define SIG0 0.04
static double ic(int which, double x) {
    double v = 0;
    if (which & 1) v += (x > 0.1 && x < 0.3) ? 1.0 : 0.0;
    if (which & 2) v += exp(-0.5 * pow((x - 0.65) / SIG0, 2));
    return v;
}
/* cell average of the initial condition (fine sampling) */
static double ic_avg(int which, double x0, double dx) {
    double s = 0; const int K = 64;
    for (int k = 0; k < K; k++) s += ic(which, x0 + (k + 0.5) * dx / K);
    return s / K;
}

/* advect for time T (a = 1); returns the solution in u (length n) */
static void advect(int scheme, int n, double C, double T, double *u) {
    double dx = 1.0 / n, dt = C * dx;
    int steps = (int)llround(T / dt);
    double *F = malloc(n * sizeof(double)), *un = malloc(n * sizeof(double));
    for (int s = 0; s < steps; s++) {
        for (int i = 0; i < n; i++) {                 /* flux at face i+1/2 */
            int im = (i - 1 + n) % n, ip = (i + 1) % n;
            double du = u[ip] - u[i], r = (fabs(du) > 1e-14) ? (u[i] - u[im]) / du : 0.0;
            F[i] = u[i] + 0.5 * (1.0 - C) * limiter(scheme, r) * du;
        }
        for (int i = 0; i < n; i++) un[i] = u[i] - C * (F[i] - F[(i - 1 + n) % n]);
        memcpy(u, un, n * sizeof(double));
    }
    free(F); free(un);
}

static void experiment_A(void) {
    /* profiles after one revolution */
    const int n = 200; const double C = 0.5;
    double *u[3];
    for (int s = 0; s < 3; s++) {
        u[s] = malloc(n * sizeof(double));
        for (int i = 0; i < n; i++) u[s][i] = ic_avg(3, (double)i / n, 1.0 / n);
        advect(s, n, C, 1.0, u[s]);
    }
    FILE *f = out("profiles.csv");
    fprintf(f, "x,exact,upwind,lax-wendroff,van-leer\n");
    for (int i = 0; i < n; i++)
        fprintf(f, "%.5f,%.6f,%.6f,%.6f,%.6f\n", (i + 0.5) / n, ic_avg(3, (double)i / n, 1.0 / n),
                u[0][i], u[1][i], u[2][i]);
    fclose(f);
    printf("A. profiles after one revolution (N = %d, C = %.1f): min / max\n", n, C);
    for (int s = 0; s < 3; s++) {
        double mn = 1e9, mx = -1e9;
        for (int i = 0; i < n; i++) { if (u[s][i] < mn) mn = u[s][i]; if (u[s][i] > mx) mx = u[s][i]; }
        printf("   %-13s min %+.4f  max %.4f\n", sname[s], mn, mx);
        free(u[s]);
    }

    /* L1 error vs grid, smooth (gaussian) and discontinuous (square) */
    f = out("convergence.csv");
    fprintf(f, "profile,scheme,N,L1\n");
    printf("\n   L1 error after one revolution, C = 0.5\n   %-8s %-13s", "profile", "scheme");
    int Ns[] = {50, 100, 200, 400, 800, 1600};
    for (int k = 0; k < 6; k++) printf(" %9d", Ns[k]);
    printf("   order\n");
    for (int which = 2; which >= 1; which--)
        for (int s = 0; s < 3; s++) {
            printf("   %-8s %-13s", which == 2 ? "gaussian" : "square", sname[s]);
            double e[6];
            for (int k = 0; k < 6; k++) {
                int nn = Ns[k]; double *v = malloc(nn * sizeof(double));
                for (int i = 0; i < nn; i++) v[i] = ic_avg(which, (double)i / nn, 1.0 / nn);
                advect(s, nn, 0.5, 1.0, v);
                e[k] = 0;
                for (int i = 0; i < nn; i++) e[k] += fabs(v[i] - ic_avg(which, (double)i / nn, 1.0 / nn)) / nn;
                fprintf(f, "%s,%s,%d,%.6e\n", which == 2 ? "gaussian" : "square", sname[s], nn, e[k]);
                printf(" %9.2e", e[k]);
                free(v);
            }
            printf("   %.2f\n", log(e[4] / e[5]) / log(2.0));
        }
    fclose(f);
}

/* ---------------- B: measured viscosity ---------------- */
static double variance(const double *u, int n) {
    double m0 = 0, m1 = 0, m2 = 0;
    for (int i = 0; i < n; i++) { double x = (i + 0.5) / n; m0 += u[i]; m1 += u[i] * x; m2 += u[i] * x * x; }
    double mu = m1 / m0; return m2 / m0 - mu * mu;
}
static void experiment_B(void) {
    FILE *f = out("viscosity.csv");
    fprintf(f, "N,C,nu_predicted,nu_measured\n");
    printf("\nB. upwind: viscosity measured from the spreading of a gaussian (T = 1)\n");
    printf("   %6s %5s %14s %14s %8s\n", "N", "C", "nu predicted", "nu measured", "ratio");
    int Ns[] = {100, 200, 400, 800};
    double Cs[] = {0.1, 0.5, 0.9, 1.0};
    for (int a = 0; a < 4; a++)
        for (int b = 0; b < 4; b++) {
            int n = Ns[a]; double C = Cs[b], dx = 1.0 / n;
            double *u = malloc(n * sizeof(double));
            for (int i = 0; i < n; i++) u[i] = ic_avg(2, (double)i / n, dx);
            double v0 = variance(u, n);
            advect(UPWIND, n, C, 1.0, u);
            double nu_m = (variance(u, n) - v0) / 2.0, nu_p = dx * (1.0 - C) / 2.0;
            fprintf(f, "%d,%.1f,%.6e,%.6e\n", n, C, nu_p, nu_m);
            if (nu_p > 0) printf("   %6d %5.1f %14.3e %14.3e %8.3f\n", n, C, nu_p, nu_m, nu_m / nu_p);
            else          printf("   %6d %5.1f %14.3e %14.3e %8s\n", n, C, nu_p, nu_m, "-");
            free(u);
        }
    fclose(f);
}

/* ---------------- C: 2D false diffusion ---------------- */
/* steady u phi_x + v phi_y = 0, u,v >= 0, upwind. Inflow: phi = 0 on y = 0,
   and on x = 0 phi = 1 above y0, 0 below (y0 = 0: the step starts at the corner).
   Upwind with positive velocities can be solved by a single sweep. */
static double *solve2d(int n, double theta, double y0) {
    double u = cos(theta), v = sin(theta);
    double *p = malloc((size_t)n * n * sizeof(double));
    for (int j = 0; j < n; j++)
        for (int i = 0; i < n; i++) {
            double w = (i > 0) ? p[j * n + i - 1] : ((j + 0.5) / n > y0 ? 1.0 : 0.0);  /* west inflow */
            double s = (j > 0) ? p[(j - 1) * n + i] : 0.0;   /* south (inflow 0) */
            p[j * n + i] = (u * w + v * s) / (u + v);
        }
    return p;
}

/* 10-90% width of the step across the 45-degree streamline through the corner,
   at distance s from the corner (cells on the anti-diagonal i + j = k) */
static double width45(const double *p, int n, int k) {
    double h = 1.0 / n, lo = NAN, hi = NAN;
    /* walk along the anti-diagonal from the bottom-right (phi ~ 0) to top-left (phi ~ 1) */
    double prev_n = 0, prev_v = -1;
    for (int i = k; i >= 0; i--) {
        int j = k - i;
        if (i >= n || j >= n) continue;
        double nn = (i - j) * h / sqrt(2.0), val = p[j * n + i];   /* normal coordinate */
        if (prev_v >= 0) {
            if (prev_v < 0.1 && val >= 0.1) lo = prev_n + (0.1 - prev_v) / (val - prev_v) * (nn - prev_n);
            if (prev_v < 0.9 && val >= 0.9) hi = prev_n + (0.9 - prev_v) / (val - prev_v) * (nn - prev_n);
        }
        prev_n = nn; prev_v = val;
    }
    return fabs(hi - lo);
}

static void experiment_C(void) {
    double thetas[] = {0.0, 22.5, 45.0};
    for (int t = 0; t < 3; t++) {
        int n = 200; double *p = solve2d(n, thetas[t] * M_PI / 180.0, thetas[t] == 0.0 ? 0.5 : 0.0);
        char name[64]; snprintf(name, sizeof name, "field2d_%02.0f.csv", thetas[t]);
        FILE *f = out(name);
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) fprintf(f, "%s%.4f", i ? "," : "", p[j * n + i]);
            fprintf(f, "\n");
        }
        fclose(f); free(p);
    }
    FILE *f = out("falsediff.csv");
    fprintf(f, "N,s,width_measured,width_predicted\n");
    printf("\nC. 2D upwind, step advected at 45 degrees: 10-90%% width of the smeared step\n");
    printf("   %6s %8s %12s %12s\n", "N", "s", "measured", "predicted");
    int Ns[] = {100, 200, 400};
    for (int a = 0; a < 3; a++) {
        int n = Ns[a]; double h = 1.0 / n;
        double *p = solve2d(n, M_PI / 4, 0.0);
        double fr[] = {0.1, 0.2, 0.4, 0.7, 1.0};
        for (int b = 0; b < 5; b++) {
            double s = fr[b];                                /* distance from corner along the flow */
            int k = (int)llround(s * sqrt(2.0) / h - 1.0);   /* cells with (i+j+1) h = s sqrt2 */
            double sm = (k + 1) * h / sqrt(2.0);
            double w = width45(p, n, k);
            /* modified equation: nu = U h / (2 sqrt2) in every direction ->
               cross-stream variance 2 nu s / U = s h / sqrt2, width = 2.563 sigma */
            double wp = 2.563 * sqrt(sm * h / sqrt(2.0));
            fprintf(f, "%d,%.4f,%.5f,%.5f\n", n, sm, w, wp);
            printf("   %6d %8.3f %12.4f %12.4f\n", n, sm, w, wp);
        }
        free(p);
    }
    /* aligned flow, step at y = 0.5: count cells with 0.01 < phi < 0.99 in the last column */
    double *p = solve2d(200, 0.0, 0.5);
    int smeared = 0;
    for (int j = 0; j < 200; j++) if (p[j * 200 + 199] > 0.01 && p[j * 200 + 199] < 0.99) smeared++;
    printf("   theta = 0 (flow along the grid): smeared cells at the outlet = %d (step stays sharp)\n", smeared);
    free(p);
    fclose(f);
}

int main(int argc, char **argv) {
    if (argc > 1) dir = argv[1];
    experiment_A();
    experiment_B();
    experiment_C();
    return 0;
}
