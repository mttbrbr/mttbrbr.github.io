/**
 * Taylor-Green vortex at Re = 1600: how much dissipation comes from the scheme?
 *
 * 3D incompressible Navier-Stokes, periodic box [0, 2pi)^3, staggered (MAC)
 * finite differences, like a finite-volume code on a Cartesian mesh.
 *  - convection in divergence form, face values from:
 *      0 central        (second order, conserves kinetic energy)
 *      1 upwind         (first order)
 *      2 linearUpwind   (second-order upwind)
 *  - diffusion: 7-point Laplacian
 *  - pressure: exact projection with FFTs (discrete Laplacian eigenvalues)
 *  - time: SSP-RK3, projection after every stage
 * Initial condition: u = sin x cos y cos z, v = -cos x sin y cos z, w = 0.
 * Output: kinetic energy, total dissipation -dE/dt and the physical part
 * nu <|grad u|^2>; their difference is the numerical dissipation of the scheme.
 *
 * Build: gcc -O3 -march=native tgv.c -o tgv -lfftw3 -lm
 * Run:   ./tgv outdir N scheme Re tend
 */
#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <time.h>

static int N, SCH; static double h, nu;
static size_t NR, NCX;
#define ID(i, j, k) ((((size_t)(((k) + N) % N)) * N + (size_t)(((j) + N) % N)) * N + (size_t)(((i) + N) % N))
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

/* face value of q between qL (upstream if a > 0) and qR, with qLL, qRR for linearUpwind */
static inline double fval(double a, double qLL, double qL, double qR, double qRR) {
    switch (SCH) {
    case 0: return 0.5 * (qL + qR);
    case 1: return a > 0 ? qL : qR;
    default: return a > 0 ? qL + 0.5 * (qL - qLL) : qR + 0.5 * (qR - qRR);
    }
}

/* convection + diffusion for all three components: R = -div(u q) + nu lap(q) */
static void rhs(const double *u, const double *v, const double *w, double *Ru, double *Rv, double *Rw) {
    double ih = 1.0 / h, ih2 = 1.0 / (h * h);
    for (int k = 0; k < N; k++) for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
        size_t c = ID(i, j, k);
        /* ---- u at (i h, (j+.5) h, (k+.5) h) ---- */
        {
            double a, fe, fw, fn, fs, ft, fb;
            a = 0.5 * (u[c] + u[ID(i + 1, j, k)]);
            fe = a * fval(a, u[ID(i - 1, j, k)], u[c], u[ID(i + 1, j, k)], u[ID(i + 2, j, k)]);
            a = 0.5 * (u[ID(i - 1, j, k)] + u[c]);
            fw = a * fval(a, u[ID(i - 2, j, k)], u[ID(i - 1, j, k)], u[c], u[ID(i + 1, j, k)]);
            a = 0.5 * (v[ID(i - 1, j + 1, k)] + v[ID(i, j + 1, k)]);
            fn = a * fval(a, u[ID(i, j - 1, k)], u[c], u[ID(i, j + 1, k)], u[ID(i, j + 2, k)]);
            a = 0.5 * (v[ID(i - 1, j, k)] + v[ID(i, j, k)]);
            fs = a * fval(a, u[ID(i, j - 2, k)], u[ID(i, j - 1, k)], u[c], u[ID(i, j + 1, k)]);
            a = 0.5 * (w[ID(i - 1, j, k + 1)] + w[ID(i, j, k + 1)]);
            ft = a * fval(a, u[ID(i, j, k - 1)], u[c], u[ID(i, j, k + 1)], u[ID(i, j, k + 2)]);
            a = 0.5 * (w[ID(i - 1, j, k)] + w[ID(i, j, k)]);
            fb = a * fval(a, u[ID(i, j, k - 2)], u[ID(i, j, k - 1)], u[c], u[ID(i, j, k + 1)]);
            double lap = (u[ID(i + 1, j, k)] + u[ID(i - 1, j, k)] + u[ID(i, j + 1, k)] + u[ID(i, j - 1, k)]
                        + u[ID(i, j, k + 1)] + u[ID(i, j, k - 1)] - 6 * u[c]) * ih2;
            Ru[c] = -((fe - fw) + (fn - fs) + (ft - fb)) * ih + nu * lap;
        }
        /* ---- v at ((i+.5) h, j h, (k+.5) h) ---- */
        {
            double a, fe, fw, fn, fs, ft, fb;
            a = 0.5 * (u[ID(i + 1, j - 1, k)] + u[ID(i + 1, j, k)]);
            fe = a * fval(a, v[ID(i - 1, j, k)], v[c], v[ID(i + 1, j, k)], v[ID(i + 2, j, k)]);
            a = 0.5 * (u[ID(i, j - 1, k)] + u[ID(i, j, k)]);
            fw = a * fval(a, v[ID(i - 2, j, k)], v[ID(i - 1, j, k)], v[c], v[ID(i + 1, j, k)]);
            a = 0.5 * (v[c] + v[ID(i, j + 1, k)]);
            fn = a * fval(a, v[ID(i, j - 1, k)], v[c], v[ID(i, j + 1, k)], v[ID(i, j + 2, k)]);
            a = 0.5 * (v[ID(i, j - 1, k)] + v[c]);
            fs = a * fval(a, v[ID(i, j - 2, k)], v[ID(i, j - 1, k)], v[c], v[ID(i, j + 1, k)]);
            a = 0.5 * (w[ID(i, j - 1, k + 1)] + w[ID(i, j, k + 1)]);
            ft = a * fval(a, v[ID(i, j, k - 1)], v[c], v[ID(i, j, k + 1)], v[ID(i, j, k + 2)]);
            a = 0.5 * (w[ID(i, j - 1, k)] + w[ID(i, j, k)]);
            fb = a * fval(a, v[ID(i, j, k - 2)], v[ID(i, j, k - 1)], v[c], v[ID(i, j, k + 1)]);
            double lap = (v[ID(i + 1, j, k)] + v[ID(i - 1, j, k)] + v[ID(i, j + 1, k)] + v[ID(i, j - 1, k)]
                        + v[ID(i, j, k + 1)] + v[ID(i, j, k - 1)] - 6 * v[c]) * ih2;
            Rv[c] = -((fe - fw) + (fn - fs) + (ft - fb)) * ih + nu * lap;
        }
        /* ---- w at ((i+.5) h, (j+.5) h, k h) ---- */
        {
            double a, fe, fw, fn, fs, ft, fb;
            a = 0.5 * (u[ID(i + 1, j, k - 1)] + u[ID(i + 1, j, k)]);
            fe = a * fval(a, w[ID(i - 1, j, k)], w[c], w[ID(i + 1, j, k)], w[ID(i + 2, j, k)]);
            a = 0.5 * (u[ID(i, j, k - 1)] + u[ID(i, j, k)]);
            fw = a * fval(a, w[ID(i - 2, j, k)], w[ID(i - 1, j, k)], w[c], w[ID(i + 1, j, k)]);
            a = 0.5 * (v[ID(i, j + 1, k - 1)] + v[ID(i, j + 1, k)]);
            fn = a * fval(a, w[ID(i, j - 1, k)], w[c], w[ID(i, j + 1, k)], w[ID(i, j + 2, k)]);
            a = 0.5 * (v[ID(i, j, k - 1)] + v[ID(i, j, k)]);
            fs = a * fval(a, w[ID(i, j - 2, k)], w[ID(i, j - 1, k)], w[c], w[ID(i, j + 1, k)]);
            a = 0.5 * (w[c] + w[ID(i, j, k + 1)]);
            ft = a * fval(a, w[ID(i, j, k - 1)], w[c], w[ID(i, j, k + 1)], w[ID(i, j, k + 2)]);
            a = 0.5 * (w[ID(i, j, k - 1)] + w[c]);
            fb = a * fval(a, w[ID(i, j, k - 2)], w[ID(i, j, k - 1)], w[c], w[ID(i, j, k + 1)]);
            double lap = (w[ID(i + 1, j, k)] + w[ID(i - 1, j, k)] + w[ID(i, j + 1, k)] + w[ID(i, j - 1, k)]
                        + w[ID(i, j, k + 1)] + w[ID(i, j, k - 1)] - 6 * w[c]) * ih2;
            Rw[c] = -((fe - fw) + (fn - fs) + (ft - fb)) * ih + nu * lap;
        }
    }
}

/* exact discrete projection: solve lap(phi) = div(u), u -= grad(phi) */
static double *divb, *eig; static fftw_complex *spec; static fftw_plan pf, pb;
static void project(double *u, double *v, double *w) {
    for (int k = 0; k < N; k++) for (int j = 0; j < N; j++) for (int i = 0; i < N; i++)
        divb[ID(i, j, k)] = (u[ID(i + 1, j, k)] - u[ID(i, j, k)] + v[ID(i, j + 1, k)] - v[ID(i, j, k)]
                           + w[ID(i, j, k + 1)] - w[ID(i, j, k)]) / h;
    fftw_execute(pf);
    for (size_t q = 0; q < NCX; q++) spec[q] = eig[q] != 0 ? spec[q] / eig[q] / (double)NR : 0;
    fftw_execute(pb);                                   /* divb now holds phi */
    for (int k = 0; k < N; k++) for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
        size_t c = ID(i, j, k);
        u[c] -= (divb[c] - divb[ID(i - 1, j, k)]) / h;
        v[c] -= (divb[c] - divb[ID(i, j - 1, k)]) / h;
        w[c] -= (divb[c] - divb[ID(i, j, k - 1)]) / h;
    }
}
static void energy(const double *u, const double *v, const double *w, double *E, double *ephys) {
    double e = 0, g = 0;
    for (int k = 0; k < N; k++) for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
        size_t c = ID(i, j, k);
        e += u[c] * u[c] + v[c] * v[c] + w[c] * w[c];
        const double *q[3] = {u, v, w};
        for (int m = 0; m < 3; m++) {
            double dx = q[m][ID(i + 1, j, k)] - q[m][c], dy = q[m][ID(i, j + 1, k)] - q[m][c], dz = q[m][ID(i, j, k + 1)] - q[m][c];
            g += dx * dx + dy * dy + dz * dz;
        }
    }
    *E = 0.5 * e / NR; *ephys = nu * g / (h * h) / NR;
}

int main(int argc, char **argv) {
    const char *dir = argc > 1 ? argv[1] : ".";
    N = argc > 2 ? atoi(argv[2]) : 64; SCH = argc > 3 ? atoi(argv[3]) : 0;
    double Re = argc > 4 ? atof(argv[4]) : 1600, tend = argc > 5 ? atof(argv[5]) : 20;
    nu = Re > 0 ? 1.0 / Re : 0; h = 2 * M_PI / N; NR = (size_t)N * N * N; NCX = (size_t)N * N * (N / 2 + 1);
    const char *names[] = {"central", "upwind", "linearUpwind"};
    double *u = malloc(NR * 8), *v = malloc(NR * 8), *w = malloc(NR * 8);
    double *u0 = malloc(NR * 8), *v0 = malloc(NR * 8), *w0 = malloc(NR * 8);
    double *Ru = malloc(NR * 8), *Rv = malloc(NR * 8), *Rw = malloc(NR * 8);
    divb = fftw_alloc_real(NR); spec = fftw_alloc_complex(NCX); eig = malloc(NCX * 8);
    pf = fftw_plan_dft_r2c_3d(N, N, N, divb, spec, FFTW_MEASURE);
    pb = fftw_plan_dft_c2r_3d(N, N, N, spec, divb, FFTW_MEASURE);
    for (int k = 0; k < N; k++) for (int j = 0; j < N; j++) for (int i = 0; i <= N / 2; i++)
        eig[((size_t)k * N + j) * (N / 2 + 1) + i] = -4 / (h * h) * (pow(sin(M_PI * i / N), 2) + pow(sin(M_PI * j / N), 2) + pow(sin(M_PI * k / N), 2));
    for (int k = 0; k < N; k++) for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
        size_t c = ID(i, j, k);
        double xc = (i + 0.5) * h, yc = (j + 0.5) * h, zc = (k + 0.5) * h, xf = i * h, yf = j * h;
        u[c] = sin(xf) * cos(yc) * cos(zc);
        v[c] = -cos(xc) * sin(yf) * cos(zc);
        w[c] = 0;
    }
    project(u, v, w);
    char path[512]; snprintf(path, sizeof path, "%s/tgv_N%d_%s_Re%g.csv", dir, N, names[SCH], Re);
    FILE *f = fopen(path, "w"); fprintf(f, "t,energy,physical_dissipation\n");
    double dt = 0.25 * h, t = 0, E, ep, t0 = now();
    int steps = (int)lround(tend / dt);
    for (int s = 0; s <= steps; s++) {
        if (s % 5 == 0) { energy(u, v, w, &E, &ep); fprintf(f, "%.5f,%.8e,%.8e\n", t, E, ep); }
        if (s == steps) break;
        memcpy(u0, u, NR * 8); memcpy(v0, v, NR * 8); memcpy(w0, w, NR * 8);
        rhs(u, v, w, Ru, Rv, Rw);
        for (size_t c = 0; c < NR; c++) { u[c] = u0[c] + dt * Ru[c]; v[c] = v0[c] + dt * Rv[c]; w[c] = w0[c] + dt * Rw[c]; }
        project(u, v, w);
        rhs(u, v, w, Ru, Rv, Rw);
        for (size_t c = 0; c < NR; c++) { u[c] = 0.75 * u0[c] + 0.25 * (u[c] + dt * Ru[c]); v[c] = 0.75 * v0[c] + 0.25 * (v[c] + dt * Rv[c]); w[c] = 0.75 * w0[c] + 0.25 * (w[c] + dt * Rw[c]); }
        project(u, v, w);
        rhs(u, v, w, Ru, Rv, Rw);
        for (size_t c = 0; c < NR; c++) { u[c] = u0[c] / 3 + 2.0 / 3 * (u[c] + dt * Ru[c]); v[c] = v0[c] / 3 + 2.0 / 3 * (v[c] + dt * Rv[c]); w[c] = w0[c] / 3 + 2.0 / 3 * (w[c] + dt * Rw[c]); }
        project(u, v, w);
        t += dt;
    }
    fclose(f);
    energy(u, v, w, &E, &ep);
    printf("N = %3d  %-12s Re = %g  t = %.1f  E = %.5f  (%d steps, %.0f s)\n", N, names[SCH], Re, tend, E, steps, now() - t0);
    return 0;
}
