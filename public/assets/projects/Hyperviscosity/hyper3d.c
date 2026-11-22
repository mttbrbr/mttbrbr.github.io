/**
 * Hyperviscosity in a pseudo-spectral DNS of 3D forced isotropic turbulence.
 *
 *   u_t = u x w - grad(p + |u|^2/2) - nu_p (-lap)^p u + f,   div u = 0
 *
 * p = 1: Navier-Stokes. p > 1: hyperviscosity.
 * Periodic box [0, 2pi)^3, N^3 points, Fourier pseudo-spectral (rotational form),
 * 2/3 dealiasing, projection onto divergence-free fields, RK4 with integrating
 * factor (dissipation integrated exactly). Forcing: constant energy injection
 * rate eps in the shells |k| <= 2.5 (f = eps u / (2 E_forced)).
 * nu_p chosen so that every p damps k_max at the same rate: nu_p k_max^(2p) = nu_1 k_max^2.
 * Spectra averaged over the statistically steady part of the run.
 *
 * Build: gcc -O3 -march=native hyper3d.c -o hyper3d -lfftw3 -lm
 * Run:   ./hyper3d outdir N p nu1 tend tavg     (e.g. ./hyper3d out 64 1 8e-3 30 15)
 */
#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <time.h>

static int N, NH, P;
static size_t NC, NR;
static double nup, kmax, EPS = 0.1;
static fftw_plan fwd, bwd;
static double *rbuf;
static fftw_complex *cbuf;
#define CI(i, j, l) (((size_t)(l) * N + (j)) * NH + (i))      /* kx = i, ky = j, kz = l */
static inline int kw(int j) { return j <= N / 2 ? j : j - N; }
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

static double *mask, *K2, *Dk;
static double *KX, *KY, *KZ;

static void to_phys(const fftw_complex *in, double *out) {
    memcpy(cbuf, in, sizeof(fftw_complex) * NC); fftw_execute(bwd);
    double s = 1.0 / (double)NR;
    for (size_t k = 0; k < NR; k++) out[k] = rbuf[k] * s;
}
static void to_spec(const double *in, fftw_complex *out) {
    memcpy(rbuf, in, sizeof(double) * NR); fftw_execute(fwd); memcpy(out, cbuf, sizeof(fftw_complex) * NC);
}
static void project(fftw_complex *a, fftw_complex *b, fftw_complex *c) {
    for (size_t k = 0; k < NC; k++) {
        if (K2[k] == 0) { a[k] = b[k] = c[k] = 0; continue; }
        double complex kd = (KX[k] * a[k] + KY[k] * b[k] + KZ[k] * c[k]) / K2[k];
        a[k] -= KX[k] * kd; b[k] -= KY[k] * kd; c[k] -= KZ[k] * kd;
    }
}

static double *u[3], *w[3], *nl[3];
static fftw_complex *tc;
/* nonlinear term: P[ (u x w) ]^, dealiased */
static void nonlinear(fftw_complex *const U[3], fftw_complex *R[3]) {
    for (int c = 0; c < 3; c++) to_phys(U[c], u[c]);
    /* vorticity w = i k x u */
    for (size_t k = 0; k < NC; k++) tc[k] = I * (KY[k] * U[2][k] - KZ[k] * U[1][k]);
    to_phys(tc, w[0]);
    for (size_t k = 0; k < NC; k++) tc[k] = I * (KZ[k] * U[0][k] - KX[k] * U[2][k]);
    to_phys(tc, w[1]);
    for (size_t k = 0; k < NC; k++) tc[k] = I * (KX[k] * U[1][k] - KY[k] * U[0][k]);
    to_phys(tc, w[2]);
    for (size_t k = 0; k < NR; k++) {
        nl[0][k] = u[1][k] * w[2][k] - u[2][k] * w[1][k];
        nl[1][k] = u[2][k] * w[0][k] - u[0][k] * w[2][k];
        nl[2][k] = u[0][k] * w[1][k] - u[1][k] * w[0][k];
    }
    for (int c = 0; c < 3; c++) { to_spec(nl[c], R[c]); for (size_t k = 0; k < NC; k++) R[c][k] *= mask[k]; }
    project(R[0], R[1], R[2]);
}
/* energy spectrum (shells), total energy, dissipation, energy in the forced shells */
static double spectrum(fftw_complex *const U[3], double *Ek, int nk, double *E, double *diss) {
    for (int k = 0; k < nk; k++) Ek[k] = 0;
    double e = 0, d = 0, ef = 0, norm = 1.0 / ((double)NR * NR);
    for (int l = 0; l < N; l++) for (int j = 0; j < N; j++) for (int i = 0; i < NH; i++) {
        size_t k = CI(i, j, l);
        if (K2[k] == 0) continue;
        double mult = (i == 0 || i == N / 2) ? 1 : 2;
        double a2 = 0; for (int c = 0; c < 3; c++) a2 += creal(U[c][k] * conj(U[c][k]));
        double ek = 0.5 * a2 * norm * mult;
        int sh = (int)lround(sqrt(K2[k]));
        if (sh < nk) Ek[sh] += ek;
        e += ek; d += 2 * Dk[k] * ek;
        if (K2[k] <= 6.25) ef += ek;
    }
    *E = e; *diss = d; return ef;
}

int main(int argc, char **argv) {
    const char *dir = argc > 1 ? argv[1] : ".";
    N = argc > 2 ? atoi(argv[2]) : 64; P = argc > 3 ? atoi(argv[3]) : 1;
    double nu1 = argc > 4 ? atof(argv[4]) : 8e-3, tend = argc > 5 ? atof(argv[5]) : 30, tavg = argc > 6 ? atof(argv[6]) : 15;
    NH = N / 2 + 1; NC = (size_t)N * N * NH; NR = (size_t)N * N * N;
    kmax = N / 3.0; nup = nu1 * pow(kmax, 2.0 - 2.0 * P);
    rbuf = fftw_alloc_real(NR); cbuf = fftw_alloc_complex(NC);
    fwd = fftw_plan_dft_r2c_3d(N, N, N, rbuf, cbuf, FFTW_MEASURE);
    bwd = fftw_plan_dft_c2r_3d(N, N, N, cbuf, rbuf, FFTW_MEASURE);
    mask = malloc(NC * 8); K2 = malloc(NC * 8); Dk = malloc(NC * 8); KX = malloc(NC * 8); KY = malloc(NC * 8); KZ = malloc(NC * 8);
    for (int l = 0; l < N; l++) for (int j = 0; j < N; j++) for (int i = 0; i < NH; i++) {
        size_t k = CI(i, j, l);
        KX[k] = i; KY[k] = kw(j); KZ[k] = kw(l);
        K2[k] = KX[k] * KX[k] + KY[k] * KY[k] + KZ[k] * KZ[k];
        mask[k] = (fabs(KX[k]) < kmax && fabs(KY[k]) < kmax && fabs(KZ[k]) < kmax) ? 1 : 0;
        Dk[k] = nup * pow(K2[k], P);
    }
    fftw_complex *U[3], *W1[3], *A[3], *B[3], *Cc[3], *D[3];
    for (int c = 0; c < 3; c++) {
        U[c] = fftw_alloc_complex(NC); W1[c] = fftw_alloc_complex(NC); A[c] = fftw_alloc_complex(NC);
        B[c] = fftw_alloc_complex(NC); Cc[c] = fftw_alloc_complex(NC); D[c] = fftw_alloc_complex(NC);
        u[c] = malloc(NR * 8); w[c] = malloc(NR * 8); nl[c] = malloc(NR * 8);
    }
    tc = fftw_alloc_complex(NC);

    /* initial condition: random solenoidal field, E(k) ~ k^4 exp(-2 (k/3)^2), same for every run */
    unsigned long long rng = 2024;
    for (size_t k = 0; k < NC; k++) {
        double kk = sqrt(K2[k]), amp = kk > 0 ? sqrt(pow(kk, 4) * exp(-2 * kk * kk / 9.0) / (4 * M_PI * kk * kk)) : 0;
        for (int c = 0; c < 3; c++) {
            rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
            double ph = 2 * M_PI * (double)(rng >> 11) / 9007199254740992.0;
            U[c][k] = amp * cexp(I * ph) * mask[k] * (double)NR;
        }
    }
    project(U[0], U[1], U[2]);
    for (int c = 0; c < 3; c++) { to_phys(U[c], u[c]); to_spec(u[c], U[c]); }   /* make it real */
    int nk = (int)kmax + 2;
    double *Ek = calloc(nk, 8), *Eavg = calloc(nk, 8), E, diss, Ef;
    spectrum(U, Ek, nk, &E, &diss);
    double sc = sqrt(0.5 / E);
    for (int c = 0; c < 3; c++) for (size_t k = 0; k < NC; k++) U[c][k] *= sc;

    char path[512];
    snprintf(path, sizeof path, "%s/history3d_N%d_p%d.csv", dir, N, P); FILE *fh = fopen(path, "w");
    fprintf(fh, "t,energy,dissipation\n");
    double dx = 2 * M_PI / N, dt = 0.3 * dx / 2.0;     /* CFL ~0.3 for |u| <~ 2 */
    int steps = (int)lround(tend / dt), navg = 0;
    double *E1 = malloc(NC * 8), *E2 = malloc(NC * 8);
    for (size_t k = 0; k < NC; k++) { E1[k] = exp(-Dk[k] * dt); E2[k] = exp(-Dk[k] * dt / 2); }
    double t0 = now(), Eacc = 0, Dacc = 0;
    for (int s = 0; s <= steps; s++) {
        double t = s * dt;
        if (s % 10 == 0) {
            spectrum(U, Ek, nk, &E, &diss);
            fprintf(fh, "%.4f,%.8e,%.8e\n", t, E, diss); fflush(fh);
            if (t >= tend - tavg) { for (int k = 0; k < nk; k++) Eavg[k] += Ek[k]; Eacc += E; Dacc += diss; navg++; }
        }
        if (s == steps) break;
        /* forcing: f = eps u / (2 E_forced) on |k| <= 2.5, constant during the step */
        Ef = spectrum(U, Ek, nk, &E, &diss);
        double fa = EPS / (2 * Ef);
        /* RK4 with integrating factor; forcing added to the nonlinear term */
        #define RHS(X, R) do { nonlinear(X, R); for (int c_ = 0; c_ < 3; c_++) for (size_t k_ = 0; k_ < NC; k_++) \
                                 if (K2[k_] <= 6.25 && K2[k_] > 0) R[c_][k_] += fa * X[c_][k_]; } while (0)
        RHS(U, A);
        for (int c = 0; c < 3; c++) for (size_t k = 0; k < NC; k++) W1[c][k] = E2[k] * (U[c][k] + 0.5 * dt * A[c][k]);
        RHS(W1, B);
        for (int c = 0; c < 3; c++) for (size_t k = 0; k < NC; k++) W1[c][k] = E2[k] * U[c][k] + 0.5 * dt * B[c][k];
        RHS(W1, Cc);
        for (int c = 0; c < 3; c++) for (size_t k = 0; k < NC; k++) W1[c][k] = E1[k] * U[c][k] + dt * E2[k] * Cc[c][k];
        RHS(W1, D);
        for (int c = 0; c < 3; c++) for (size_t k = 0; k < NC; k++)
            U[c][k] = E1[k] * U[c][k] + dt / 6 * (E1[k] * A[c][k] + 2 * E2[k] * (B[c][k] + Cc[c][k]) + D[c][k]);
    }
    double el = now() - t0;
    fclose(fh);
    /* averaged spectrum and derived quantities */
    double Em = Eacc / navg, epsm = Dacc / navg, urms = sqrt(2 * Em / 3);
    double lam = 0, num = 0;
    for (int k = 1; k < nk; k++) { Eavg[k] /= navg; num += k * k * Eavg[k]; }
    lam = sqrt(15 * 0 + 5 * Em / num);                 /* Taylor microscale: lambda^2 = 5 E / int k^2 E */
    double nu_eff = epsm / (2 * num);                   /* viscosity that would give the same eps */
    double Rel = urms * lam / nu_eff, eta = pow(nu_eff * nu_eff * nu_eff / epsm, 0.25);
    snprintf(path, sizeof path, "%s/spectrum3d_N%d_p%d.csv", dir, N, P); FILE *fs = fopen(path, "w");
    fprintf(fs, "k,E,D\n");
    for (int k = 1; k < nk; k++) fprintf(fs, "%d,%.8e,%.8e\n", k, Eavg[k], 2 * nup * pow((double)k, 2.0 * P) * Eavg[k]);
    fclose(fs);
    printf("N = %3d p = %d nu_p = %.3e | <E> = %.4f <eps> = %.4f (injected %.3f) | u_rms %.3f Re_lambda %.0f "
           "(with nu_eff = %.2e) k_max eta_eff = %.2f | %d steps, %.0f s\n",
           N, P, nup, Em, epsm, EPS, urms, Rel, nu_eff, kmax * eta, steps, el);
    return 0;
}
