/**
 * 2D turbulence, pseudo-spectral DNS: decaying, or forced at an intermediate
 * scale k_f (energy injection rate eps) with linear drag alpha at large scales.
 * (Derived from the hyperviscosity code; use p = 1 for ordinary viscosity.)
 *
 * Run: ./turb2d outdir N p nu tend [kf eps alpha]
 *
 *   w_t + u . grad(w) = -nu_p (-lap)^p w,   lap(psi) = -w,  u = (psi_y, -psi_x)
 *
 * p = 1 is ordinary viscosity (Navier-Stokes); p > 1 is hyperviscosity, which
 * concentrates the dissipation at the smallest resolved scales.
 * Periodic box [0, 2pi)^2, N x N points, Fourier pseudo-spectral method with
 * 2/3 dealiasing, RK4 with an integrating factor (the dissipation term is
 * integrated exactly, so even p = 8 does not limit the time step).
 *
 * The same random initial condition (energy peaked at k = 10) is used in every run.
 * nu_p is chosen so that every run damps the largest resolved wavenumber
 * k_max at the same rate: nu_p k_max^(2p) = nu_1 k_max^2.
 *
 * Build: gcc -O3 -march=native hyper.c -o hyper -lfftw3 -lm
 * Run:   ./hyper outdir N p nu1 tend      (e.g. ./hyper out 256 4 5e-4 30)
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
static double nup, kmax;
static fftw_plan fwd, bwd;
static double *rbuf;
static fftw_complex *cbuf;
#define CI(i, j) ((size_t)(j) * NH + (i))       /* complex index: kx = i, ky = j */
#define RI(i, j) ((size_t)(j) * N + (i))
static inline int ky_of(int j) { return j <= N / 2 ? j : j - N; }
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

static double *mask;   /* 2/3 dealiasing */
static double *Dk;     /* dissipation rate nu_p k^(2p) */

static void to_phys(const fftw_complex *in, double *out) {
    memcpy(cbuf, in, sizeof(fftw_complex) * N * NH);
    fftw_execute(bwd);
    double s = 1.0 / ((double)N * N);
    for (size_t k = 0; k < (size_t)N * N; k++) out[k] = rbuf[k] * s;
}
static void to_spec(const double *in, fftw_complex *out) {
    memcpy(rbuf, in, sizeof(double) * N * N);
    fftw_execute(fwd);
    memcpy(out, cbuf, sizeof(fftw_complex) * N * NH);
}

/* nonlinear term  -u.grad(w), dealiased */
static double *ux, *uy, *wx, *wy, *nl;
static fftw_complex *tmp;
static double g_kf = 0, g_fc = 0;   /* forcing shell and amplitude for the current step */
static void nonlinear_raw(const fftw_complex *w, fftw_complex *out);
static void nonlinear(const fftw_complex *w, fftw_complex *out) {
    nonlinear_raw(w, out);
    if (g_kf <= 0) return;
    for (int j = 0; j < N; j++) for (int i = 0; i < NH; i++) {
        double k = sqrt((double)i * i + (double)ky_of(j) * ky_of(j));
        if (fabs(k - g_kf) <= 1.0) out[CI(i, j)] += g_fc * w[CI(i, j)];
    }
}
static void nonlinear_raw(const fftw_complex *w, fftw_complex *out) {
    for (int j = 0; j < N; j++) for (int i = 0; i < NH; i++) {
        double kx = i, ky = ky_of(j), k2 = kx * kx + ky * ky;
        double complex psi = k2 > 0 ? w[CI(i, j)] / k2 : 0;
        tmp[CI(i, j)] = I * ky * psi;            /* u =  psi_y */
    }
    to_phys(tmp, ux);
    for (int j = 0; j < N; j++) for (int i = 0; i < NH; i++) {
        double kx = i, ky = ky_of(j), k2 = kx * kx + ky * ky;
        double complex psi = k2 > 0 ? w[CI(i, j)] / k2 : 0;
        tmp[CI(i, j)] = -I * kx * psi;           /* v = -psi_x */
    }
    to_phys(tmp, uy);
    for (int j = 0; j < N; j++) for (int i = 0; i < NH; i++) tmp[CI(i, j)] = I * (double)i * w[CI(i, j)];
    to_phys(tmp, wx);
    for (int j = 0; j < N; j++) for (int i = 0; i < NH; i++) tmp[CI(i, j)] = I * (double)ky_of(j) * w[CI(i, j)];
    to_phys(tmp, wy);
    for (size_t k = 0; k < (size_t)N * N; k++) nl[k] = -(ux[k] * wx[k] + uy[k] * wy[k]);
    to_spec(nl, out);
    for (size_t k = 0; k < (size_t)N * NH; k++) out[k] *= mask[k];
}

/* energy spectrum E(k) (shells of width 1), total energy and enstrophy */
static void spectrum(const fftw_complex *w, double *Ek, int nk, double *E, double *Z, double *eps) {
    for (int k = 0; k < nk; k++) Ek[k] = 0;
    double e = 0, z = 0, d = 0, norm = 1.0 / ((double)N * N * N * N);
    for (int j = 0; j < N; j++) for (int i = 0; i < NH; i++) {
        double kx = i, ky = ky_of(j), k2 = kx * kx + ky * ky;
        if (k2 == 0) continue;
        double mult = (i == 0 || i == N / 2) ? 1 : 2;   /* Hermitian half */
        double a2 = creal(w[CI(i, j)] * conj(w[CI(i, j)])) * norm * mult;
        double ek = 0.5 * a2 / k2;                         /* |u|^2/2 = |w|^2/(2 k^2) */
        int sh = (int)lround(sqrt(k2));
        if (sh < nk) Ek[sh] += ek;
        e += ek; z += 0.5 * a2; d += 2 * Dk[CI(i, j)] * ek;
    }
    *E = e; *Z = z; *eps = d;
}

int main(int argc, char **argv) {
    const char *dir = argc > 1 ? argv[1] : ".";
    N = argc > 2 ? atoi(argv[2]) : 256;
    P = argc > 3 ? atoi(argv[3]) : 1;
    double nu1 = argc > 4 ? atof(argv[4]) : 5e-4, tend = argc > 5 ? atof(argv[5]) : 30;
    double kf = argc > 6 ? atof(argv[6]) : 0, feps = argc > 7 ? atof(argv[7]) : 0, alpha = argc > 8 ? atof(argv[8]) : 0;
    NH = N / 2 + 1;
    kmax = N / 3.0;                                  /* largest wavenumber kept by the 2/3 rule */
    nup = nu1 * pow(kmax, 2.0 - 2.0 * P);            /* same damping at k_max for every p */

    size_t nc = (size_t)N * NH, nr = (size_t)N * N;
    rbuf = fftw_alloc_real(nr); cbuf = fftw_alloc_complex(nc);
    fwd = fftw_plan_dft_r2c_2d(N, N, rbuf, cbuf, FFTW_MEASURE);
    bwd = fftw_plan_dft_c2r_2d(N, N, cbuf, rbuf, FFTW_MEASURE);
    fftw_complex *w = fftw_alloc_complex(nc), *w1 = fftw_alloc_complex(nc), *a = fftw_alloc_complex(nc),
                 *b = fftw_alloc_complex(nc), *c = fftw_alloc_complex(nc), *d = fftw_alloc_complex(nc);
    tmp = fftw_alloc_complex(nc);
    ux = malloc(nr * 8); uy = malloc(nr * 8); wx = malloc(nr * 8); wy = malloc(nr * 8); nl = malloc(nr * 8);
    mask = malloc(nc * 8); Dk = malloc(nc * 8);
    for (int j = 0; j < N; j++) for (int i = 0; i < NH; i++) {
        double kx = i, ky = ky_of(j);
        mask[CI(i, j)] = (fabs(kx) < kmax && fabs(ky) < kmax) ? 1 : 0;
        Dk[CI(i, j)] = nup * pow(kx * kx + ky * ky, P);
    }

    /* initial condition: random phases, E(k) ~ k^4 exp(-2 (k/k0)^2), k0 = 10, same for every run */
    unsigned long long rng = 12345;
    for (int j = 0; j < N; j++) for (int i = 0; i < NH; i++) {
        rng = rng * 6364136223846793005ULL + 1442695040888963407ULL;
        double ph = 2 * M_PI * (double)(rng >> 11) / 9007199254740992.0;
        double kx = i, ky = ky_of(j), k = sqrt(kx * kx + ky * ky), k0 = kf > 0 ? kf : 10;
        double Ek = k > 0 ? pow(k, 4) * exp(-2 * (k / k0) * (k / k0)) : 0;
        /* |w|^2 spread over the shell: amplitude ~ sqrt(E(k) k^2 / (pi k)) */
        double amp = k > 0 ? sqrt(Ek * k * k / (M_PI * k)) : 0;
        w[CI(i, j)] = amp * cexp(I * ph) * mask[CI(i, j)] * (double)N * N;
    }
    /* make it exactly real and normalize the energy to 0.5 */
    double *wp = malloc(nr * 8);
    to_phys(w, wp); to_spec(wp, w);
    int nk = (int)kmax + 2;
    double *Ek = calloc(nk, 8), E, Z, eps;
    spectrum(w, Ek, nk, &E, &Z, &eps);
    double sc = sqrt((kf > 0 ? 1e-4 : 0.5) / E);
    for (size_t k = 0; k < nc; k++) w[k] *= sc;

    char path[512];
    snprintf(path, sizeof path, "%s/history_N%d_p%d%s.csv", dir, N, P, kf > 0 ? "_forced" : ""); FILE *fh = fopen(path, "w");
    fprintf(fh, "t,energy,enstrophy,dissipation\n");
    snprintf(path, sizeof path, "%s/spectra_N%d_p%d%s.csv", dir, N, P, kf > 0 ? "_forced" : ""); FILE *fs = fopen(path, "w");
    fprintf(fs, "t,k,E,D\n");

    double dt = 0.5 * (2 * M_PI / N) / 3.0;          /* CFL ~0.5 for |u| <~ 3 */
    int steps = (int)lround(tend / dt);
    double *E1 = malloc(nc * 8), *E2 = malloc(nc * 8);
    for (size_t k = 0; k < nc; k++) { E1[k] = exp(-(Dk[k] + alpha) * dt); E2[k] = exp(-(Dk[k] + alpha) * dt / 2); }
    double t0 = now();
    for (int s = 0; s <= steps; s++) {
        double t = s * dt;
        if (s % 50 == 0) {
            spectrum(w, Ek, nk, &E, &Z, &eps);
            fprintf(fh, "%.4f,%.8e,%.8e,%.8e\n", t, E, Z, eps);
        }
        int sample = 0;
        double ts[] = {5, 10, 20, 30, 50};
        for (int q = 0; q < 5; q++) if (fabs(t - ts[q]) < dt / 2) sample = 1;
        if (kf > 0 && s > 0 && fmod(t + dt / 2, 5.0) < dt) sample = 1;
        if (sample || s == 0) {               /* vorticity snapshots at the sampling times (and t = 0) */
            static double *snap = NULL; if (!snap) snap = malloc(nr * 8);
            to_phys(w, snap);
            snprintf(path, sizeof path, "%s/vort_N%d_p%d%s_t%02d.csv", dir, N, P, kf > 0 ? "_forced" : "", (int)lround(t));
            FILE *fq = fopen(path, "w");
            for (int jj = 0; jj < N; jj += (N >= 512 ? 2 : 1)) { for (int ii = 0; ii < N; ii += (N >= 512 ? 2 : 1)) fprintf(fq, "%s%.4e", ii ? "," : "", snap[RI(ii, jj)]); fprintf(fq, "\n"); }
            fclose(fq);
        }
        if (sample) {
            spectrum(w, Ek, nk, &E, &Z, &eps);
            for (int k = 1; k < nk; k++) {
                double Dsh = 0;   /* dissipation spectrum 2 nu_p k^(2p) E(k) */
                Dsh = 2 * nup * pow((double)k, 2.0 * P) * Ek[k];
                fprintf(fs, "%.1f,%d,%.8e,%.8e\n", t, k, Ek[k], Dsh);
            }
        }
        if (s == steps) break;
        /* forcing: f = eps w / (2 E_f) on the shell |k - kf| <= 1 (energy injection rate eps) */
        if (kf > 0) {
            double Ef = 0, norm = 1.0 / ((double)N * N * N * N);
            for (int j = 0; j < N; j++) for (int i = 0; i < NH; i++) {
                double kk2 = (double)i * i + (double)ky_of(j) * ky_of(j), k = sqrt(kk2);
                if (kk2 > 0 && fabs(k - kf) <= 1.0)
                    Ef += 0.5 * creal(w[CI(i, j)] * conj(w[CI(i, j)])) * norm * ((i == 0 || i == N / 2) ? 1 : 2) / kk2;
            }
            g_kf = kf; g_fc = feps / (2 * Ef);
            if (g_fc > 2.0) g_fc = 2.0;      /* cap while the forced shell is still almost empty */
        }
        /* RK4 with integrating factor */
        nonlinear(w, a);
        for (size_t k = 0; k < nc; k++) w1[k] = E2[k] * (w[k] + 0.5 * dt * a[k]);
        nonlinear(w1, b);
        for (size_t k = 0; k < nc; k++) w1[k] = E2[k] * w[k] + 0.5 * dt * b[k];
        nonlinear(w1, c);
        for (size_t k = 0; k < nc; k++) w1[k] = E1[k] * w[k] + dt * E2[k] * c[k];
        nonlinear(w1, d);
        for (size_t k = 0; k < nc; k++)
            w[k] = E1[k] * w[k] + dt / 6 * (E1[k] * a[k] + 2 * E2[k] * (b[k] + c[k]) + d[k]);
    }
    double el = now() - t0;
    fclose(fh); fclose(fs);
    /* final vorticity field */
    to_phys(w, wp);
    snprintf(path, sizeof path, "%s/vorticity_N%d_p%d%s.csv", dir, N, P, kf > 0 ? "_forced" : ""); FILE *fv = fopen(path, "w");
    for (int j = 0; j < N; j++) { for (int i = 0; i < N; i++) fprintf(fv, "%s%.5e", i ? "," : "", wp[RI(i, j)]); fprintf(fv, "\n"); }
    fclose(fv);
    spectrum(w, Ek, nk, &E, &Z, &eps);
    printf("N = %d  p = %d  nu_p = %.3e  t = %.1f  energy %.5f  enstrophy %.4f  (%d steps, %.1f s)\n",
           N, P, nup, tend, E, Z, steps, el);
    return 0;
}
