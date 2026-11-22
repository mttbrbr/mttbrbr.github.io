/**
 * Time-stepping schemes: stable is not the same as accurate.
 *
 * Space is treated exactly (Fourier modes, DFT), so every error below comes
 * from the time integration alone. Periodic domain [0, 1), N points.
 *
 *  A. Advection u_t + a u_x = 0 of a Gaussian, one full period, Courant
 *     numbers C = a dt / dx from 0.1 to 10.
 *  B. Diffusion u_t = nu u_xx of a step, 3-point Laplacian (as in a finite-volume
 *     code), large time step (nu dt / dx^2 = 20): which schemes oscillate.
 *
 * Schemes (for each Fourier mode, u' = lambda u):
 *   euler-explicit, rk4,
 *   euler-implicit           (OpenFOAM: Euler)
 *   bdf2                     (OpenFOAM: backward)
 *   crank-nicolson theta=0.5 (OpenFOAM: CrankNicolson 1)
 *   theta = 0.55             (slightly off-centred, like CrankNicolson 0.9)
 *
 * Build: gcc -O3 -march=native timestep.c -o timestep -lm   Run: ./timestep outdir
 */
#define _DEFAULT_SOURCE
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#define N 256
enum { EE, RK4, EI, BDF2, CN, TH55, NS };
static const char *name[NS] = {"euler-explicit", "rk4", "euler-implicit", "bdf2", "crank-nicolson", "theta-0.55"};

static void dft(const double *u, double complex *U) {
    for (int k = 0; k < N; k++) { double complex s = 0; for (int j = 0; j < N; j++) s += u[j] * cexp(-2 * M_PI * I * k * j / N); U[k] = s / N; }
}
static void idft(const double complex *U, double *u) {
    for (int j = 0; j < N; j++) { double complex s = 0; for (int k = 0; k < N; k++) s += U[k] * cexp(2 * M_PI * I * k * j / N); u[j] = creal(s); }
}
static int wavenumber(int k) { return k <= N / 2 ? k : k - N; }

/* advance mode amplitude: z = lambda dt; for BDF2 keep the previous value */
static double complex step(int s, double complex z, double complex u, double complex *uprev, int first) {
    switch (s) {
    case EE:   return u * (1 + z);
    case RK4:  return u * (1 + z + z * z / 2 + z * z * z / 6 + z * z * z * z / 24);
    case EI:   return u / (1 - z);
    case CN:   return u * (1 + 0.5 * z) / (1 - 0.5 * z);
    case TH55: return u * (1 + 0.45 * z) / (1 - 0.55 * z);
    case BDF2: {
        double complex un;
        if (first) un = u / (1 - z);                       /* start with one implicit Euler step */
        else un = (4 * u - *uprev) / (3 - 2 * z);
        *uprev = u; return un;
    }
    }
    return u;
}
static void run(int s, double complex *U0, int nsteps, double complex (*lam)(int), double dt, double *out) {
    double complex U[N], Up[N];
    memcpy(U, U0, sizeof U);
    for (int n = 0; n < nsteps; n++)
        for (int k = 0; k < N; k++) U[k] = step(s, lam(wavenumber(k)) * dt, U[k], &Up[k], n == 0);
    idft(U, out);
}
static double a = 1.0, nu = 1e-3;
static double complex lam_adv(int k) { return -I * 2 * M_PI * k * a; }
/* diffusion: eigenvalues of the 3-point finite-volume Laplacian, as in a real code */
static double complex lam_dif(int k) { double dx = 1.0 / N; return -4 * nu / (dx * dx) * pow(sin(M_PI * k / N), 2); }

int main(int argc, char **argv) {
    const char *dir = argc > 1 ? argv[1] : ".";
    char path[512];
    double dx = 1.0 / N, u0[N], u[N];
    double complex U0[N];

    /* A. advection of a Gaussian, one period */
    for (int j = 0; j < N; j++) u0[j] = exp(-0.5 * pow((j * dx - 0.5) / 0.04, 2));
    dft(u0, U0);
    snprintf(path, sizeof path, "%s/advection.csv", dir); FILE *f = fopen(path, "w");
    fprintf(f, "courant,scheme,peak,phase_shift,l2_error\n");
    snprintf(path, sizeof path, "%s/advection_profiles.csv", dir); FILE *fp = fopen(path, "w");
    fprintf(fp, "courant,scheme,x,u\n");
    double Cs[] = {0.1, 0.25, 0.5, 1, 2, 5, 10};
    printf("A. advection, one period: peak height (exact 1) / phase error (fraction of domain)\n%8s", "Courant");
    for (int s = 0; s < NS; s++) printf(" %22s", name[s]);
    printf("\n");
    for (int c = 0; c < 7; c++) {
        double C = Cs[c], dt = C * dx / a; int nsteps = (int)lround(1.0 / dt);
        dt = 1.0 / nsteps;                                  /* exactly one period */
        printf("%8.2f", C);
        for (int s = 0; s < NS; s++) {
            run(s, U0, nsteps, lam_adv, dt, u);
            double mx = -1e300, xm = 0, e2 = 0, n2 = 0; int jm = 0;
            for (int j = 0; j < N; j++) { if (u[j] > mx) { mx = u[j]; jm = j; } e2 += (u[j] - u0[j]) * (u[j] - u0[j]); n2 += u0[j] * u0[j]; }
            /* sub-cell peak position by parabolic fit */
            double ym = u[(jm - 1 + N) % N], yp = u[(jm + 1) % N], den = ym - 2 * mx + yp;
            xm = (jm + (den != 0 ? 0.5 * (ym - yp) / den : 0)) * dx;
            double ph = xm - 0.5;
            int ok = isfinite(mx) && mx < 1e6;
            if (ok) printf("   %8.3f / %+9.4f", mx, ph); else printf(" %22s", "unstable");
            fprintf(f, "%g,%s,%.6e,%.6e,%.6e\n", C, name[s], ok ? mx : NAN, ok ? ph : NAN, ok ? sqrt(e2 / n2) : NAN);
            if (ok && (C == 0.5 || C == 2 || C == 5))
                for (int j = 0; j < N; j++) fprintf(fp, "%g,%s,%.5f,%.6e\n", C, name[s], j * dx, u[j]);
        }
        printf("\n");
    }
    fclose(f); fclose(fp);

    /* B. diffusion of a step, large time step */
    for (int j = 0; j < N; j++) u0[j] = (j * dx > 0.25 && j * dx < 0.75) ? 1.0 : 0.0;
    dft(u0, U0);
    double r = 20, dt = r * dx * dx / nu;                 /* nu dt / dx^2 = 20 */
    snprintf(path, sizeof path, "%s/diffusion_profiles.csv", dir); fp = fopen(path, "w");
    fprintf(fp, "steps,scheme,x,u\n");
    printf("\nB. diffusion of a step, nu dt / dx^2 = %g: max error vs exact | total variation (exact in brackets)\n", r);
    snprintf(path, sizeof path, "%s/diffusion.csv", dir); FILE *fd = fopen(path, "w");
    fprintf(fd, "steps,scheme,max_error,total_variation,tv_exact\n");
    int nsts[] = {1, 2, 5, 20};
    for (int q = 0; q < 4; q++) {
        int nst = nsts[q];
        double ue[N]; double complex Ue[N];
        for (int k = 0; k < N; k++) Ue[k] = U0[k] * cexp(-nu * 4 * M_PI * M_PI * wavenumber(k) * wavenumber(k) * dt * nst);
        idft(Ue, ue);
        double tve = 0; for (int j = 0; j < N; j++) tve += fabs(ue[(j + 1) % N] - ue[j]);
        for (int j = 0; j < N; j++) fprintf(fp, "%d,exact,%.5f,%.6e\n", nst, j * dx, ue[j]);
        printf("  after %2d steps (exact TV %.3f):", nst, tve);
        for (int s = EI; s < NS; s++) {
            run(s, U0, nst, lam_dif, dt, u);
            double em = 0, tv = 0;
            for (int j = 0; j < N; j++) { if (fabs(u[j] - ue[j]) > em) em = fabs(u[j] - ue[j]); tv += fabs(u[(j + 1) % N] - u[j]); }
            printf("  %s %.3f | %.3f", name[s], em, tv);
            fprintf(fd, "%d,%s,%.6e,%.6e,%.6e\n", nst, name[s], em, tv, tve);
            for (int j = 0; j < N; j++) fprintf(fp, "%d,%s,%.5f,%.6e\n", nst, name[s], j * dx, u[j]);
        }
        printf("\n");
    }
    fclose(fd);
    fclose(fp);
    return 0;
}
