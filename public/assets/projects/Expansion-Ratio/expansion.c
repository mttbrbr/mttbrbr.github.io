/**
 * Mesh expansion ratio: when do changes in cell size matter?
 *
 * 1D cell-centred finite volumes on non-uniform meshes. The mesh goes from a
 * fine region (cell size h) to a coarse one (cell size R h, R = 2 or 4) either
 * abruptly (a single jump, like a 2:1 refinement interface) or gradually, with
 * a constant expansion ratio r between neighbouring cells.
 *
 *  A. Diffusion: -u'' = f on [0,1], exact u = sin(3 pi x) + x, Dirichlet ends.
 *     Global error and local truncation error vs N for each mesh type.
 *  B. Waves: u_t + u_x = 0, periodic, central (linear-interpolation) fluxes,
 *     RK4. A wave packet crosses the transition from fine to coarse; part of
 *     it is reflected as spurious waves. Measured for several resolutions
 *     (fine cells per wavelength) after one crossing.
 *
 * Build: gcc -O3 -march=native expansion.c -o expansion -lm   Run: ./expansion outdir
 */
#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* build cell faces on [0, L): fine region, transition, coarse region.
   r <= 0: abrupt jump; r > 1: geometric grading with ratio r until size R h.
   Returns the number of cells; faces in xf[0..n]. */
static int build_mesh(double L, double h, double R, double r, double xjump, double *xf, int maxn) {
    int n = 0; double x = 0, sz = h;
    xf[0] = 0;
    while (x < xjump - 1e-12) { x += h; xf[++n] = x; }
    if (r > 1) { while (sz * r < R * h * (1 + 1e-9)) { sz *= r; x += sz; xf[++n] = x; } }
    sz = R * h;
    while (x + sz <= L + 1e-9 && n < maxn - 1) { x += sz; xf[++n] = x; }
    return n;
}

/* ---------- A. diffusion ---------- */
static double uex(double x) { return sin(3 * M_PI * x) + x; }
static void tridiag(int n, double *a, double *b, double *c, double *d) {   /* in place, solution in d */
    for (int i = 1; i < n; i++) { double m = a[i] / b[i - 1]; b[i] -= m * c[i - 1]; d[i] -= m * d[i - 1]; }
    d[n - 1] /= b[n - 1];
    for (int i = n - 2; i >= 0; i--) d[i] = (d[i] - c[i] * d[i + 1]) / b[i];
}
/* returns L2 error; *trunc = max local truncation error (residual of the exact solution / cell size) */
static double diffusion(int n, const double *xf, double *trunc, double *uerrmax) {
    double *xc = malloc(n * 8), *a = calloc(n, 8), *b = calloc(n, 8), *c = calloc(n, 8), *d = calloc(n, 8);
    for (int i = 0; i < n; i++) xc[i] = 0.5 * (xf[i] + xf[i + 1]);
    for (int i = 0; i < n; i++) {
        double V = xf[i + 1] - xf[i];
        /* exact cell integral of f */
        d[i] = (-3 * M_PI) * (cos(3 * M_PI * xf[i + 1]) - cos(3 * M_PI * xf[i])) / 1.0 * (1.0);
        d[i] = 3 * M_PI * (cos(3 * M_PI * xf[i]) - cos(3 * M_PI * xf[i + 1]));   /* int f dx */
        (void)V;
        double dw = i > 0 ? xc[i] - xc[i - 1] : xc[0] - xf[0];
        double de = i < n - 1 ? xc[i + 1] - xc[i] : xf[n] - xc[n - 1];
        double aw = 1 / dw, ae = 1 / de;
        b[i] = aw + ae;
        if (i > 0) a[i] = -aw; else d[i] += aw * uex(xf[0]);
        if (i < n - 1) c[i] = -ae; else d[i] += ae * uex(xf[n]);
    }
    /* truncation error: plug the exact solution into the discrete equations */
    double tmax = 0;
    for (int i = 0; i < n; i++) {
        double r = b[i] * uex(xc[i]) - d[i] + (i > 0 ? a[i] * uex(xc[i - 1]) : 0) + (i < n - 1 ? c[i] * uex(xc[i + 1]) : 0);
        double tr = fabs(r) / (xf[i + 1] - xf[i]);
        if (tr > tmax) tmax = tr;
    }
    tridiag(n, a, b, c, d);
    double e2 = 0, emax = 0;
    for (int i = 0; i < n; i++) { double e = fabs(d[i] - uex(xc[i])); e2 += e * e * (xf[i + 1] - xf[i]); if (e > emax) emax = e; }
    *trunc = tmax; *uerrmax = emax;
    free(xc); free(a); free(b); free(c); free(d);
    return sqrt(e2);
}

/* ---------- B. waves ---------- */
static int NW; static double *xfW, *VW;
static void rhs_adv(const double *u, double *r) {
    /* flux at face i+1/2 (between cell i and i+1, periodic) by linear interpolation at the face */
    for (int i = 0; i < NW; i++) {
        int ip = (i + 1) % NW;
        double xci = 0.5 * (xfW[i] + xfW[i + 1]);
        double xcp = ip ? 0.5 * (xfW[ip] + xfW[ip + 1]) : 0.5 * (xfW[0] + xfW[1]) + xfW[NW];
        double xface = xfW[i + 1], w = (xcp - xface) / (xcp - xci);
        double F = w * u[i] + (1 - w) * u[ip];
        r[i] -= F / VW[i]; r[ip] += F / VW[ip];
    }
}
static double packet(double x, double x0, double lam) { return exp(-0.5 * pow((x - x0) / (2 * lam), 2)) * sin(2 * M_PI * (x - x0) / lam); }

int main(int argc, char **argv) {
    const char *dir = argc > 1 ? argv[1] : ".";
    char path[512];
    int maxn = 400000; double *xf = malloc((maxn + 1) * 8);
    struct { const char *name; double R, r; } M[] = {
        {"uniform", 1, 0}, {"jump 2:1", 2, 0}, {"graded r=1.2 to 2x", 2, 1.2}, {"graded r=1.1 to 2x", 2, 1.1},
        {"jump 4:1", 4, 0}, {"graded r=1.2 to 4x", 4, 1.2}};
    int nm = 6;

    /* A */
    snprintf(path, sizeof path, "%s/diffusion.csv", dir); FILE *f = fopen(path, "w");
    fprintf(f, "mesh,cells,h_fine,l2_error,max_error,max_truncation\n");
    printf("A. diffusion: L2 error / max truncation error\n");
    for (int m = 0; m < nm; m++) {
        printf("  %-22s", M[m].name);
        for (int k = 0; k < 6; k++) {
            double h = 1.0 / (40 << k);
            int n = build_mesh(1.0, h, M[m].R, M[m].r, 0.4, xf, maxn);
            xf[n] = 1.0;                                           /* close the domain exactly */
            double tr, em, e = diffusion(n, xf, &tr, &em);
            fprintf(f, "%s,%d,%.6e,%.6e,%.6e,%.6e\n", M[m].name, n, h, e, em, tr);
            printf("  %5d: %.1e/%.1e", n, e, tr);
        }
        printf("\n");
    }
    fclose(f);

    /* B */
    snprintf(path, sizeof path, "%s/waves.csv", dir); f = fopen(path, "w");
    fprintf(f, "mesh,cells_per_wavelength_fine,reflected_fraction,transmitted_fraction,energy_kept\n");
    snprintf(path, sizeof path, "%s/wave_profiles.csv", dir); FILE *fp = fopen(path, "w");
    fprintf(fp, "mesh,ppw,x,u\n");
    printf("\nB. waves: fraction of the packet energy reflected backwards after crossing the transition\n");
    int ppws[] = {16, 32};
    for (int q = 0; q < 2; q++)
    for (int m = 0; m < nm; m++) {
        int ppw = ppws[q];
        double h = 1.0 / 400, lam = ppw * h;                     /* ppw fine cells per wavelength */
        double L = 6.0, xj = 2.0;
        NW = build_mesh(L, h, M[m].R, M[m].r, xj, xf, maxn);
        /* make the domain periodic: coarse region ends at L; go back to fine at the end abruptly is avoided:
           mirror the transition by appending a reversed copy */
        int n1 = NW; double Lr = xf[n1];
        xfW = malloc((2 * n1 + 1) * 8);
        for (int i = 0; i <= n1; i++) xfW[i] = xf[i];
        for (int i = 1; i <= n1; i++) xfW[n1 + i] = Lr + (Lr - xf[n1 - i]);
        NW = 2 * n1; double Lt = xfW[NW];
        VW = malloc(NW * 8);
        for (int i = 0; i < NW; i++) VW[i] = xfW[i + 1] - xfW[i];
        double *u = malloc(NW * 8), *k1 = malloc(NW * 8), *k2 = malloc(NW * 8), *k3 = malloc(NW * 8), *k4 = malloc(NW * 8), *tmp = malloc(NW * 8);
        double x0 = xj - 8 * lam;                                  /* packet starts in the fine region */
        for (int i = 0; i < NW; i++) u[i] = packet(0.5 * (xfW[i] + xfW[i + 1]), x0, lam);
        double e0 = 0; for (int i = 0; i < NW; i++) e0 += u[i] * u[i] * VW[i];
        double dt = 0.4 * h, T = 16 * lam;                         /* the packet ends 8 wavelengths past the jump */
        int steps = (int)(T / dt);
        for (int s = 0; s < steps; s++) {
            #define RHS(in, out) do { memset(out, 0, NW * 8); rhs_adv(in, out); } while (0)
            RHS(u, k1); for (int i = 0; i < NW; i++) tmp[i] = u[i] + 0.5 * dt * k1[i];
            RHS(tmp, k2); for (int i = 0; i < NW; i++) tmp[i] = u[i] + 0.5 * dt * k2[i];
            RHS(tmp, k3); for (int i = 0; i < NW; i++) tmp[i] = u[i] + dt * k3[i];
            RHS(tmp, k4); for (int i = 0; i < NW; i++) u[i] += dt / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
        }
        /* reflected energy: everything left of the transition (x < 0.8) after the packet has moved on */
        double er = 0, et = 0;
        for (int i = 0; i < NW; i++) {
            double xc = 0.5 * (xfW[i] + xfW[i + 1]);
            if (xc < xj - 4 * lam) er += u[i] * u[i] * VW[i]; else et += u[i] * u[i] * VW[i];
            if (xc < xj + 12 * lam && xc > xj - 14 * lam) fprintf(fp, "%s,%d,%.5f,%.6e\n", M[m].name, ppw, xc, u[i]);
        }
        printf("  %2d cells/wavelength  %-22s reflected %.2e  transmitted %.3f\n", ppw, M[m].name, er / e0, et / e0);
        fprintf(f, "%s,%d,%.6e,%.6e,%.6e\n", M[m].name, ppw, er / e0, et / e0, (er + et) / e0);
        (void)Lt;
        free(u); free(k1); free(k2); free(k3); free(k4); free(tmp); free(xfW); free(VW);
    }
    fclose(f); fclose(fp);
    return 0;
}
