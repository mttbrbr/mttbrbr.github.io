/**
 * Mesh non-orthogonality and the non-orthogonal correction.
 *
 * Poisson equation -lap(u) = f on a periodic mesh of parallelograms:
 *   x = xi + t * eta, y = eta, t = tan(theta), N x N cells of size h.
 * Every face is non-orthogonal by exactly theta (the angle between the
 * face normal S and the vector d joining the two cell centres), the number
 * checkMesh reports as "non-orthogonality".
 * Manufactured solution u = sin(2 pi xi) cos(2 pi eta), periodic in xi and eta.
 *
 * Face gradient (OpenFOAM snGrad schemes):
 *   orthogonal  : (u_N - u_P) |S| / |d|               (ignores the angle)
 *   uncorrected : (u_N - u_P) |Delta| / |d|, Delta = d (S.S)/(d.S)   (over-relaxed, implicit part only)
 *   corrected   : uncorrected + k . grad(u)_f,  k = S - Delta    (explicit correction,
 *                 updated nCorr times = nNonOrthogonalCorrectors)
 *
 * Build: gcc -O3 -march=native nonortho.c -o nonortho -lm     Run: ./nonortho outdir
 */
#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int N; static double h, t;
#define C(i, j) ((((j) + N) % N) * N + (((i) + N) % N))

static double uex(double xi, double eta) { return sin(2 * M_PI * xi) * cos(2 * M_PI * eta); }
static double fsrc(double xi, double eta) {   /* f = -lap(u) in physical coordinates */
    double a = 2 * M_PI * xi, b = 2 * M_PI * eta;
    return 4 * M_PI * M_PI * (2 + t * t) * sin(a) * cos(b) - 8 * M_PI * M_PI * t * cos(a) * sin(b);
}

/* geometry of the two face families (identical for every cell on this mesh)
   east face (between (i,j) and (i+1,j)): S_e = h (1, -t),  d_e = h (1, 0)
   north face (between (i,j) and (i,j+1)): S_n = h (0, 1),  d_n = h (t, 1)   */
static double Se[2], Sn[2], de[2], dn[2], ke[2], kn[2], cE, cN;

static void geometry(int scheme) {
    Se[0] = h; Se[1] = -t * h; de[0] = h; de[1] = 0;
    Sn[0] = 0; Sn[1] = h;      dn[0] = t * h; dn[1] = h;
    double SSe = Se[0] * Se[0] + Se[1] * Se[1], dSe = de[0] * Se[0] + de[1] * Se[1];
    double SSn = Sn[0] * Sn[0] + Sn[1] * Sn[1], dSn = dn[0] * Sn[0] + dn[1] * Sn[1];
    double lde = hypot(de[0], de[1]), ldn = hypot(dn[0], dn[1]);
    if (scheme == 0) { cE = sqrt(SSe) / lde; cN = sqrt(SSn) / ldn; ke[0] = ke[1] = kn[0] = kn[1] = 0; }   /* orthogonal */
    else {
        cE = SSe / dSe; cN = SSn / dSn;              /* |Delta| / |d| */
        ke[0] = Se[0] - de[0] * SSe / dSe; ke[1] = Se[1] - de[1] * SSe / dSe;
        kn[0] = Sn[0] - dn[0] * SSn / dSn; kn[1] = Sn[1] - dn[1] * SSn / dSn;
    }
}

/* A u = sum over faces of c (u_P - u_N): symmetric positive semi-definite */
static void amul(const double *u, double *y) {
    for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
        int c = C(i, j);
        y[c] = cE * (2 * u[c] - u[C(i + 1, j)] - u[C(i - 1, j)]) + cN * (2 * u[c] - u[C(i, j + 1)] - u[C(i, j - 1)]);
    }
}
static void cg(double *x, const double *b) {
    int n = N * N;
    double *r = malloc(n * 8), *p = malloc(n * 8), *q = malloc(n * 8), mb = 0;
    for (int k = 0; k < n; k++) mb += b[k];
    mb /= n;
    amul(x, q);
    double rr = 0, bb = 0;
    for (int k = 0; k < n; k++) { r[k] = b[k] - mb - q[k]; p[k] = r[k]; rr += r[k] * r[k]; bb += b[k] * b[k]; }
    for (int it = 0; it < 20 * N && rr > 1e-28 * bb; it++) {
        amul(p, q);
        double pq = 0; for (int k = 0; k < n; k++) pq += p[k] * q[k];
        double al = rr / pq, rn = 0;
        for (int k = 0; k < n; k++) { x[k] += al * p[k]; r[k] -= al * q[k]; rn += r[k] * r[k]; }
        double be = rn / rr; rr = rn;
        for (int k = 0; k < n; k++) p[k] = r[k] + be * p[k];
    }
    double m = 0; for (int k = 0; k < n; k++) m += x[k];
    m /= n;
    for (int k = 0; k < n; k++) x[k] -= m;
    free(r); free(p); free(q);
}
/* Gauss gradient with linear face interpolation; faces: e, w, n, s with outward S */
static void gradient(const double *u, double *gx, double *gy) {
    double V = h * h;
    for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
        int c = C(i, j);
        double ue_ = 0.5 * (u[c] + u[C(i + 1, j)]), uw = 0.5 * (u[c] + u[C(i - 1, j)]);
        double un = 0.5 * (u[c] + u[C(i, j + 1)]), us = 0.5 * (u[c] + u[C(i, j - 1)]);
        gx[c] = ((ue_ - uw) * Se[0] + (un - us) * Sn[0]) / V;
        gy[c] = ((ue_ - uw) * Se[1] + (un - us) * Sn[1]) / V;
    }
}
static double l2err(const double *u) {
    double e2 = 0, n2 = 0, me = 0, mx = 0;
    for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) { me += u[C(i, j)]; mx += uex((i + .5) * h, (j + .5) * h); }
    me /= N * N; mx /= N * N;
    for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
        double ex = uex((i + .5) * h, (j + .5) * h) - mx, d = u[C(i, j)] - me - ex;
        e2 += d * d; n2 += ex * ex;
    }
    return sqrt(e2 / n2);
}

/* scheme: 0 orthogonal, 1 uncorrected, 2 corrected with nCorr correctors (-1 = until converged)
   returns relative L2 error; *used = correctors actually performed */
static double solve(int n, double theta_deg, int scheme, int nCorr, int *used, FILE *hist) {
    N = n; h = 1.0 / N; t = tan(theta_deg * M_PI / 180);
    geometry(scheme);
    int nc = N * N;
    double *u = calloc(nc, 8), *b = malloc(nc * 8), *gx = calloc(nc, 8), *gy = calloc(nc, 8), *uold = malloc(nc * 8);
    int maxc = scheme == 2 ? (nCorr < 0 ? 200 : nCorr) : 0, k;
    for (k = 0; k <= maxc; k++) {
        /* rhs: source + explicit non-orthogonal correction from the current solution */
        if (scheme == 2) gradient(u, gx, gy);
        for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
            int c = C(i, j);
            double s = fsrc((i + .5) * h, (j + .5) * h) * h * h;
            if (scheme == 2) {   /* sum over faces of k . grad_f (outward); w and s faces use -k */
                double ge = (ke[0] * 0.5 * (gx[c] + gx[C(i + 1, j)]) + ke[1] * 0.5 * (gy[c] + gy[C(i + 1, j)]));
                double gw = (ke[0] * 0.5 * (gx[c] + gx[C(i - 1, j)]) + ke[1] * 0.5 * (gy[c] + gy[C(i - 1, j)]));
                double gn = (kn[0] * 0.5 * (gx[c] + gx[C(i, j + 1)]) + kn[1] * 0.5 * (gy[c] + gy[C(i, j + 1)]));
                double gs = (kn[0] * 0.5 * (gx[c] + gx[C(i, j - 1)]) + kn[1] * 0.5 * (gy[c] + gy[C(i, j - 1)]));
                s += (ge - gw) + (gn - gs);
            }
            b[c] = s;
        }
        memcpy(uold, u, nc * 8);
        cg(u, b);
        double du = 0, nu_ = 0;
        for (int c = 0; c < nc; c++) { du += (u[c] - uold[c]) * (u[c] - uold[c]); nu_ += u[c] * u[c]; }
        if (hist) fprintf(hist, "%g,%d,%d,%.6e,%.6e\n", theta_deg, n, k, sqrt(du / nu_), l2err(u));
        if (!isfinite(du) || sqrt(du / nu_) > 1e3) { k = -1; break; }
        if (nCorr < 0 && k > 0 && sqrt(du / nu_) < 1e-10) break;
    }
    if (used) *used = k;
    double e = (k == -1) ? NAN : l2err(u);
    free(u); free(b); free(gx); free(gy); free(uold);
    return e;
}

int main(int argc, char **argv) {
    const char *dir = argc > 1 ? argv[1] : ".";
    char path[512]; FILE *f, *hist;
    snprintf(path, sizeof path, "%s/correctors.csv", dir); hist = fopen(path, "w");
    fprintf(hist, "theta,N,corrector,change,error\n");

    /* A. error vs angle, N = 64 */
    snprintf(path, sizeof path, "%s/angle.csv", dir); f = fopen(path, "w");
    fprintf(f, "theta,orthogonal,uncorrected,corrected_1,corrected_2,corrected_conv,correctors_needed\n");
    printf("A. relative L2 error vs non-orthogonality, 64 x 64 cells\n");
    printf("%7s %12s %12s %12s %12s %12s %8s\n", "theta", "orthogonal", "uncorrected", "1 corrector", "2 correctors", "converged", "needed");
    double th[] = {0, 10, 20, 30, 40, 50, 60, 70, 75, 80, 85};
    for (int a = 0; a < 11; a++) {
        int used;
        double eo = solve(64, th[a], 0, 0, NULL, NULL), eu = solve(64, th[a], 1, 0, NULL, NULL);
        double e1 = solve(64, th[a], 2, 1, NULL, NULL), e2 = solve(64, th[a], 2, 2, NULL, NULL);
        double ec = solve(64, th[a], 2, -1, &used, hist);
        printf("%7.0f %12.3e %12.3e %12.3e %12.3e %12.3e %8d\n", th[a], eo, eu, e1, e2, ec, used);
        fprintf(f, "%g,%.6e,%.6e,%.6e,%.6e,%.6e,%d\n", th[a], eo, eu, e1, e2, ec, used);
    }
    fclose(f);

    /* B. error vs mesh size at 0, 30 and 60 degrees */
    snprintf(path, sizeof path, "%s/grid.csv", dir); f = fopen(path, "w");
    fprintf(f, "theta,N,orthogonal,uncorrected,corrected_conv\n");
    printf("\nB. relative L2 error vs mesh size\n%7s %6s %12s %12s %12s\n", "theta", "N", "orthogonal", "uncorrected", "corrected");
    double tb[] = {0, 30, 60};
    for (int a = 0; a < 3; a++) for (int n = 16; n <= 256; n *= 2) {
        int used;
        double eo = solve(n, tb[a], 0, 0, NULL, NULL), eu = solve(n, tb[a], 1, 0, NULL, NULL), ec = solve(n, tb[a], 2, -1, &used, NULL);
        printf("%7.0f %6d %12.3e %12.3e %12.3e\n", tb[a], n, eo, eu, ec);
        fprintf(f, "%g,%d,%.6e,%.6e,%.6e\n", tb[a], n, eo, eu, ec);
    }
    fclose(f); fclose(hist);
    return 0;
}
