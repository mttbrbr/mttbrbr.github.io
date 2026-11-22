/**
 * The core of an OpenFOAM pressure solve, in ~250 lines of C.
 *
 * - Mesh: n x n x n hexahedral cells, but stored like an unstructured mesh:
 *   a list of internal faces with owner/neighbour cell indices.
 * - Matrix: OpenFOAM's LDU format: one diagonal coefficient per cell, one
 *   off-diagonal coefficient per face (symmetric: upper == lower).
 *   Discretized Laplacian (pressure-like equation), Dirichlet walls.
 * - Solver: preconditioned conjugate gradient with diagonal (Jacobi)
 *   preconditioner: Amul over faces, dot products, axpy. This is where an
 *   incompressible OpenFOAM run spends most of its time.
 *
 * Precision (compile time), named after OpenFOAM's WM_PRECISION_OPTION:
 *   -DPREC_DP   (default) matrix and solver vectors in double
 *   -DPREC_SP   matrix and solver vectors in float
 *   -DPREC_SPDP matrix stored in float, solver vectors in double
 *
 * Run: ./ldu_cg n mode
 *   mode = time  : ~1e8 cell-iterations (restarted every 50), prints ns per cell per iteration
 *   mode = solve : converge to 1e-12 (or stall), prints residual history
 *   mode = fields: face interpolation + divergence + field update, in storage precision
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#if defined(PREC_SP)
typedef float scalar;  typedef float  solveScalar; static const char *PREC = "SP";
#elif defined(PREC_SPDP)
typedef float scalar;  typedef double solveScalar; static const char *PREC = "SPDP";
#else
typedef double scalar; typedef double solveScalar; static const char *PREC = "DP";
#endif

#ifndef FLAGS
#define FLAGS "?"
#endif

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

static int nCells, nFaces;
static int *own, *nei;
static scalar *diag, *upper;

static double face_coef(int f) { return 1.0 + 0.3 * sin(0.37 * f); }

static void build(int n) {
    nCells = n * n * n;
    nFaces = 3 * n * n * (n - 1);
    own = malloc(nFaces * sizeof(int)); nei = malloc(nFaces * sizeof(int));
    diag = calloc(nCells, sizeof(scalar)); upper = malloc(nFaces * sizeof(scalar));
    int f = 0;
    /* faces ordered by owner, neighbour ascending: OpenFOAM's "upper triangular" order */
    for (int k = 0; k < n; k++) for (int j = 0; j < n; j++) for (int i = 0; i < n; i++) {
        int c = (k * n + j) * n + i;
        if (i + 1 < n) { own[f] = c; nei[f] = c + 1;     f++; }
        if (j + 1 < n) { own[f] = c; nei[f] = c + n;     f++; }
        if (k + 1 < n) { own[f] = c; nei[f] = c + n * n; f++; }
    }
    /* face coefficients |S_f| / d_f of a slightly non-uniform mesh: not round numbers */
    double *dd = calloc(nCells, sizeof(double));
    for (f = 0; f < nFaces; f++) {
        double a = face_coef(f);
        upper[f] = (scalar)(-a);
        dd[own[f]] += a; dd[nei[f]] += a;
    }
    /* diagonal: sum of neighbour coefficients + 2 for every wall face (Dirichlet) */
    for (int c = 0; c < nCells; c++) {
        int i = c % n, j = (c / n) % n, k = c / (n * n);
        int walls = (i == 0) + (i == n - 1) + (j == 0) + (j == n - 1) + (k == 0) + (k == n - 1);
        diag[c] = (scalar)(dd[c] + 2.0 * walls);
    }
    free(dd);
}

/* y = A x, OpenFOAM's lduMatrix::Amul */
static void Amul(solveScalar *restrict y, const solveScalar *restrict x) {
    for (int c = 0; c < nCells; c++) y[c] = diag[c] * x[c];
    for (int f = 0; f < nFaces; f++) {
        y[own[f]] += upper[f] * x[nei[f]];
        y[nei[f]] += upper[f] * x[own[f]];
    }
}
static solveScalar dot(const solveScalar *a, const solveScalar *b) {
    solveScalar s = 0; for (int c = 0; c < nCells; c++) s += a[c] * b[c]; return s;
}

/* exact residual in double with the double-precision matrix (the "true" problem) */
static double true_residual(const solveScalar *x, const double *b) {
    double r2 = 0, b2 = 0; int n = (int)lround(cbrt(nCells));
    double *y = calloc(nCells, sizeof(double));
    for (int c = 0; c < nCells; c++) {
        int i = c % n, j = (c / n) % n, k = c / (n * n);
        int walls = (i == 0) + (i == n - 1) + (j == 0) + (j == n - 1) + (k == 0) + (k == n - 1);
        y[c] = 2.0 * walls * (double)x[c];
    }
    for (int f = 0; f < nFaces; f++) {                  /* matrix rebuilt in double */
        double a = face_coef(f);
        y[own[f]] += a * ((double)x[own[f]] - (double)x[nei[f]]);
        y[nei[f]] += a * ((double)x[nei[f]] - (double)x[own[f]]);
    }
    for (int c = 0; c < nCells; c++) { double r = b[c] - y[c]; r2 += r * r; b2 += b[c] * b[c]; }
    free(y); return sqrt(r2 / b2);
}

/* "Field operations": the rest of a finite-volume iteration, in storage precision
   (scalar): interpolate a vector field to the faces, compute the face flux
   phi = Sf . U_f, accumulate its divergence in the cells, update the field.
   Like fvc::interpolate + fvc::div + a field update in OpenFOAM. */
static int fields(int n) {
    int nc = nCells, nf = nFaces;
    scalar *U = malloc(3 * (size_t)nc * sizeof(scalar)), *Sf = malloc(3 * (size_t)nf * sizeof(scalar));
    scalar *w = malloc(nf * sizeof(scalar)), *phi = malloc(nf * sizeof(scalar)), *dv = malloc(nc * sizeof(scalar));
    for (int c = 0; c < nc; c++) for (int d = 0; d < 3; d++) U[3 * c + d] = (scalar)(1.0 + 0.1 * sin(0.01 * c + d));
    for (int f = 0; f < nf; f++) {
        int dir = (nei[f] - own[f] == 1) ? 0 : (nei[f] - own[f] == n ? 1 : 2);
        for (int d = 0; d < 3; d++) Sf[3 * f + d] = (scalar)(d == dir ? 1.0 : 0.0);
        w[f] = (scalar)0.5;
    }
    int maxit = (int)(1e8 / nc); if (maxit < 20) maxit = 20; if (maxit > 2000) maxit = 2000;
    double t0 = now();
    for (int it = 0; it < maxit; it++) {
        for (int f = 0; f < nf; f++) {
            const scalar *uo = U + 3 * (size_t)own[f], *un = U + 3 * (size_t)nei[f], *s = Sf + 3 * (size_t)f;
            scalar wf = w[f];
            phi[f] = s[0] * (wf * uo[0] + (1 - wf) * un[0]) + s[1] * (wf * uo[1] + (1 - wf) * un[1])
                   + s[2] * (wf * uo[2] + (1 - wf) * un[2]);
        }
        for (int c = 0; c < nc; c++) dv[c] = 0;
        for (int f = 0; f < nf; f++) { dv[own[f]] += phi[f]; dv[nei[f]] -= phi[f]; }
        for (int c = 0; c < nc; c++) for (int d = 0; d < 3; d++) U[3 * c + d] -= (scalar)1e-6 * dv[c];
    }
    double el = now() - t0, chk = 0;
    for (int c = 0; c < nc; c++) chk += U[3 * c];
    double mb = (nc * 4.0 * sizeof(scalar) + nf * (5.0 * sizeof(scalar) + 8.0)) / 1e6;
    printf("fields,%s,\"%s\",%d,%d,%.2f,%.4f,%.10e\n", PREC, FLAGS, n, nc, mb, el * 1e9 / ((double)maxit * nc), chk);
    return 0;
}

int main(int argc, char **argv) {
    int n = argc > 1 ? atoi(argv[1]) : 64;
    const char *mode = argc > 2 ? argv[2] : "time";
    build(n);
    if (!strcmp(mode, "fields")) return fields(n);
    solveScalar *x = calloc(nCells, sizeof *x), *r = malloc(nCells * sizeof *r), *z = malloc(nCells * sizeof *z),
                *p = malloc(nCells * sizeof *p), *q = malloc(nCells * sizeof *q);
    double *bd = malloc(nCells * sizeof(double));
    for (int c = 0; c < nCells; c++) bd[c] = 1.0 + 0.5 * sin(0.001 * c);   /* source term */
    for (int c = 0; c < nCells; c++) r[c] = (solveScalar)bd[c];
    int timing = !strcmp(mode, "time");
    int maxit = 5000;
    if (timing) { maxit = (int)(1e8 / nCells); maxit = (maxit + 49) / 50 * 50; if (maxit < 50) maxit = 50; if (maxit > 2000) maxit = 2000; }
    double tol = 1e-12;

    /* PCG, diagonal preconditioner */
    solveScalar rz = 0, r0 = sqrt(dot(r, r));
    for (int c = 0; c < nCells; c++) { z[c] = r[c] / diag[c]; p[c] = z[c]; }
    rz = dot(r, z);
    double t0 = now(); int it;
    for (it = 1; it <= maxit; it++) {
        if (timing && it % 50 == 0) {          /* restart: keep timing away from converged (tiny) numbers */
            for (int c = 0; c < nCells; c++) { x[c] = 0; r[c] = (solveScalar)bd[c]; z[c] = r[c] / diag[c]; p[c] = z[c]; }
            rz = dot(r, z);
        }
        Amul(q, p);
        solveScalar al = rz / dot(p, q);
        for (int c = 0; c < nCells; c++) { x[c] += al * p[c]; r[c] -= al * q[c]; }
        solveScalar rn = sqrt(dot(r, r));
        if (!timing && (it % 10 == 0 || it == 1))
            printf("hist,%s,%d,%.6e,%.6e\n", PREC, it, (double)(rn / r0), true_residual(x, bd));
        if (!timing && rn / r0 < tol) break;
        for (int c = 0; c < nCells; c++) z[c] = r[c] / diag[c];
        solveScalar rz2 = dot(r, z), be = rz2 / rz; rz = rz2;
        for (int c = 0; c < nCells; c++) p[c] = z[c] + be * p[c];
    }
    double el = now() - t0;
    if (it > maxit) it = maxit;
    double chk = 0; for (int c = 0; c < nCells; c++) chk += x[c];
    double mb = (nCells * (sizeof(scalar) + 5.0 * sizeof(solveScalar)) + nFaces * (sizeof(scalar) + 8.0)) / 1e6;
    if (timing)
        printf("time,%s,\"%s\",%d,%d,%.2f,%.4f,%.10e\n", PREC, FLAGS, n, nCells, mb, el * 1e9 / ((double)it * nCells), chk);
    else
    {
        printf("final,%s,%d,%d,%.6e,%.10e\n", PREC, n, it, true_residual(x, bd), chk);
        char name[64]; snprintf(name, sizeof name, "x_%s_%d.bin", PREC, n);   /* solution, as double */
        FILE *o = fopen(name, "wb");
        for (int c = 0; c < nCells; c++) { double v = x[c]; fwrite(&v, sizeof v, 1, o); }
        fclose(o);
    }
    return 0;
}
