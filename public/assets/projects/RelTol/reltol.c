/**
 * How accurately should the inner linear systems be solved? (fvSolution relTol)
 * Derived from the segregated-vs-coupled cavity code: same discretization,
 * same SIMPLE and coupled solvers, with the inner tolerances as parameters.
 *
 * Original description:
 * Segregated (SIMPLE) vs coupled pressure-velocity solution of the steady 2D
 * lid-driven cavity, on the same discretization.
 *
 *  - Staggered (MAC) grid, N x N cells, unit square, lid u = 1 at the top.
 *  - Finite volumes, first-order upwind convection, central diffusion.
 *  - Non-linearity handled the same way in both methods: Picard iteration
 *    (face mass fluxes taken from the previous outer iteration).
 *
 *  SIMPLE   : per outer iteration solve u*, then v*, then a pressure
 *             correction p', correct u, v, p. Needs under-relaxation (alpha_u, alpha_p).
 *  COUPLED  : per outer iteration assemble ONE linear system with u, v and p
 *             together (saddle-point matrix) and solve it.
 *
 *  Every linear system is solved exactly with the same banded LU (LAPACK
 *  dgbtrf/dgbtrs), so the only difference is how pressure and velocity are
 *  coupled. Convergence is measured with the SAME residual for both: the full
 *  non-linear residual of momentum + continuity.
 *
 * Build: gcc -O3 -march=native cavity.c -o cavity -llapack -lm
 * Run:   ./cavity outdir
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* ---------------- LAPACK banded LU ---------------- */
extern void dgbtrf_(int *m, int *n, int *kl, int *ku, double *ab, int *ldab, int *ipiv, int *info);
extern void dgbtrs_(char *tr, int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab, int *ipiv,
                    double *b, int *ldb, int *info, size_t len);

typedef struct { int n, kl, ku, ldab; double *ab; int *ipiv; } Band;
static Band band_new(int n, int kl, int ku) {
    Band B = {n, kl, ku, 2 * kl + ku + 1, NULL, NULL};
    B.ab = calloc((size_t)B.ldab * n, sizeof(double));
    B.ipiv = malloc(n * sizeof(int));
    return B;
}
static void band_clear(Band *B) { memset(B->ab, 0, (size_t)B->ldab * B->n * sizeof(double)); }
static inline void band_add(Band *B, int r, int c, double v) { B->ab[(B->kl + B->ku + r - c) + (size_t)c * B->ldab] += v; }
static void band_solve(Band *B, double *rhs) {
    int info, one = 1; char t = 'N';
    dgbtrf_(&B->n, &B->n, &B->kl, &B->ku, B->ab, &B->ldab, B->ipiv, &info);
    if (info) { fprintf(stderr, "dgbtrf failed (%d)\n", info); exit(1); }
    dgbtrs_(&t, &B->n, &B->kl, &B->ku, &one, B->ab, &B->ldab, B->ipiv, rhs, &B->n, &info, 1);
}
static void band_free(Band *B) { free(B->ab); free(B->ipiv); }
static double band_mb(const Band *B) { return (double)B->ldab * B->n * 8 / 1e6; }

static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

/* ---------------- grid and fields ---------------- */
static int N; static double h, nu;
/* u[i][j]: vertical face i (x = i h), row j, i = 0..N (0 and N are walls)
   v[i][j]: horizontal face j (y = j h), column i, j = 0..N
   p[i][j]: cell centre */
#define U(i, j) u[(i) * N + (j)]
#define V(i, j) v[(i) * (N + 1) + (j)]
#define P(i, j) p[(i) * N + (j)]

/* momentum coefficients (Picard: fluxes from current fields) */
typedef struct { double aP, aE, aW, aN, aS, b; } Coef;   /* b: boundary source (lid) */

static Coef coef_u(const double *u, const double *v, int i, int j) {
    Coef c = {0};
    double Fe = 0.5 * h * (U(i, j) + U(i + 1, j)), Fw = 0.5 * h * (U(i - 1, j) + U(i, j));
    double Fn = 0.5 * h * (V(i - 1, j + 1) + V(i, j + 1)), Fs = 0.5 * h * (V(i - 1, j) + V(i, j));
    double D = nu;
    c.aE = D + fmax(-Fe, 0); c.aW = D + fmax(Fw, 0);
    double sum = c.aE + c.aW;
    if (j == N - 1) { sum += 2 * nu; c.b += 2 * nu * 1.0; c.aN = 0; }   /* moving lid, u = 1 */
    else { c.aN = D + fmax(-Fn, 0); sum += c.aN; }
    if (j == 0) { sum += 2 * nu; c.aS = 0; }                         /* bottom wall, u = 0 */
    else { c.aS = D + fmax(Fs, 0); sum += c.aS; }
    if (i + 1 == N) c.aE = 0;                                         /* wall face, u = 0 */
    if (i - 1 == 0) c.aW = 0;
    c.aP = sum + (Fe - Fw + Fn - Fs);
    return c;
}
static Coef coef_v(const double *u, const double *v, int i, int j) {
    Coef c = {0};
    double Fn = 0.5 * h * (V(i, j) + V(i, j + 1)), Fs = 0.5 * h * (V(i, j - 1) + V(i, j));
    double Fe = 0.5 * h * (U(i + 1, j - 1) + U(i + 1, j)), Fw = 0.5 * h * (U(i, j - 1) + U(i, j));
    double D = nu;
    c.aN = D + fmax(-Fn, 0); c.aS = D + fmax(Fs, 0);
    double sum = c.aN + c.aS;
    if (i == N - 1) { sum += 2 * nu; c.aE = 0; } else { c.aE = D + fmax(-Fe, 0); sum += c.aE; }
    if (i == 0) { sum += 2 * nu; c.aW = 0; } else { c.aW = D + fmax(Fw, 0); sum += c.aW; }
    if (j + 1 == N) c.aN = 0;
    if (j - 1 == 0) c.aS = 0;
    c.aP = sum + (Fe - Fw + Fn - Fs);
    return c;
}

/* full non-linear residual: momentum (u, v) + continuity, L2 norm */
static double residual(const double *u, const double *v, const double *p) {
    double r2 = 0;
    for (int i = 1; i < N; i++) for (int j = 0; j < N; j++) {
        Coef c = coef_u(u, v, i, j);
        double r = c.aP * U(i, j) - c.aE * U(i + 1, j) - c.aW * U(i - 1, j)
                 - (j < N - 1 ? c.aN * U(i, j + 1) : 0) - (j > 0 ? c.aS * U(i, j - 1) : 0)
                 - c.b - h * (P(i - 1, j) - P(i, j));
        r2 += r * r;
    }
    for (int i = 0; i < N; i++) for (int j = 1; j < N; j++) {
        Coef c = coef_v(u, v, i, j);
        double r = c.aP * V(i, j) - c.aN * V(i, j + 1) - c.aS * V(i, j - 1)
                 - (i < N - 1 ? c.aE * V(i + 1, j) : 0) - (i > 0 ? c.aW * V(i - 1, j) : 0)
                 - h * (P(i, j - 1) - P(i, j));
        r2 += r * r;
    }
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        double r = h * (U(i + 1, j) - U(i, j) + V(i, j + 1) - V(i, j));
        r2 += r * r;
    }
    return sqrt(r2);
}

typedef struct { int iters, converged; double seconds, mb, umin, ymin; } Result;

static void centreline(const double *u, Result *R) {   /* min u on x = 0.5 (face i = N/2) */
    R->umin = 1e9;
    for (int j = 0; j < N; j++) if (U(N / 2, j) < R->umin) { R->umin = U(N / 2, j); R->ymin = (j + 0.5) * h; }
}

#define TOL 1e-8
#define MAXIT 20000

/* ---------------- SIMPLE ---------------- */
static Result simple(double au, double ap, FILE *hist, const char *tag) {
    Result R = {0};
    double *u = calloc((N + 1) * N, 8), *v = calloc(N * (N + 1), 8), *p = calloc(N * N, 8);
    double *du = calloc((N + 1) * N, 8), *dv = calloc(N * (N + 1), 8), *pc = malloc(N * N * 8);
    int nu_ = (N - 1) * N, nv_ = N * (N - 1), np_ = N * N;
    Band Bu = band_new(nu_, N - 1, N - 1), Bv = band_new(nv_, N, N), Bp = band_new(np_, N, N);
    double *ru = malloc(nu_ * 8), *rv = malloc(nv_ * 8), *rp = malloc(np_ * 8);
    double *u2 = malloc((N + 1) * N * 8), *v2 = malloc(N * (N + 1) * 8);
    R.mb = fmax(fmax(band_mb(&Bu), band_mb(&Bv)), band_mb(&Bp));
    double r0 = residual(u, v, p), t0 = now();
    for (int it = 1; it <= MAXIT; it++) {
        /* u* */
        band_clear(&Bu);
        for (int i = 1; i < N; i++) for (int j = 0; j < N; j++) {
            Coef c = coef_u(u, v, i, j); int k = j * (N - 1) + (i - 1);
            double aPr = c.aP / au;
            band_add(&Bu, k, k, aPr);
            if (i + 1 < N) band_add(&Bu, k, k + 1, -c.aE);
            if (i - 1 > 0) band_add(&Bu, k, k - 1, -c.aW);
            if (j < N - 1) band_add(&Bu, k, k + (N - 1), -c.aN);
            if (j > 0) band_add(&Bu, k, k - (N - 1), -c.aS);
            ru[k] = c.b + h * (P(i - 1, j) - P(i, j)) + (1 - au) * aPr * U(i, j);
            du[i * N + j] = h / aPr;
        }
        /* v* (with the same old fields) */
        band_clear(&Bv);
        for (int i = 0; i < N; i++) for (int j = 1; j < N; j++) {
            Coef c = coef_v(u, v, i, j); int k = (j - 1) * N + i;
            double aPr = c.aP / au;
            band_add(&Bv, k, k, aPr);
            if (i + 1 < N) band_add(&Bv, k, k + 1, -c.aE);
            if (i > 0) band_add(&Bv, k, k - 1, -c.aW);
            if (j + 1 < N) band_add(&Bv, k, k + N, -c.aN);
            if (j - 1 > 0) band_add(&Bv, k, k - N, -c.aS);
            rv[k] = c.b + h * (P(i, j - 1) - P(i, j)) + (1 - au) * aPr * V(i, j);
            dv[i * (N + 1) + j] = h / aPr;
        }
        band_solve(&Bu, ru); band_solve(&Bv, rv);
        memcpy(u2, u, (N + 1) * N * 8); memcpy(v2, v, N * (N + 1) * 8);
        for (int i = 1; i < N; i++) for (int j = 0; j < N; j++) u2[i * N + j] = ru[j * (N - 1) + (i - 1)];
        for (int i = 0; i < N; i++) for (int j = 1; j < N; j++) v2[i * (N + 1) + j] = rv[(j - 1) * N + i];
        /* pressure correction */
        band_clear(&Bp);
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
            int k = j * N + i;
            if (k == 0) { band_add(&Bp, 0, 0, 1.0); rp[0] = 0; continue; }   /* pin p'(0,0) */
            double aE = (i + 1 < N) ? h * du[(i + 1) * N + j] : 0, aW = (i > 0) ? h * du[i * N + j] : 0;
            double aN = (j + 1 < N) ? h * dv[i * (N + 1) + j + 1] : 0, aS = (j > 0) ? h * dv[i * (N + 1) + j] : 0;
            band_add(&Bp, k, k, aE + aW + aN + aS);
            if (i + 1 < N) band_add(&Bp, k, k + 1, -aE);
            if (i > 0) band_add(&Bp, k, k - 1, -aW);
            if (j + 1 < N) band_add(&Bp, k, k + N, -aN);
            if (j > 0) band_add(&Bp, k, k - N, -aS);
            rp[k] = -h * (u2[(i + 1) * N + j] - u2[i * N + j] + v2[i * (N + 1) + j + 1] - v2[i * (N + 1) + j]);
        }
        band_solve(&Bp, rp);
        for (int k = 0; k < np_; k++) pc[k] = rp[k];
        /* corrections */
        for (int i = 1; i < N; i++) for (int j = 0; j < N; j++)
            U(i, j) = u2[i * N + j] + du[i * N + j] * (pc[j * N + i - 1] - pc[j * N + i]);
        for (int i = 0; i < N; i++) for (int j = 1; j < N; j++)
            V(i, j) = v2[i * (N + 1) + j] + dv[i * (N + 1) + j] * (pc[(j - 1) * N + i] - pc[j * N + i]);
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) P(i, j) += ap * pc[j * N + i];

        double r = residual(u, v, p) / r0;
        if (hist) fprintf(hist, "%s,%d,%.6e\n", tag, it, r);
        if (!isfinite(r) || r > 1e6) { R.iters = it; break; }
        if (r < TOL) { R.iters = it; R.converged = 1; break; }
        R.iters = it;
    }
    R.seconds = now() - t0;
    centreline(u, &R);
    band_free(&Bu); band_free(&Bv); band_free(&Bp);
    free(u); free(v); free(p); free(du); free(dv); free(pc); free(ru); free(rv); free(rp); free(u2); free(v2);
    return R;
}

/* ---------------- COUPLED ---------------- */
static int *iu, *iv, *ip;   /* global unknown index for each u, v, p */
static int pass_kl, pass_ku;
static Band *curB;
static void put(int r, int c, double val) {
    if (!curB) { if (r - c > pass_kl) pass_kl = r - c; if (c - r > pass_ku) pass_ku = c - r; }
    else band_add(curB, r, c, val);
}
static void assemble_coupled(const double *u, const double *v, double *rhs) {
    for (int i = 0; i <= N; i++) for (int j = 0; j < N; j++) {
        if (i == 0 || i == N) continue;
        Coef c = coef_u(u, v, i, j); int r = iu[i * N + j];
        put(r, r, c.aP);
        if (i + 1 < N) put(r, iu[(i + 1) * N + j], -c.aE);
        if (i - 1 > 0) put(r, iu[(i - 1) * N + j], -c.aW);
        if (j < N - 1) put(r, iu[i * N + j + 1], -c.aN);
        if (j > 0) put(r, iu[i * N + j - 1], -c.aS);
        put(r, ip[i * N + j], h); put(r, ip[(i - 1) * N + j], -h);        /* + h (p_E - p_W) */
        if (rhs) rhs[r] = c.b;
    }
    for (int i = 0; i < N; i++) for (int j = 1; j < N; j++) {
        Coef c = coef_v(u, v, i, j); int r = iv[i * (N + 1) + j];
        put(r, r, c.aP);
        if (j + 1 < N) put(r, iv[i * (N + 1) + j + 1], -c.aN);
        if (j - 1 > 0) put(r, iv[i * (N + 1) + j - 1], -c.aS);
        if (i + 1 < N) put(r, iv[(i + 1) * (N + 1) + j], -c.aE);
        if (i > 0) put(r, iv[(i - 1) * (N + 1) + j], -c.aW);
        put(r, ip[i * N + j], h); put(r, ip[i * N + j - 1], -h);
        if (rhs) rhs[r] = c.b;
    }
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        int r = ip[i * N + j];
        if (i == 0 && j == 0) { put(r, r, 1.0); if (rhs) rhs[r] = 0; continue; }   /* pin p(0,0) = 0 */
        if (i + 1 < N) put(r, iu[(i + 1) * N + j], h);
        if (i > 0) put(r, iu[i * N + j], -h);
        if (j + 1 < N) put(r, iv[i * (N + 1) + j + 1], h);
        if (j > 0) put(r, iv[i * (N + 1) + j], -h);
        if (rhs) rhs[r] = 0;
    }
}
static Result coupled(double alpha, FILE *hist, const char *tag) {
    Result R = {0};
    double *u = calloc((N + 1) * N, 8), *v = calloc(N * (N + 1), 8), *p = calloc(N * N, 8);
    iu = malloc((N + 1) * N * sizeof(int)); iv = malloc(N * (N + 1) * sizeof(int)); ip = malloc(N * N * sizeof(int));
    int n = 0;                                   /* interleave unknowns cell by cell: small bandwidth */
    for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
        if (i >= 1) iu[i * N + j] = n++;
        if (j >= 1) iv[i * (N + 1) + j] = n++;
        ip[i * N + j] = n++;
    }
    pass_kl = pass_ku = 0; curB = NULL; assemble_coupled(u, v, NULL);   /* measure bandwidth */
    Band B = band_new(n, pass_kl, pass_ku);
    R.mb = band_mb(&B);
    double *x = malloc(n * 8), r0 = residual(u, v, p), t0 = now();
    for (int it = 1; it <= MAXIT; it++) {
        band_clear(&B); curB = &B; assemble_coupled(u, v, x); curB = NULL;
        band_solve(&B, x);
        for (int i = 1; i < N; i++) for (int j = 0; j < N; j++) U(i, j) = alpha * x[iu[i * N + j]] + (1 - alpha) * U(i, j);
        for (int i = 0; i < N; i++) for (int j = 1; j < N; j++) V(i, j) = alpha * x[iv[i * (N + 1) + j]] + (1 - alpha) * V(i, j);
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) P(i, j) = alpha * x[ip[i * N + j]] + (1 - alpha) * P(i, j);
        double r = residual(u, v, p) / r0;
        if (hist) fprintf(hist, "%s,%d,%.6e\n", tag, it, r);
        R.iters = it;
        if (!isfinite(r) || r > 1e6) break;
        if (r < TOL) { R.converged = 1; break; }
    }
    R.seconds = now() - t0;
    centreline(u, &R);
    band_free(&B); free(x); free(u); free(v); free(p); free(iu); free(iv); free(ip);
    return R;
}

/* =====================================================================
 *  ITERATIVE VERSIONS (no LU): what real codes do
 * ===================================================================== */
static Coef *cu, *cv;        /* momentum coefficients per face */
static double *du_, *dv_;    /* h / aP (relaxed or not): pressure-velocity link */
static long work_gs, work_cg, work_mv;   /* counters: GS sweeps, CG iterations, matvecs */

static void build_coefs(const double *u, const double *v, double alpha) {
    for (int i = 1; i < N; i++) for (int j = 0; j < N; j++) {
        cu[i * N + j] = coef_u(u, v, i, j);
        du_[i * N + j] = h / (cu[i * N + j].aP / alpha);
    }
    for (int i = 0; i < N; i++) for (int j = 1; j < N; j++) {
        cv[i * (N + 1) + j] = coef_v(u, v, i, j);
        dv_[i * (N + 1) + j] = h / (cv[i * (N + 1) + j].aP / alpha);
    }
}
/* one Gauss-Seidel sweep on (aP/alpha) x - sum a_nb x_nb = r, for u and v */
static void gs_u(double *x, const double *r, double alpha) {
    for (int i = 1; i < N; i++) for (int j = 0; j < N; j++) {
        Coef c = cu[i * N + j]; int k = i * N + j;
        double s = r[k] + c.aE * x[k + N] + c.aW * x[k - N];
        if (j < N - 1) s += c.aN * x[k + 1];
        if (j > 0) s += c.aS * x[k - 1];
        x[k] = s / (c.aP / alpha);
    }
}
static void gs_v(double *x, const double *r, double alpha) {
    for (int i = 0; i < N; i++) for (int j = 1; j < N; j++) {
        Coef c = cv[i * (N + 1) + j]; int k = i * (N + 1) + j;
        double s = r[k] + c.aN * x[k + 1] + c.aS * x[k - 1];
        if (i < N - 1) s += c.aE * x[k + (N + 1)];
        if (i > 0) s += c.aW * x[k - (N + 1)];
        x[k] = s / (c.aP / alpha);
    }
}
static double res_u(const double *x, const double *r, double alpha) {
    double s2 = 0;
    for (int i = 1; i < N; i++) for (int j = 0; j < N; j++) {
        Coef c = cu[i * N + j]; int k = i * N + j;
        double a = r[k] + c.aE * x[k + N] + c.aW * x[k - N] + (j < N - 1 ? c.aN * x[k + 1] : 0) + (j > 0 ? c.aS * x[k - 1] : 0)
                 - c.aP / alpha * x[k];
        s2 += a * a;
    }
    return sqrt(s2);
}
static double res_v(const double *x, const double *r, double alpha) {
    double s2 = 0;
    for (int i = 0; i < N; i++) for (int j = 1; j < N; j++) {
        Coef c = cv[i * (N + 1) + j]; int k = i * (N + 1) + j;
        double a = r[k] + c.aN * x[k + 1] + c.aS * x[k - 1] + (i < N - 1 ? c.aE * x[k + N + 1] : 0) + (i > 0 ? c.aW * x[k - N - 1] : 0)
                 - c.aP / alpha * x[k];
        s2 += a * a;
    }
    return sqrt(s2);
}
/* GS until the residual drops by relTol (or maxSweeps) */
static void mom_solve(double *xu, const double *ru, double *xv, const double *rv, double alpha, double relTol, int maxSweeps) {
    double r0u = res_u(xu, ru, alpha), r0v = res_v(xv, rv, alpha);
    for (int s = 0; s < maxSweeps; s++) {
        gs_u(xu, ru, alpha); gs_v(xv, rv, alpha); work_gs++;
        if (res_u(xu, ru, alpha) <= relTol * r0u && res_v(xv, rv, alpha) <= relTol * r0v) break;
    }
}
/* pressure (correction) operator M x = sum a_nb (x_P - x_nb), a = h * d; cell 0 is fixed (Dirichlet) */
static void pmat(const double *x, double *y) {
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        int k = i * N + j;
        if (k == 0) { y[0] = 0; continue; }
        double s = 0;
        if (i + 1 < N) { double a = h * du_[(i + 1) * N + j]; s += a * (x[k] - ((i + 1) * N + j ? x[(i + 1) * N + j] : 0)); }
        if (i > 0)     { double a = h * du_[i * N + j];       s += a * (x[k] - ((i - 1) * N + j ? x[(i - 1) * N + j] : 0)); }
        if (j + 1 < N) { double a = h * dv_[i * (N + 1) + j + 1]; s += a * (x[k] - x[k + 1]); }
        if (j > 0)     { double a = h * dv_[i * (N + 1) + j]; s += a * (x[k] - (k - 1 ? x[k - 1] : 0)); }
        y[k] = s;
    }
}
static double pdiag(int i, int j) {
    double d = 0;
    if (i + 1 < N) d += h * du_[(i + 1) * N + j];
    if (i > 0) d += h * du_[i * N + j];
    if (j + 1 < N) d += h * dv_[i * (N + 1) + j + 1];
    if (j > 0) d += h * dv_[i * (N + 1) + j];
    return d;
}
/* Jacobi-preconditioned CG on M x = b with x[0] = x0 fixed (moved to the rhs) */
static void pcg(double *x, const double *b_in, double x0, double relTol, int maxit) {
    int n = N * N;
    double *b = malloc(n * 8), *r = malloc(n * 8), *z = malloc(n * 8), *pp = malloc(n * 8), *q = malloc(n * 8);
    memcpy(b, b_in, n * 8); b[0] = 0;
    /* neighbours of cell 0 (i=0,j=0): cell (1,0) via u face 1, cell (0,1) via v face 1 */
    if (N > 1) { b[1 * N + 0] += h * du_[1 * N + 0] * x0; b[0 * N + 1] += h * dv_[0 * (N + 1) + 1] * x0; }
    x[0] = 0;
    pmat(x, q);
    double rz = 0;
    for (int k = 1; k < n; k++) r[k] = b[k] - q[k];
    r[0] = 0;
    for (int k = 0; k < n; k++) { z[k] = k ? r[k] / pdiag(k / N, k % N) : 0; pp[k] = z[k]; rz += r[k] * z[k]; }
    double r0 = 0; for (int k = 1; k < n; k++) r0 += r[k] * r[k]; r0 = sqrt(r0);
    for (int it = 0; it < maxit && r0 > 1e-300; it++) {
        pmat(pp, q); work_cg++;
        double pq = 0; for (int k = 1; k < n; k++) pq += pp[k] * q[k];
        double al = rz / pq, rn = 0;
        for (int k = 1; k < n; k++) { x[k] += al * pp[k]; r[k] -= al * q[k]; rn += r[k] * r[k]; }
        if (sqrt(rn) <= relTol * r0) break;
        double rz2 = 0;
        for (int k = 1; k < n; k++) { z[k] = r[k] / pdiag(k / N, k % N); rz2 += r[k] * z[k]; }
        double be = rz2 / rz; rz = rz2;
        for (int k = 1; k < n; k++) pp[k] = z[k] + be * pp[k];
    }
    x[0] = x0;
    free(b); free(r); free(z); free(pp); free(q);
}

/* ---- SIMPLE with iterative inner solvers ---- */
static Result simple_iter(double au, double ap, double relU, double relP, FILE *hist, const char *tag) {
    Result R = {0};
    int nU = (N + 1) * N, nV = N * (N + 1), nP = N * N;
    double *u = calloc(nU, 8), *v = calloc(nV, 8), *p = calloc(nP, 8), *pc = calloc(nP, 8);
    double *ru = calloc(nU, 8), *rv = calloc(nV, 8), *rp = calloc(nP, 8);
    cu = calloc(nU, sizeof(Coef)); cv = calloc(nV, sizeof(Coef)); du_ = calloc(nU, 8); dv_ = calloc(nV, 8);
    R.mb = (3.0 * nU + 3.0 * nV + 8.0 * nP) * 8 / 1e6 + (nU + nV) * sizeof(Coef) / 1e6;
    work_gs = work_cg = work_mv = 0;
    double r0 = residual(u, v, p), t0 = now();
    for (int it = 1; it <= MAXIT; it++) {
        build_coefs(u, v, au);
        for (int i = 1; i < N; i++) for (int j = 0; j < N; j++)
            ru[i * N + j] = cu[i * N + j].b + h * (P(i - 1, j) - P(i, j)) + (1 - au) * cu[i * N + j].aP / au * U(i, j);
        for (int i = 0; i < N; i++) for (int j = 1; j < N; j++)
            rv[i * (N + 1) + j] = cv[i * (N + 1) + j].b + h * (P(i, j - 1) - P(i, j)) + (1 - au) * cv[i * (N + 1) + j].aP / au * V(i, j);
        mom_solve(u, ru, v, rv, au, relU, 1000);                    /* u*, v* (in place) */
        for (int i = 0; i < N; i++) for (int j = 0; j < N; j++)
            rp[i * N + j] = -h * (U(i + 1, j) - U(i, j) + V(i, j + 1) - V(i, j));
        memset(pc, 0, nP * 8);
        pcg(pc, rp, 0.0, relP, 20000);
        for (int i = 1; i < N; i++) for (int j = 0; j < N; j++) U(i, j) += du_[i * N + j] * (pc[(i - 1) * N + j] - pc[i * N + j]);
        for (int i = 0; i < N; i++) for (int j = 1; j < N; j++) V(i, j) += dv_[i * (N + 1) + j] * (pc[i * N + j - 1] - pc[i * N + j]);
        for (int k = 0; k < nP; k++) p[k] += ap * pc[k];
        double r = residual(u, v, p) / r0;
        if (hist) fprintf(hist, "%s,%d,%.6e\n", tag, it, r);
        R.iters = it;
        if (!isfinite(r) || r > 1e6) break;
        if (r < TOL) { R.converged = 1; break; }
    }
    R.seconds = now() - t0;
    centreline(u, &R);
    free(u); free(v); free(p); free(pc); free(ru); free(rv); free(rp); free(cu); free(cv); free(du_); free(dv_);
    return R;
}

/* ---- coupled with FGMRES + SIMPLE-type block preconditioner ---- */
static int OU, OV, OP_, LEN;
static void kmat(const double *x, double *y) {          /* y = K x, full saddle-point operator */
    const double *xu = x + OU, *xv = x + OV, *xp = x + OP_;
    double *yu = y + OU, *yv = y + OV, *yp = y + OP_;
    work_mv++;
    for (int k = 0; k < (N + 1) * N; k++) yu[k] = xu[k];
    for (int k = 0; k < N * (N + 1); k++) yv[k] = xv[k];
    for (int i = 1; i < N; i++) for (int j = 0; j < N; j++) {
        Coef c = cu[i * N + j]; int k = i * N + j;
        yu[k] = c.aP * xu[k] - c.aE * xu[k + N] - c.aW * xu[k - N] - (j < N - 1 ? c.aN * xu[k + 1] : 0)
              - (j > 0 ? c.aS * xu[k - 1] : 0) + h * (xp[i * N + j] - xp[(i - 1) * N + j]);
    }
    for (int i = 0; i < N; i++) for (int j = 1; j < N; j++) {
        Coef c = cv[i * (N + 1) + j]; int k = i * (N + 1) + j;
        yv[k] = c.aP * xv[k] - c.aN * xv[k + 1] - c.aS * xv[k - 1] - (i < N - 1 ? c.aE * xv[k + N + 1] : 0)
              - (i > 0 ? c.aW * xv[k - N - 1] : 0) + h * (xp[i * N + j] - xp[i * N + j - 1]);
    }
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        int k = i * N + j;
        yp[k] = k ? h * (xu[(i + 1) * N + j] - xu[i * N + j] + xv[i * (N + 1) + j + 1] - xv[i * (N + 1) + j]) : xp[0];
    }
}
static double *pre_tmp;
static void kprec(const double *r, double *z) {          /* z = P^-1 r (SIMPLE-type block preconditioner) */
    double *zu = z + OU, *zv = z + OV, *zp = z + OP_;
    memset(z, 0, LEN * 8);
    for (int s = 0; s < 2; s++) { gs_u(zu, r + OU, 1.0); gs_v(zv, r + OV, 1.0); work_gs++; }   /* z* ~ A^-1 r_u */
    double *rhs = pre_tmp;
    for (int i = 0; i < N; i++) for (int j = 0; j < N; j++) {
        int k = i * N + j;
        rhs[k] = r[OP_ + k] - h * (zu[(i + 1) * N + j] - zu[i * N + j] + zv[i * (N + 1) + j + 1] - zv[i * (N + 1) + j]);
    }
    pcg(zp, rhs, r[OP_], 0.1, 500);                                   /* M z_p = r_p - D z* */
    for (int i = 1; i < N; i++) for (int j = 0; j < N; j++) zu[i * N + j] -= du_[i * N + j] * (zp[i * N + j] - zp[(i - 1) * N + j]);
    for (int i = 0; i < N; i++) for (int j = 1; j < N; j++) zv[i * (N + 1) + j] -= dv_[i * (N + 1) + j] * (zp[i * N + j] - zp[i * N + j - 1]);
}
#define RESTART 30
static void fgmres(double *x, const double *b, double relTol, int maxit) {
    int m = RESTART;
    double *V_ = malloc((size_t)(m + 1) * LEN * 8), *Z = malloc((size_t)m * LEN * 8), *w = malloc(LEN * 8);
    double H[RESTART + 1][RESTART], cs[RESTART], sn[RESTART], g[RESTART + 1];
    double bn = 0; for (int k = 0; k < LEN; k++) bn += b[k] * b[k]; bn = sqrt(bn);
    double rinit = -1;
    for (int outer = 0; outer < maxit; outer++) {
        kmat(x, w);
        double beta = 0; for (int k = 0; k < LEN; k++) { V_[k] = b[k] - w[k]; beta += V_[k] * V_[k]; }
        beta = sqrt(beta);
        if (rinit < 0) rinit = beta;
        if (beta <= relTol * rinit || beta < 1e-14 * bn) break;
        for (int k = 0; k < LEN; k++) V_[k] /= beta;
        memset(g, 0, sizeof g); g[0] = beta;
        int jj;
        for (jj = 0; jj < m; jj++) {
            double *vj = V_ + (size_t)jj * LEN, *zj = Z + (size_t)jj * LEN, *vn = V_ + (size_t)(jj + 1) * LEN;
            kprec(vj, zj); kmat(zj, w);
            for (int i = 0; i <= jj; i++) {                           /* modified Gram-Schmidt */
                double *vi = V_ + (size_t)i * LEN, d = 0;
                for (int k = 0; k < LEN; k++) d += w[k] * vi[k];
                H[i][jj] = d;
                for (int k = 0; k < LEN; k++) w[k] -= d * vi[k];
            }
            double nw = 0; for (int k = 0; k < LEN; k++) nw += w[k] * w[k]; nw = sqrt(nw);
            H[jj + 1][jj] = nw;
            for (int k = 0; k < LEN; k++) vn[k] = nw > 0 ? w[k] / nw : 0;
            for (int i = 0; i < jj; i++) {                            /* apply old rotations */
                double t = cs[i] * H[i][jj] + sn[i] * H[i + 1][jj];
                H[i + 1][jj] = -sn[i] * H[i][jj] + cs[i] * H[i + 1][jj]; H[i][jj] = t;
            }
            double d = hypot(H[jj][jj], H[jj + 1][jj]);
            cs[jj] = H[jj][jj] / d; sn[jj] = H[jj + 1][jj] / d;
            H[jj][jj] = d; H[jj + 1][jj] = 0;
            g[jj + 1] = -sn[jj] * g[jj]; g[jj] = cs[jj] * g[jj];
            if (fabs(g[jj + 1]) <= relTol * rinit) { jj++; break; }
        }
        double y[RESTART];                                            /* back substitution */
        for (int i = jj - 1; i >= 0; i--) {
            double s = g[i];
            for (int k = i + 1; k < jj; k++) s -= H[i][k] * y[k];
            y[i] = s / H[i][i];
        }
        for (int i = 0; i < jj; i++) { double *zi = Z + (size_t)i * LEN; for (int k = 0; k < LEN; k++) x[k] += y[i] * zi[k]; }
        if (fabs(g[jj]) <= relTol * rinit) break;
    }
    free(V_); free(Z); free(w);
}
static Result coupled_iter(double relLin, FILE *hist, const char *tag) {
    Result R = {0};
    int nU = (N + 1) * N, nV = N * (N + 1), nP = N * N;
    OU = 0; OV = nU; OP_ = nU + nV; LEN = nU + nV + nP;
    double *x = calloc(LEN, 8), *b = calloc(LEN, 8);
    cu = calloc(nU, sizeof(Coef)); cv = calloc(nV, sizeof(Coef)); du_ = calloc(nU, 8); dv_ = calloc(nV, 8);
    pre_tmp = malloc(nP * 8);
    R.mb = (2.0 * RESTART + 4) * LEN * 8 / 1e6 + (nU + nV) * sizeof(Coef) / 1e6;
    work_gs = work_cg = work_mv = 0;
    double *u = x + OU, *v = x + OV, *p = x + OP_;
    double r0 = residual(u, v, p), t0 = now();
    for (int it = 1; it <= MAXIT; it++) {
        build_coefs(u, v, 1.0);
        memset(b, 0, LEN * 8);
        for (int i = 1; i < N; i++) for (int j = 0; j < N; j++) b[OU + i * N + j] = cu[i * N + j].b;
        fgmres(x, b, relLin, 200);
        double r = residual(u, v, p) / r0;
        if (hist) fprintf(hist, "%s,%d,%.6e\n", tag, it, r);
        R.iters = it;
        if (!isfinite(r) || r > 1e6) break;
        if (r < TOL) { R.converged = 1; break; }
    }
    R.seconds = now() - t0;
    centreline(u, &R);
    free(x); free(b); free(cu); free(cv); free(du_); free(dv_); free(pre_tmp);
    return R;
}

static void print(const char *name, Result R) {
    if (R.converged)
        printf("   %-26s %6d %9.3f %9.2f %10.5f\n", name, R.iters, R.seconds, R.mb, R.umin);
    else
        printf("   %-26s %6s %9.3f %9.2f %10s\n", name, R.iters >= MAXIT ? "max" : "DIVERGED", R.seconds, R.mb, "-");
}

int main(int argc, char **argv) {
    const char *dir = argc > 1 ? argv[1] : ".";
    N = argc > 2 ? atoi(argv[2]) : 64; h = 1.0 / N; nu = 1.0 / 100;
    char path[512];
    snprintf(path, sizeof path, "%s/simple_reltol_N%d.csv", dir, N); FILE *fs = fopen(path, "w");
    fprintf(fs, "relU,relP,outer,seconds,gs_sweeps,cg_iters,converged\n");
    printf("SIMPLE (0.9 / 0.2), N = %d, Re = 100: outer iterations | time (s) | GS sweeps | CG iterations\n", N);
    double rU[] = {0.5, 0.1, 0.01}, rP[] = {0.9, 0.5, 0.2, 0.1, 0.05, 0.01, 1e-3, 1e-6};
    for (int a = 0; a < 3; a++) for (int b = 0; b < 8; b++) {
        Result R = simple_iter(0.9, 0.2, rU[a], rP[b], NULL, "");
        printf("  relU %-5g relP %-6g  %6d %8.3f %8ld %9ld %s\n", rU[a], rP[b], R.iters, R.seconds, work_gs, work_cg, R.converged ? "" : "NOT CONVERGED");
        fprintf(fs, "%g,%g,%d,%.4f,%ld,%ld,%d\n", rU[a], rP[b], R.iters, R.seconds, work_gs, work_cg, R.converged);
        fflush(stdout);
    }
    fclose(fs);
    snprintf(path, sizeof path, "%s/coupled_reltol_N%d.csv", dir, N); FILE *fc = fopen(path, "w");
    fprintf(fc, "relLin,outer,seconds,matvecs,cg_iters,converged\n");
    printf("coupled FGMRES, N = %d: outer | time (s) | matvecs | CG iterations\n", N);
    double rL[] = {0.9, 0.5, 0.2, 0.1, 0.05, 0.01, 1e-3, 1e-6};
    for (int b = 0; b < 8; b++) {
        Result R = coupled_iter(rL[b], NULL, "");
        printf("  relTol %-6g  %6d %8.3f %8ld %9ld %s\n", rL[b], R.iters, R.seconds, work_mv, work_cg, R.converged ? "" : "NOT CONVERGED");
        fprintf(fc, "%g,%d,%.4f,%ld,%ld,%d\n", rL[b], R.iters, R.seconds, work_mv, work_cg, R.converged);
        fflush(stdout);
    }
    fclose(fc);
    return 0;
}
