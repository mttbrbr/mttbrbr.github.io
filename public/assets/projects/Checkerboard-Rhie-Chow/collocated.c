/**
 * Pressure checkerboarding on a collocated grid, and how Rhie-Chow removes it.
 *
 * Steady 2D lid-driven cavity, Re = 100, N x N cells. Velocity AND pressure
 * stored at the cell centres (collocated, like OpenFOAM). SIMPLE algorithm,
 * first-order upwind convection, Gauss-Seidel inner solvers.
 *
 * The only difference between the two runs is how the face velocity (and so
 * the mass flux) is obtained from the cell values:
 *   linear     : u_f = average of the two cells
 *   Rhie-Chow  : u_f = average - D_f [ (p_N - p_P)/h - average of cell gradients ]
 *
 * Build: gcc -O3 -march=native collocated.c -o collocated -lm
 * Run:   ./collocated outdir [N]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static int N; static double h, nu = 0.01, au = 0.7, ap = 0.3;
#define C(i, j) ((j) * N + (i))

/* cell fields */
static double *u, *v, *p, *pc, *aPu, *Du, *gx, *gy, *bu, *bv;
static double *aE, *aW, *aN, *aS;          /* momentum neighbour coefficients (same for u and v) */
/* face velocities: ue[i + (N+1) j] on vertical faces x = i h, vn[i + N j] on horizontal faces y = j h */
static double *ue, *vn;
#define UE(i, j) ue[(i) + (N + 1) * (j)]
#define VN(i, j) vn[(i) + N * (j)]

static int rhiechow;          /* 0 linear, 1 Rhie-Chow, 2 Rhie-Chow + balanced body force */
static double lid = 1.0, g = 0.0;
static double fface(int j) { double y = j * h; return y < 0.5 - 1e-12 ? g : (y > 0.5 + 1e-12 ? -g : 0.0); }   /* force at face y = j h */
static double fcell(int j) {
    if (rhiechow == 2) return 0.5 * ((j + 1 < N ? fface(j + 1) : 0) + (j > 0 ? fface(j) : 0));   /* same reconstruction as grad p */
    double y = (j + 0.5) * h; return y < 0.5 ? g : -g;
}

static void pressure_gradient(const double *q, double *gxo, double *gyo) {
    /* Gauss gradient with linear face interpolation, zero-gradient at walls */
    for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
        double qe = i < N - 1 ? 0.5 * (q[C(i, j)] + q[C(i + 1, j)]) : q[C(i, j)];
        double qw = i > 0 ? 0.5 * (q[C(i, j)] + q[C(i - 1, j)]) : q[C(i, j)];
        double qn = j < N - 1 ? 0.5 * (q[C(i, j)] + q[C(i, j + 1)]) : q[C(i, j)];
        double qs = j > 0 ? 0.5 * (q[C(i, j)] + q[C(i, j - 1)]) : q[C(i, j)];
        gxo[C(i, j)] = (qe - qw) / h; gyo[C(i, j)] = (qn - qs) / h;
    }
}

static void momentum_coefficients(void) {
    for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
        int c = C(i, j);
        double Fe = h * UE(i + 1, j), Fw = h * UE(i, j), Fn = h * VN(i, j + 1), Fs = h * VN(i, j);
        double D = nu, sum = 0;
        aE[c] = i < N - 1 ? D + fmax(-Fe, 0) : 0;  aW[c] = i > 0 ? D + fmax(Fw, 0) : 0;
        aN[c] = j < N - 1 ? D + fmax(-Fn, 0) : 0;  aS[c] = j > 0 ? D + fmax(Fs, 0) : 0;
        sum = aE[c] + aW[c] + aN[c] + aS[c];
        double wall = 0; bu[c] = 0; bv[c] = 0;
        if (i == 0) wall += 2 * nu;
        if (i == N - 1) wall += 2 * nu;
        if (j == 0) wall += 2 * nu;
        if (j == N - 1) { wall += 2 * nu; bu[c] += 2 * nu * lid; }     /* lid moves with u = lid */
        aPu[c] = (sum + wall + (Fe - Fw + Fn - Fs)) / au;              /* relaxed */
        Du[c] = h * h / aPu[c];
    }
}

/* Gauss-Seidel on  aPu x = sum a_nb x_nb + rhs ; until residual / initial < rel */
static void solve_momentum(double *x, const double *rhs, double rel) {
    double r0 = -1;
    for (int s = 0; s < 100; s++) {
        double r2 = 0;
        for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
            int c = C(i, j);
            double nb = (i < N - 1 ? aE[c] * x[c + 1] : 0) + (i > 0 ? aW[c] * x[c - 1] : 0)
                      + (j < N - 1 ? aN[c] * x[c + N] : 0) + (j > 0 ? aS[c] * x[c - N] : 0);
            double xn = (nb + rhs[c]) / aPu[c];
            r2 += (xn - x[c]) * (xn - x[c]) * aPu[c] * aPu[c];
            x[c] = xn;
        }
        if (r0 < 0) r0 = sqrt(r2);
        if (sqrt(r2) <= rel * r0) break;
    }
}

/* face velocity from cell values */
static void face_velocities(void) {
    pressure_gradient(p, gx, gy);
    for (int j = 0; j < N; j++) for (int i = 1; i < N; i++) {
        int P = C(i - 1, j), E = C(i, j);
        double f = 0.5 * (u[P] + u[E]);
        if (rhiechow) {
            double Df = 0.5 * (Du[P] + Du[E]);
            f -= Df * ((p[E] - p[P]) / h - 0.5 * (gx[P] + gx[E]));
        }
        UE(i, j) = f;
    }
    for (int j = 1; j < N; j++) for (int i = 0; i < N; i++) {
        int P = C(i, j - 1), Nn = C(i, j);
        double f = 0.5 * (v[P] + v[Nn]);
        if (rhiechow) {
            double Df = 0.5 * (Du[P] + Du[Nn]);
            f -= Df * ((p[Nn] - p[P]) / h - 0.5 * (gy[P] + gy[Nn]));
        }

        if (rhiechow == 2) {   /* balanced force: treat the body force like the pressure gradient */
            double Df = 0.5 * (Du[P] + Du[Nn]);
            f += Df * (fface(j) - 0.5 * (fcell(j - 1) + fcell(j)));
        }
        VN(i, j) = f;
    }
    /* walls: zero normal velocity */
    for (int j = 0; j < N; j++) { UE(0, j) = 0; UE(N, j) = 0; }
    for (int i = 0; i < N; i++) { VN(i, 0) = 0; VN(i, N) = 0; }
}

/* pressure-correction operator: L p' = continuity of the face-velocity corrections.
   Rhie-Chow : u'_f = -D_f (p'_N - p'_P) / h          (compact 5-point stencil)
   linear    : u'_f = average of u'_P = -D_P (grad p')_P  (wide stencil, every other cell) */
static double *tgx, *tgy;
static void apply_L(const double *q, double *out) {
    if (rhiechow) {
        for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
            int c = C(i, j); double s = 0;
            if (i < N - 1) s += 0.5 * (Du[c] + Du[c + 1]) * (q[c] - q[c + 1]);
            if (i > 0)     s += 0.5 * (Du[c] + Du[c - 1]) * (q[c] - q[c - 1]);
            if (j < N - 1) s += 0.5 * (Du[c] + Du[c + N]) * (q[c] - q[c + N]);
            if (j > 0)     s += 0.5 * (Du[c] + Du[c - N]) * (q[c] - q[c - N]);
            out[c] = s;                                   /* = -h * div(u') */
        }
    } else {
        pressure_gradient(q, tgx, tgy);
        for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
            int c = C(i, j); double s = 0;
            if (i < N - 1) s += 0.5 * (Du[c] * tgx[c] + Du[c + 1] * tgx[c + 1]) * h;
            if (i > 0)     s -= 0.5 * (Du[c] * tgx[c] + Du[c - 1] * tgx[c - 1]) * h;
            if (j < N - 1) s += 0.5 * (Du[c] * tgy[c] + Du[c + N] * tgy[c + N]) * h;
            if (j > 0)     s -= 0.5 * (Du[c] * tgy[c] + Du[c - N] * tgy[c - N]) * h;
            out[c] = -s;                                  /* = -h * div(u') with u' = -D grad p' */
        }
    }
}
static double mass_imbalance(double *b) {   /* b = -h * div(u_f), returns L2 norm */
    double s2 = 0;
    for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
        double d = h * (UE(i + 1, j) - UE(i, j) + VN(i, j + 1) - VN(i, j));
        b[C(i, j)] = -d; s2 += d * d;
    }
    return sqrt(s2);
}
/* conjugate gradient on L p' = b (L is symmetric positive semi-definite in both cases;
   b sums to zero, so the system is consistent). Mean of p' removed. */
static void solve_pc(double *x, const double *b, double rel, int maxit) {
    int n = N * N;
    double *r = malloc(n * 8), *pp = malloc(n * 8), *q = malloc(n * 8);
    memset(x, 0, n * 8); memcpy(r, b, n * 8); memcpy(pp, b, n * 8);
    double rr = 0; for (int k = 0; k < n; k++) rr += r[k] * r[k];
    double r0 = sqrt(rr);
    for (int it = 0; it < maxit && r0 > 0; it++) {
        apply_L(pp, q);
        double pq = 0; for (int k = 0; k < n; k++) pq += pp[k] * q[k];
        if (pq <= 0) break;
        double al = rr / pq, rn = 0;
        for (int k = 0; k < n; k++) { x[k] += al * pp[k]; r[k] -= al * q[k]; rn += r[k] * r[k]; }
        if (sqrt(rn) <= rel * r0) break;
        double be = rn / rr; rr = rn;
        for (int k = 0; k < n; k++) pp[k] = r[k] + be * pp[k];
    }
    double m = 0; for (int k = 0; k < n; k++) m += x[k]; m /= n;
    for (int k = 0; k < n; k++) x[k] -= m;
    free(r); free(pp); free(q);
}

/* checkerboard content: with central differences on a collocated grid, the
   pressure modes (-1)^i, (-1)^j and (-1)^(i+j) have zero gradient in the
   interior. Amplitude = RMS of the projections of p on these three modes. */
static double checker_amplitude(void) {
    double m[3] = {0, 0, 0};
    for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
        double q = p[C(i, j)];
        m[0] += q * ((i & 1) ? -1 : 1); m[1] += q * ((j & 1) ? -1 : 1); m[2] += q * (((i + j) & 1) ? -1 : 1);
    }
    return sqrt(m[0] * m[0] + m[1] * m[1] + m[2] * m[2]) / (N * N);
}

int main(int argc, char **argv) {
    const char *dir = argc > 1 ? argv[1] : ".";
    N = argc > 2 ? atoi(argv[2]) : 32; h = 1.0 / N;
    double amp0 = argc > 3 ? atof(argv[3]) : 0.0;   /* initial checkerboard amplitude in p */
    g = argc > 4 ? atof(argv[4]) : 0.0;             /* body force in y: +g below y = 0.5, -g above */
    lid = argc > 5 ? atof(argv[5]) : 1.0;
    int n = N * N;
    char path[512];
    snprintf(path, sizeof path, "%s/history_N%d_a%g_g%g.csv", dir, N, amp0, g); FILE *hf = fopen(path, "w");
    fprintf(hf, "method,iter,mass_imbalance\n");
    for (rhiechow = 0; rhiechow <= 2; rhiechow++) {
        u = calloc(n, 8); v = calloc(n, 8); p = calloc(n, 8); pc = calloc(n, 8);
        aPu = calloc(n, 8); Du = calloc(n, 8); gx = calloc(n, 8); gy = calloc(n, 8);
        bu = calloc(n, 8); bv = calloc(n, 8); tgx = calloc(n, 8); tgy = calloc(n, 8);
        aE = calloc(n, 8); aW = calloc(n, 8); aN = calloc(n, 8); aS = calloc(n, 8);
        ue = calloc((N + 1) * N, 8); vn = calloc(N * (N + 1), 8);
        double *rhs = malloc(n * 8), *b = malloc(n * 8);
        const char *name = rhiechow == 0 ? "linear" : (rhiechow == 1 ? "rhie-chow" : "balanced");
        for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) p[C(i, j)] = amp0 * (((i + j) & 1) ? -1 : 1);
        double ca0 = checker_amplitude();
        double m0 = -1; int it;
        for (it = 1; it <= 20000; it++) {
            momentum_coefficients();
            pressure_gradient(p, gx, gy);
            for (int c = 0; c < n; c++) rhs[c] = bu[c] - h * h * gx[c] + (1 - au) * aPu[c] * u[c];
            solve_momentum(u, rhs, 0.1);
            for (int c = 0; c < n; c++)
                rhs[c] = bv[c] - h * h * gy[c] + h * h * fcell(c / N) + (1 - au) * aPu[c] * v[c];
            solve_momentum(v, rhs, 0.1);
            face_velocities();
            double m = mass_imbalance(b);
            if (m0 < 0) m0 = m;
            fprintf(hf, "%s,%d,%.6e\n", name, it, m / m0);
            if (m / m0 < 1e-7 && it > 10) break;
            solve_pc(pc, b, 0.01, 4 * n);
            /* correct face velocities (consistent with the operator), cell velocities and pressure */
            if (rhiechow) {   /* Rhie-Chow and balanced: same correction */
                for (int j = 0; j < N; j++) for (int i = 1; i < N; i++)
                    UE(i, j) -= 0.5 * (Du[C(i - 1, j)] + Du[C(i, j)]) * (pc[C(i, j)] - pc[C(i - 1, j)]) / h;
                for (int j = 1; j < N; j++) for (int i = 0; i < N; i++)
                    VN(i, j) -= 0.5 * (Du[C(i, j - 1)] + Du[C(i, j)]) * (pc[C(i, j)] - pc[C(i, j - 1)]) / h;
            }
            pressure_gradient(pc, tgx, tgy);
            for (int c = 0; c < n; c++) { u[c] -= Du[c] * tgx[c]; v[c] -= Du[c] * tgy[c]; p[c] += ap * pc[c]; }
            if (!rhiechow) {
                for (int j = 0; j < N; j++) for (int i = 1; i < N; i++) UE(i, j) = 0.5 * (u[C(i - 1, j)] + u[C(i, j)]);
                for (int j = 1; j < N; j++) for (int i = 0; i < N; i++) VN(i, j) = 0.5 * (v[C(i, j - 1)] + v[C(i, j)]);
            }
        }
        /* outputs: centreline u (x = 0.5, average of the two central columns), pressure field */
        double umin = 1e9;
        for (int j = 0; j < N; j++) { double uc = 0.5 * (u[C(N / 2 - 1, j)] + u[C(N / 2, j)]); if (uc < umin) umin = uc; }
        double mean = 0, pmin = 1e9, pmax = -1e9;
        for (int c = 0; c < n; c++) mean += p[c];
        mean /= n;
        for (int c = 0; c < n; c++) { if (p[c] - mean < pmin) pmin = p[c] - mean; if (p[c] - mean > pmax) pmax = p[c] - mean; }
        double vmax = 0, perr = 0, pex_mean = 0;
        for (int c = 0; c < n; c++) { double a = hypot(u[c], v[c]); if (a > vmax) vmax = a; }
        for (int j = 0; j < N; j++) { double y = (j + 0.5) * h; pex_mean += N * (y < 0.5 ? g * y : g * (1 - y)); }
        pex_mean /= n;
        for (int c = 0; c < n; c++) { double y = (c / N + 0.5) * h, pe = (y < 0.5 ? g * y : g * (1 - y)) - pex_mean;
                                      double e = fabs(p[c] - mean - pe); if (e > perr) perr = e; }
        printf("%-10s N = %3d  start %.2f  outer %5d  u_min = %.5f  p [%.3f, %.3f]  checker %.4f -> %.4f  |u|max %.2e  p err %.2e\n",
               name, N, amp0, it, umin, pmin, pmax, ca0, checker_amplitude(), vmax, perr);
        snprintf(path, sizeof path, "%s/p_%s_N%d_a%g_g%g.csv", dir, name, N, amp0, g); FILE *f = fopen(path, "w");
        for (int j = 0; j < N; j++) { for (int i = 0; i < N; i++) fprintf(f, "%s%.6e", i ? "," : "", p[C(i, j)] - mean); fprintf(f, "\n"); }
        fclose(f);
        snprintf(path, sizeof path, "%s/vel_%s_N%d_a%g_g%g.csv", dir, name, N, amp0, g); f = fopen(path, "w");
        for (int j = 0; j < N; j++) { for (int i = 0; i < N; i++) fprintf(f, "%s%.6e", i ? "," : "", hypot(u[C(i, j)], v[C(i, j)])); fprintf(f, "\n"); }
        fclose(f);
        snprintf(path, sizeof path, "%s/u_%s_N%d_a%g_g%g.csv", dir, name, N, amp0, g); f = fopen(path, "w");
        for (int j = 0; j < N; j++) fprintf(f, "%.5f,%.6e\n", (j + 0.5) * h, 0.5 * (u[C(N / 2 - 1, j)] + u[C(N / 2, j)]));
        fclose(f);
        free(u); free(v); free(p); free(pc); free(aPu); free(Du); free(gx); free(gy); free(bu); free(bv);
        free(tgx); free(tgy); free(aE); free(aW); free(aN); free(aS); free(ue); free(vn); free(rhs); free(b);
    }
    fclose(hf);
    return 0;
}
