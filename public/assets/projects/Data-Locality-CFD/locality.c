/**
 * Data locality in a tiny CFD solver (2D lid-driven cavity, fractional step,
 * Jacobi pressure solver). Derived from ../solver_previo/solver.c.
 *
 * The *numerics* are identical in every build. Only two things change:
 *   - LAYOUT_AOS : fields stored as Array of Structures (one struct per cell)
 *                  instead of Structure of Arrays (one array per field)
 *   - ORDER_JI   : loops nested j-outer / i-inner, i.e. walking the grid
 *                  "against" the row-major storage (stride = N doubles)
 *
 * Build: see Makefile (make all). Run: ./bin_soa_ij N [steps]
 */
#define _POSIX_C_SOURCE 200112L
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define PRESSURE_ITERS 40   /* even: Jacobi ping-pong ends in p */

static size_t N;            /* grid is N x N, row-major: k = i*N + j */
#define IDX(i, j) ((size_t)(i) * N + (size_t)(j))

/* ---------------- memory layout ---------------- */
#ifdef LAYOUT_AOS
typedef struct { double u, v, us, vs, p, pn; } Cell;   /* 48 B per cell */
static Cell *C;
#define U(k)  C[k].u
#define V(k)  C[k].v
#define US(k) C[k].us
#define VS(k) C[k].vs
#define P(k)  C[k].p
#define PN(k) C[k].pn
#define LAYOUT_NAME "AoS"
#else
static double *restrict u_, *restrict v_, *restrict us_, *restrict vs_,
              *restrict p_, *restrict pn_;
#define U(k)  u_[k]
#define V(k)  v_[k]
#define US(k) us_[k]
#define VS(k) vs_[k]
#define P(k)  p_[k]
#define PN(k) pn_[k]
#ifdef PAD
#define LAYOUT_NAME "SoA+pad"
#else
#define LAYOUT_NAME "SoA"
#endif
#endif

/* ---------------- loop order ---------------- */
#ifdef ORDER_JI
#define FOR_INTERIOR for (size_t j = 1; j < N - 1; ++j) for (size_t i = 1; i < N - 1; ++i)
#define ORDER_NAME "j-i (strided)"
#else
#define FOR_INTERIOR for (size_t i = 1; i < N - 1; ++i) for (size_t j = 1; j < N - 1; ++j)
#define ORDER_NAME "i-j (contiguous)"
#endif

static double now(void) {
    struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec * 1e-9;
}

static void *xalloc(size_t bytes) {
    void *ptr = NULL;
    if (posix_memalign(&ptr, 64, bytes)) { perror("posix_memalign"); exit(1); }
    return ptr;
}

int main(int argc, char **argv) {
    N = (argc > 1) ? strtoul(argv[1], NULL, 10) : 256;
    const size_t cells = N * N;
    /* default: keep total work ~constant (~6e8 cell-updates) across N */
    int steps = (argc > 2) ? atoi(argv[2]) : (int)(6e8 / (42.0 * cells));
    if (steps < 1) steps = 1;

    const double dx = 1.0 / (N - 1), dy = dx, rho = 1.0, nu = 0.01;
    const double dt = fmin(0.2 * dx * dx / nu, 0.2 * dx);   /* stable for every N */

#ifdef LAYOUT_AOS
    C = xalloc(cells * sizeof(Cell));
#elif defined(PAD)
    /* One pool; each array starts a few cache lines after the previous one ends,
       so the 6 fields never sit at the same offset modulo a large power of two. */
    {
        const size_t stride = cells + 8 * 5;          /* +5 cache lines per array */
        double *pool = xalloc(6 * stride * sizeof(double));
        u_ = pool; v_ = pool + stride; us_ = pool + 2 * stride;
        vs_ = pool + 3 * stride; p_ = pool + 4 * stride; pn_ = pool + 5 * stride;
    }
#else
    u_ = xalloc(cells * 8); v_ = xalloc(cells * 8); us_ = xalloc(cells * 8);
    vs_ = xalloc(cells * 8); p_ = xalloc(cells * 8); pn_ = xalloc(cells * 8);
#endif
    for (size_t k = 0; k < cells; ++k) U(k) = V(k) = US(k) = VS(k) = P(k) = PN(k) = 0.0;

    const double inv2dx = 1.0 / (2 * dx), inv2dy = 1.0 / (2 * dy);
    const double idx2 = 1.0 / (dx * dx), idy2 = 1.0 / (dy * dy);
    const double denom = 2.0 * (idx2 + idy2);

    double t0 = now();
    for (int s = 0; s < steps; ++s) {
        /* 1. boundary conditions: lid at i = N-1 moves with u = 1 */
        for (size_t j = 0; j < N; ++j) {
            U(IDX(0, j)) = 0; U(IDX(N - 1, j)) = 1; V(IDX(0, j)) = 0; V(IDX(N - 1, j)) = 0;
        }
        for (size_t i = 0; i < N; ++i) {
            U(IDX(i, 0)) = 0; U(IDX(i, N - 1)) = 0; V(IDX(i, 0)) = 0; V(IDX(i, N - 1)) = 0;
        }

        /* 2. predictor: convection + diffusion */
        FOR_INTERIOR {
            size_t c = IDX(i, j), up = c + N, dn = c - N, rt = c + 1, lt = c - 1;
            double dudx = (U(up) - U(dn)) * inv2dx, dudy = (U(rt) - U(lt)) * inv2dy;
            double dvdx = (V(up) - V(dn)) * inv2dx, dvdy = (V(rt) - V(lt)) * inv2dy;
            double lapu = (U(up) - 2 * U(c) + U(dn)) * idx2 + (U(rt) - 2 * U(c) + U(lt)) * idy2;
            double lapv = (V(up) - 2 * V(c) + V(dn)) * idx2 + (V(rt) - 2 * V(c) + V(lt)) * idy2;
            US(c) = U(c) + dt * (-(U(c) * dudx + V(c) * dudy) + nu * lapu);
            VS(c) = V(c) + dt * (-(U(c) * dvdx + V(c) * dvdy) + nu * lapv);
        }

        /* 3. pressure Poisson, Jacobi with ping-pong p -> pn -> p */
        for (int it = 0; it < PRESSURE_ITERS; it += 2) {
            for (int half = 0; half < 2; ++half) {
#define JACOBI(SRC, DST)                                                               \
                for (size_t j = 0; j < N; ++j) {                                       \
                    SRC(IDX(0, j)) = SRC(IDX(1, j)); SRC(IDX(N-1, j)) = SRC(IDX(N-2, j)); } \
                for (size_t i = 0; i < N; ++i) {                                       \
                    SRC(IDX(i, 0)) = SRC(IDX(i, 1)); SRC(IDX(i, N-1)) = SRC(IDX(i, N-2)); } \
                FOR_INTERIOR {                                                         \
                    size_t c = IDX(i, j), up = c + N, dn = c - N, rt = c + 1, lt = c - 1; \
                    double div = (US(up) - US(dn)) * inv2dx + (VS(rt) - VS(lt)) * inv2dy; \
                    DST(c) = ((SRC(up) + SRC(dn)) * idx2 + (SRC(rt) + SRC(lt)) * idy2     \
                              - (rho / dt) * div) / denom;                             \
                }
                if (half == 0) { JACOBI(P, PN) } else { JACOBI(PN, P) }
            }
        }

        /* 4. projection: make velocity divergence-free */
        FOR_INTERIOR {
            size_t c = IDX(i, j), up = c + N, dn = c - N, rt = c + 1, lt = c - 1;
            U(c) = US(c) - (dt / rho) * (P(up) - P(dn)) * inv2dx;
            V(c) = VS(c) - (dt / rho) * (P(rt) - P(lt)) * inv2dy;
        }
    }
    double el = now() - t0;

    double su = 0, sp = 0;
    for (size_t k = 0; k < cells; ++k) { su += fabs(U(k)); sp += fabs(P(k)); }

    double updates = (double)steps * (PRESSURE_ITERS + 2) * (double)(N - 2) * (N - 2);
    double mb = cells * 6.0 * 8.0 / (1024.0 * 1024.0);
    /* CSV: layout,order,N,working_set_MB,steps,seconds,ns_per_cell_update,checksum_u,checksum_p */
    printf("%s,%s,%zu,%.2f,%d,%.4f,%.4f,%.10e,%.10e\n", LAYOUT_NAME, ORDER_NAME, N, mb, steps,
           el, el * 1e9 / updates, su, sp);
    return 0;
}
