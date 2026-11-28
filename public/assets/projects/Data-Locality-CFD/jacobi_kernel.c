/* The Jacobi pressure sweep of locality.c, isolated in one function so that
 * the compiler output is easy to read. One row of the grid is updated.
 *   gcc -O3 -march=native -S jacobi_kernel.c -o jacobi_soa.s
 *   gcc -O3 -march=native -S -DAOS jacobi_kernel.c -o jacobi_aos.s
 * Add -fopt-info-vec to see which loop is vectorized and with which width. */
#include <stddef.h>
#ifdef AOS
typedef struct { double u, v, us, vs, p, pn; } Cell;
#define P(k) C[k].p
#define PN(k) C[k].pn
#define US(k) C[k].us
#define VS(k) C[k].vs
void jacobi_row(Cell *C, size_t N, size_t i, double inv2dx, double inv2dy,
  double idx2, double idy2, double rho_dt, double denom) {
#else
#define P(k) p[k]
#define PN(k) pn[k]
#define US(k) us[k]
#define VS(k) vs[k]
void jacobi_row(double *restrict us, double *restrict vs, double *restrict p, double *restrict pn,
  size_t N, size_t i, double inv2dx, double inv2dy, double idx2, double idy2, double rho_dt, double denom) {
#endif
  for (size_t j = 1; j < N - 1; ++j) {
    size_t c = i*N + j, up = c + N, dn = c - N, rt = c + 1, lt = c - 1;
    double div = (US(up) - US(dn)) * inv2dx + (VS(rt) - VS(lt)) * inv2dy;
    PN(c) = ((P(up) + P(dn)) * idx2 + (P(rt) + P(lt)) * idy2 - rho_dt * div) / denom;
  }
}
