/**
 * Two things -ffast-math is allowed to break.
 *  1. Compensated (Kahan) summation: the compiler may simplify the
 *     correction term away, because algebraically it is zero.
 *  2. NaN detection: -ffast-math implies -ffinite-math-only, so the compiler
 *     may assume NaN never happens and fold isnan() to false.
 * Build twice and compare:
 *   gcc -O3 fastmath_traps.c -o traps_safe -lm
 *   gcc -O3 -ffast-math fastmath_traps.c -o traps_fast -lm
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

__attribute__((noinline)) static double naive_sum(const double *a, long n) {
    double s = 0; for (long i = 0; i < n; i++) s += a[i]; return s;
}
__attribute__((noinline)) static double kahan_sum(const double *a, long n) {
    double s = 0, c = 0;
    for (long i = 0; i < n; i++) { double y = a[i] - c, t = s + y; c = (t - s) - y; s = t; }
    return s;
}
__attribute__((noinline)) static int diverged(double residual) { return isnan(residual) || isinf(residual); }

int main(int argc, char **argv) {
    long n = 10000000;
    double *a = malloc(n * sizeof(double));
    /* one big value followed by many small ones: the classic case for compensation */
    a[0] = 1e8; for (long i = 1; i < n; i++) a[i] = 1e-8 * (1 + (i % 7));
    double exact = 1e8;                       /* exact sum of the small terms, computed carefully */
    { long double e = 0; for (long i = 1; i < n; i++) e += (long double)a[i]; exact += (double)e; }
    double sn = naive_sum(a, n), sk = kahan_sum(a, n);
    printf("sum:  exact %.10f\n      naive %.10f (error %.2e)\n      kahan %.10f (error %.2e)\n",
           exact, sn, fabs(sn - exact), sk, fabs(sk - exact));
    double zero = argc > 100 ? 1.0 : 0.0;     /* not known at compile time */
    double r = zero / zero;                   /* NaN: a "diverged" residual */
    printf("residual = %f -> diverged() says %s\n", r, diverged(r) ? "YES" : "no (!)");
    free(a);
    return 0;
}
