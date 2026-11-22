/**
 * Amdahl's law, measured, and the limit it does not include.
 *
 * A fixed amount of work is split into a serial part (fraction s) and a part
 * that runs in parallel with OpenMP on P threads.
 *   kind = compute : each work item is a chain of floating-point operations on
 *                    registers (no memory traffic): the ideal case for Amdahl
 *   kind = memory  : each work item streams through large arrays (a triad,
 *                    a[i] = b[i] + s c[i]): limited by memory bandwidth
 * Output: time and speedup for s = 0, 5%, 20%, 50% and P = 1..8.
 *
 * Build: gcc -O3 -march=native -fopenmp amdahl.c -o amdahl -lm
 * Run:   ./amdahl outdir
 */
#define _DEFAULT_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define NITEMS 64          /* work items */
#define MEMN (1 << 22)     /* doubles per array per item group: 32 MB per array */

static double sink;
/* compute-bound work item: dependent multiply-adds on registers */
static double work_compute(int item) {
    double x = 1.0 + item * 1e-3, y = 0.5;
    for (long k = 0; k < 20000000; k++) { x = x * 0.9999999 + y; y = y * 1.0000001 - 1e-9; }
    return x + y;
}
/* memory-bound work item: triad over a slice of large arrays */
static double *A, *B, *Cc;
static double work_memory(int item) {
    long n = (long)MEMN * 4 / NITEMS;            /* each item owns a contiguous block (2 MB per array) */
    double *a = A + item * n, *b = B + item * n, *c = Cc + item * n;
    for (long k = 0; k < n; k++) a[k] = b[k] + 0.5 * c[k];   /* one pass: no reuse, data always from RAM */
    return a[0];
}

int main(int argc, char **argv) {
    const char *dir = argc > 1 ? argv[1] : ".";
    char path[512]; snprintf(path, sizeof path, "%s/amdahl.csv", dir);
    FILE *f = fopen(path, "w"); fprintf(f, "kind,serial_fraction,threads,seconds,speedup\n");
    A = malloc(MEMN * 4L * 8); B = malloc(MEMN * 4L * 8); Cc = malloc(MEMN * 4L * 8);
    #pragma omp parallel for
    for (long i = 0; i < MEMN * 4L; i++) { A[i] = 0; B[i] = 1; Cc[i] = 2; }
    double sf[] = {0, 0.05, 0.2, 0.5};
    for (int kind = 0; kind < 2; kind++)
        for (int q = 0; q < 4; q++) {
            int nser = (int)lround(sf[q] * NITEMS);
            double t1 = 0;
            for (int P = 1; P <= 8; P++) {
                double best = 1e30;
                for (int rep = 0; rep < 6; rep++) {          /* first repetition = warm-up (frequency ramp), best of the rest */
                    double t0 = omp_get_wtime(), acc = 0;
                    for (int R = 0; R < (kind ? 20 : 1); R++) {   /* memory kind: stream the 384 MB of arrays 20 times */
                        for (int it = 0; it < nser; it++) acc += kind ? work_memory(it) : work_compute(it);   /* serial part */
                        #pragma omp parallel for num_threads(P) reduction(+:acc) schedule(static)
                        for (int it = nser; it < NITEMS; it++) acc += kind ? work_memory(it) : work_compute(it);
                    }
                    sink += acc;
                    double t = omp_get_wtime() - t0; if (rep > 0 && t < best) best = t;
                }
                if (P == 1) t1 = best;
                printf("%-8s s = %4.0f%%  P = %d  %7.3f s  speedup %.2f\n", kind ? "memory" : "compute", 100 * sf[q], P, best, t1 / best);
                fprintf(f, "%s,%g,%d,%.5f,%.4f\n", kind ? "memory" : "compute", sf[q], P, best, t1 / best);
                fflush(stdout);
            }
        }
    fclose(f);
    printf("(checksum %g)\n", sink);
    return 0;
}
