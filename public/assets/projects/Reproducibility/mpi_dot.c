/**
 * The same dot product with 1, 2, 3, 4... MPI processes. Each process sums its
 * own block, MPI_Allreduce adds the partial sums. Mathematically identical,
 * bitwise different.
 * Build: mpicc -O2 mpi_dot.c -o mpi_dot -lm
 * Run:   for p in 1 2 3 4; do mpirun -np $p ./mpi_dot; done
 */
#include <mpi.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);
    const long n = 1000000;
    long lo = n * rank / size, hi = n * (rank + 1) / size;
    double s = 0;
    for (long i = lo; i < hi; i++) { double a = sin(0.001 * i) + 1e-3 * cos(7.0 * i), b = cos(0.002 * i); s += a * b; }
    double tot;
    MPI_Allreduce(&s, &tot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (!rank) printf("%d,%.17g,%a\n", size, tot, tot);
    MPI_Finalize();
    return 0;
}
