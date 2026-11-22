/**
 * MPI ping-pong between two ranks: one-way latency and bandwidth vs message size.
 * Build: mpicc -O3 pingpong.c -o pingpong
 * Run:   mpirun -np 2 --cpu-set 0,1 --bind-to hwthread ./pingpong   (two physical cores)
 *        mpirun -np 2 --cpu-set 0,4 --bind-to hwthread ./pingpong   (same core, two hyperthreads)
 * Output CSV: bytes,one_way_us,bandwidth_GBps
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2) { if (!rank) fprintf(stderr, "run with exactly 2 ranks\n"); MPI_Finalize(); return 1; }

    const size_t maxb = 64u << 20;                  /* 64 MB */
    char *buf = malloc(maxb);
    memset(buf, rank + 1, maxb);
    if (!rank) printf("bytes,one_way_us,bandwidth_GBps\n");

    for (size_t b = 8; b <= maxb; b *= 2) {
        int reps = (int)(2e8 / (b + 1000));         /* ~same total traffic per size */
        if (reps > 20000) reps = 20000;
        if (reps < 10) reps = 10;
        double best = 1e30;
        for (int trial = 0; trial < 3; trial++) {    /* best of 3 */
            MPI_Barrier(MPI_COMM_WORLD);
            double t0 = MPI_Wtime();
            for (int r = 0; r < reps; r++) {
                if (!rank) {
                    MPI_Send(buf, (int)b, MPI_BYTE, 1, 0, MPI_COMM_WORLD);
                    MPI_Recv(buf, (int)b, MPI_BYTE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    MPI_Recv(buf, (int)b, MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(buf, (int)b, MPI_BYTE, 0, 0, MPI_COMM_WORLD);
                }
            }
            double t = (MPI_Wtime() - t0) / (2.0 * reps);   /* one way */
            if (t < best) best = t;
        }
        if (!rank) printf("%zu,%.3f,%.3f\n", b, best * 1e6, b / best / 1e9);
    }
    free(buf);
    MPI_Finalize();
    return 0;
}
