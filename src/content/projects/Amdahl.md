---
title: "Amdahl's Law, Measured: Why Twice the Cores Is Not Half the Time"
date: 2026-11-22
tech: ["C", "OpenMP", "MPI", "HPC"]
excerpt: "Amdahl's law predicts the speedup of a program with a serial part, and on compute-bound work it is accurate to a few percent. A memory-bound kernel with no serial part at all stops at about 1.2x, and a CFD halo-exchange code shows 55% strong-scaling efficiency on four. Measuring what limits parallel speedup on a real machine."
---

## Technical Index
1. [**Amdahl and Gustafson**](#laws) — Two laws for two questions
2. [**The experiments**](#experiment) — OpenMP kernels and an MPI solver
3. [**Amdahl, verified**](#amdahl) — Compute-bound work with a serial part
4. [**The limit Amdahl does not include**](#memory) — Memory bandwidth
5. [**A CFD code: strong and weak scaling**](#mpi) — The halo-exchange solver
6. [**What this means in practice**](#practice) — Cores, cells per core and expectations

<a name="laws"></a>
## 1. Amdahl and Gustafson

**Amdahl's law** (1967) answers the question: if a fraction $s$ of a program is serial and the rest can be spread over $P$ processors, how much faster does a *fixed* problem run?

$$
S(P) = \frac{1}{s + \dfrac{1-s}{P}} \;\le\; \frac{1}{s}
$$

Even with infinitely many processors, the speedup cannot exceed $1/s$: with 5% of serial code, never more than 20 times. This is **strong scaling**: same problem, more processors.

**Gustafson's law** (1988) asks a different question: if the problem grows with the number of processors, does the time stay the same? In CFD this is the usual way to use a bigger machine, a bigger mesh. Its measure is the **weak-scaling efficiency**: the time on one processor with one share of the problem, divided by the time on $P$ processors with $P$ shares.

Both laws count only computation and serial code. A real machine adds two more limits: communication between processes, and the hardware that cores share. This article measures all of them on a 4-core laptop (Intel Core i7-1165G7, 8 threads with hyperthreading).

<a name="experiment"></a>
## 2. The experiments

- [amdahl.c](/assets/projects/Amdahl/amdahl.c) splits a fixed amount of work into 64 items. A fraction $s$ of them runs serially, the rest in parallel with OpenMP on 1 to 8 threads. Two kinds of item:
  - **compute-bound**: a chain of floating-point operations on values in registers, with no memory traffic;
  - **memory-bound**: a pass over large arrays ($a_i = b_i + 0.5\,c_i$, 384 MB in total, far larger than the caches), so every byte comes from RAM.
- [halo.c](/assets/projects/Amdahl/halo.c), the MPI halo-exchange solver of the [article on moving data between cores](/projects/data-between-cores/): a mesh split into subdomains (with the METIS-like partitioner), halo exchange, a sweep over the cells and a global reduction per iteration, as in a CFD code. Strong scaling on a 931 000-cell mesh with 1 to 4 processes; weak scaling with 233 000 cells per process. Each configuration was run 3 times and the median is reported.

Times are the best of 5 repetitions after a warm-up run, because the CPU changes its clock frequency with load and temperature.

<a name="amdahl"></a>
## 3. Amdahl, verified

![Compute-bound speedup vs threads, measured and predicted](/assets/projects/Amdahl/amdahl_compute.svg)

| Serial fraction | 2 threads | 4 threads | Amdahl (4) | 8 threads | Amdahl (8) |
|---|---|---|---|---|---|
| 0% | 1.84 | 3.90 | 4.00 | 7.63 | 8.00 |
| 5% | 1.64 | 3.29 | 3.48 | 5.53 | 5.93 |
| 20% | 1.64 | 2.46 | 2.50 | 3.18 | 3.33 |
| 50% | 1.33 | 1.53 | 1.60 | 1.63 | 1.78 |

*Speedup relative to the best single-thread time. Threads 5–8 run on the hyperthreads of the 4 physical cores.*

On compute-bound work Amdahl's law is accurate to within a few percent. With 20% serial code, 4 cores give 2.46 against 2.50 predicted; with 50%, the speedup is stuck around 1.6 whatever the number of threads. The small shortfall is the clock frequency: one active core runs faster than four.

Hyperthreading looks surprisingly effective here (7.6 times on 4 cores) because each work item is a chain of dependent operations: a single thread spends most of its time waiting for the previous result, and a second thread on the same core fills those gaps. On code that already keeps the core busy, hyperthreading brings much less.

<a name="memory"></a>
## 4. The limit Amdahl does not include

Now the same experiment with no serial part at all, but with work that streams through memory:

![Compute-bound vs memory-bound speedup, no serial part](/assets/projects/Amdahl/amdahl_memory.svg)

| Threads | Compute-bound speedup | Memory-bound speedup |
|---|---|---|
| 2 | 1.84 | 1.19 |
| 4 | 3.90 | 1.20 |
| 8 | 7.63 | 1.01 |

The memory-bound kernel, which according to Amdahl should scale perfectly, gains **about 20% on two cores and nothing after that**. One core alone already uses most of the memory bandwidth of this laptop; more cores only compete for it. In terms of Amdahl's law it behaves as if nearly 80% of it were serial, although not a single line of it is.

This is the situation of most CFD solvers. As shown in the [data locality](/projects/data-locality-cfd/) and [compiler and precision](/projects/openfoam-socket-optimization/) articles, a finite-volume solver spends its time moving data. Server processors have far more memory channels per core than a laptop, so they scale much further, but the principle is the same: **cores are added faster than bandwidth**.

<a name="mpi"></a>
## 5. A CFD code: strong and weak scaling

![Strong and weak scaling of the halo-exchange solver](/assets/projects/Amdahl/mpi_scaling.svg)

| Processes | Strong (ms) | Speedup | Weak (ms) | Efficiency | Exchange (ms) |
|---|---|---|---|---|---|
| 1 | 3.18 | 0.87 | 0.77 | 100% | – |
| 2 | 1.87 | 1.48 | 0.83 | 92% | 0.01 |
| 3 | 1.35 | 2.05 | 0.97 | 79% | 0.01–0.02 |
| 4 | 1.27 | 2.18 | 1.27 | 60% | 0.02 |

*Times per iteration in milliseconds. Strong scaling: 931 000 cells, speedup relative to the serial code (2.78 ms per iteration). Weak scaling: 233 000 cells per process. METIS-like decomposition, median of 3 runs.*

- **The parallel code is slower than the serial one on 1 process** (0.87): it carries the bookkeeping of halos and global reductions even when there is nothing to exchange.
- **Strong scaling reaches 2.18 on 4 cores**, an efficiency of 55%.
- **Weak scaling loses 40% of its efficiency at 4 processes**, although every process does the same amount of work.
- **Communication is not the reason.** The halo exchange takes 0.01–0.02 ms per iteration, about 1% of the time. The loss comes from the cores sharing the memory bus, exactly as in the memory-bound kernel, plus the lower clock with all cores active.

If we read these numbers through Amdahl's law (the *Karp–Flatt metric*, $e = (1/S - 1/P)/(1 - 1/P)$), the strong-scaling results correspond to an "effective serial fraction" of 23–35%. The code has essentially no serial part: what Amdahl's law calls serial here is the shared memory bandwidth. The metric is useful precisely because it makes such hidden limits visible: an effective serial fraction that grows with $P$ points to contention or communication, not to code.

<a name="practice"></a>
## 6. What this means in practice

- **Measure before buying cores.** For a memory-bound CFD solver, the useful number of cores per node is set by the memory bandwidth, not by the core count. On many servers the time per iteration stops improving well before all cores are used; sometimes using fewer cores per node and more nodes is faster.
- **Strong scaling has a sweet spot.** As cells per process fall, communication and synchronization grow relative to computation (see the [halo exchange](/projects/data-between-cores/) and [multigrid](/projects/multigrid/) articles). There is a minimum number of cells per core below which adding processes stops helping; it depends on the code, the mesh and the network, so measure it on your own case.
- **Weak scaling is the more natural use of a cluster in CFD.** Bigger machines are usually used for bigger meshes, and weak-scaling efficiency is the number to look at. Here it fell to 60% on 4 laptop cores, but because of shared memory bandwidth, not communication; on a cluster, where each node brings its own memory, weak scaling across nodes is usually much better than within one.
- **Amdahl's law is still the right first question.** Any truly serial step (file I/O, a serial mesh operation, a coarse multigrid level solved on one process) caps the speedup at $1/s$, and becomes the bottleneck as the rest scales. Profiling to find it is always worth it.
- **Hyperthreading rarely helps a CFD solver** in the way it helped the compute-bound kernel here; for memory-bound codes, one process per physical core is the usual choice.

### Reproduce it
All the files are available: [amdahl.c](/assets/projects/Amdahl/amdahl.c), [halo.c](/assets/projects/Amdahl/halo.c) and [decompose.c](/assets/projects/Amdahl/decompose.c), [plot_am.py](/assets/projects/Amdahl/plot_am.py) with [svgchart.py](/assets/projects/Amdahl/svgchart.py), and the data [amdahl.csv](/assets/projects/Amdahl/amdahl.csv), [mpi_scaling.csv](/assets/projects/Amdahl/mpi_scaling.csv).

```bash
gcc -O3 -march=native -fopenmp amdahl.c -o amdahl -lm
mkdir -p out && ./amdahl out                              # ~3 min
mpicc -O3 -march=native -DNX=1600 -DNY=800 -DITERS=200 halo.c -o halo_big -lm
for p in 1 2 3 4; do mpirun -np $p --bind-to core ./halo_big; done      # strong scaling
mpicc -O3 -march=native -DNX=800 -DNY=400 -DITERS=200 halo.c -o halo_w1 -lm
mpirun -np 1 --bind-to core ./halo_w1                     # weak scaling: also 1132x566 (2), 1386x694 (3), 1600x800 (4)
mkdir -p figures && python3 plot_am.py out figures
```
