---
title: "Moving Data Between Cores"
date: 2026-11-22
tech: ["C", "MPI", "CFD", "Linux"]
excerpt: "A parallel CFD run is as fast as its slowest core and as cheap as the data it exchanges. Measuring MPI messages, halo exchange and domain decomposition on a 4-core laptop."
---

## Technical Index
1. [**Why cores have to talk**](#why) — Domain decomposition and halo cells
2. [**The price of a message**](#pingpong) — Latency and bandwidth between two cores
3. [**Where the messages come from**](#decomp) — Four ways to split a mesh
4. [**A real halo exchange**](#halo) — Compute, exchange and waiting, measured
5. [**What this means for real CFD**](#realcfd) — From a laptop to `decomposeParDict`

In the [previous article](/projects/data-locality-cfd/) I looked at how data moves *inside* one core: caches, cache lines, layouts. Here the question is the next one: what happens when the data has to move *between* cores.

<a name="why"></a>
## 1. Why cores have to talk

A parallel CFD solver (OpenFOAM, Fluent, SU2…) splits the mesh into subdomains, one per process. Each process updates only its own cells, but the discretization couples every cell to its neighbours: the flux through a face needs the values on both sides. When the two sides of a face belong to different processes, someone has to send that value across.

The standard solution is the **halo** (or *ghost cells*): every process keeps a read-only copy of the neighbouring cells owned by other processes. At every iteration, before computing, it receives fresh values for its halo and sends the values its neighbours need. In OpenFOAM these are the `processor` boundary patches that `decomposePar` creates.

So every iteration of a parallel solver is made of three parts:

1. **compute**: update the owned cells, using owned and halo values;
2. **exchange**: send and receive halo values (MPI messages);
3. **synchronize**: global operations like the residual norm or a dot product in the linear solver (`MPI_Allreduce`) need *every* process, so the fast ones wait for the slowest one.

Only the first part is useful work. The other two are the price of running in parallel, and both depend on how the mesh was cut.

<a name="pingpong"></a>
## 2. The price of a message

First, the raw cost. [pingpong.c](/assets/projects/Data-Between-Cores/pingpong.c) is the classic MPI ping-pong: two processes, pinned on two different physical cores, send a message back and forth thousands of times, for sizes from 8 bytes to 64 MB. On a single machine Open MPI does not use a network: the message is copied through shared memory.

### Test configuration
- CPU: Intel Core i7-1165G7, 4 physical cores (8 threads), L2 1.25 MB per core, L3 12 MB shared
- RAM: 16 GB
- Open MPI 5.0, GCC 16.2, `-O3 -march=native`, processes pinned with `--bind-to core`
- OS: CachyOS Linux, default `powersave` governor: single numbers are noisy, so every result below is a median or a best of several runs

![One-way time of an MPI message between two cores](/assets/projects/Data-Between-Cores/pingpong_latency.svg)

![Bandwidth between two cores vs message size](/assets/projects/Data-Between-Cores/pingpong_bandwidth.svg)

| Message size | One-way time (µs) | Bandwidth (GB/s) |
|---|---|---|
| 8 B | 0.29 | 0.03 |
| 1 KB | 0.65 | 1.6 |
| 4 KB | 2.04 | 2.0 |
| 64 KB | 8.2 | 8.0 |
| 4 MB | 259 | 16.2 |
| 64 MB | 8 124 | 8.3 |

*One-way time: half of a round trip, in microseconds. Bandwidth: message size divided by the one-way time, in gigabytes per second.*

Three things stand out:

- **Small messages cost a fixed price.** Up to about 64 bytes the time does not depend on size: about **0.3 µs**, roughly a thousand clock cycles, just to deliver a message from one core to the next. That is the *latency*. Sending 8 bytes uses 0.2% of the available bandwidth.
- **There is a step at 4 KB.** Below it, MPI copies small messages eagerly into a shared buffer; above it, it switches to a different protocol (a handshake, then the copy), which costs an extra round trip. Many small halo messages pay the latency over and over; one bigger message pays it once.
- **Bandwidth falls after the L3 cache.** It peaks at **16 GB/s** for 4 MB messages, then drops to **8–10 GB/s** from 8 MB on. A message travels through memory: as long as sender buffer, receiver buffer and the copy fit in the shared L3, the copy runs at cache speed; beyond that, it goes through RAM. Data locality again, this time between cores.

And these are the *best* numbers you will ever get: two cores on the same chip. Between two machines over a network, the latency is typically several times higher.

<a name="decomp"></a>
## 3. Where the messages come from

How many values have to be exchanged, and with how many neighbours, is decided before the run, when the mesh is cut. To see the effect I wrote [decompose.c](/assets/projects/Data-Between-Cores/decompose.c), which builds a test mesh and splits it into 4 subdomains in four ways.

- **Mesh:** an L-shaped domain with a cylinder, 1600 × 800 bounding box, **931 032 cells** and 1 859 280 internal faces. It is treated only as a graph (cells plus face neighbours), like a solver reading an unstructured mesh, and the cells are numbered in random order, as a mesh generator could produce them.
- **Methods:** four ways to split the mesh. **The asterisk (scotch\*, metis\*) marks my own simplified versions, not the real libraries.** Scotch and METIS are real partitioning libraries with much more sophisticated algorithms; the starred versions are a few dozen lines of C each that keep only the core idea of each one, so their numbers say something about the *approach*, not about the performance of the actual tools.
  - **random**: every cell goes to a random subdomain.
  - **simple**: cut the bounding box into 2 × 2 equal blocks, like OpenFOAM's `simple`. The only method that uses coordinates.
  - **scotch\***: recursive bisection of the graph. Start a BFS from a far-away cell, cut the BFS order in half, repeat on each half, then refine the boundary.
  - **metis\***: direct 4-way split. Place 4 seeds far apart, let the smallest subdomain always grow next, then refine the boundary.

![The four decompositions of the L-shaped mesh into 4 subdomains](/assets/projects/Data-Between-Cores/decomposition_4.png)

**Communication**

| Method | Cut faces | Cut (% faces) | Halo (KB/iter) | Max nbrs |
|---|---|---|---|---|
| random | 1 393 754 | 75.0 | 21 777 | 3 |
| simple | 1 008 | 0.05 | 16 | 2 |
| scotch\* | 3 050 | 0.16 | 48 | 3 |
| metis\* | 2 989 | 0.16 | 47 | 3 |

*Cut faces: faces whose two cells belong to different subdomains (count and % of all faces). Halo: values exchanged by all subdomains together per iteration for one field in double precision, one 8-byte value per cut face in each direction, in KB. Max nbrs: the largest number of other subdomains one subdomain exchanges with.*

**Load balance**

| Method | Largest (%) | Empty | Time (ms) |
|---|---|---|---|
| random | +0.3 | 0 | 5 |
| simple | +37.5 | 1 | 0.4 |
| scotch\* | +0.0 | 0 | 378 |
| metis\* | +0.0 | 0 | 819 |

*Largest: cells in the biggest subdomain, in % above the perfect share (931 032 / 4). Empty: subdomains with no cells. Time: time to compute the decomposition, in milliseconds.*

On paper:
- **random** is perfectly balanced but cuts three faces out of four: almost the whole field would have to travel at every iteration.
- **simple** cuts the fewest faces, but it cuts the *bounding box*, not the mesh. On the L-shape the top-right block is empty: one process gets nothing, and another one gets 37.5% more cells than it should.
- The **graph methods** are balanced and cut very little, without knowing anything about geometry. The toy versions also show their limits: the scotch\* cuts follow diagonal BFS fronts, and one of its subdomains ends up in two pieces. The real libraries avoid this with *multilevel* partitioning: they coarsen the graph many times, split the small graph and refine the cut while projecting it back.

The question is how these numbers turn into time.

<a name="halo"></a>
## 4. A real halo exchange

[halo.c](/assets/projects/Data-Between-Cores/halo.c) runs the four decompositions with real MPI processes, one per core. It is not a CFD solver, but it has exactly the structure of one. Every process keeps its own cells (numbered in BFS order, for good locality) plus a halo, and iterates:

```c
/* 1-3: exchange: receive the halo, send what the neighbours need */
for (q = 0; q < P; q++) if (recv_count[q])
    MPI_Irecv(rbuf + roff[q], recv_count[q], MPI_DOUBLE, q, it, MPI_COMM_WORLD, &req[nr++]);
for (s = 0; s < nsend; s++) sbuf[s] = x[send_list[s]];          /* pack   */
for (q = 0; q < P; q++) if (send_count[q])
    MPI_Isend(sbuf + soff[q], send_count[q], MPI_DOUBLE, q, it, MPI_COMM_WORLD, &req[nr++]);
MPI_Waitall(nr, req, MPI_STATUSES_IGNORE);
memcpy(x + nown, rbuf, nhalo * sizeof(double));                  /* unpack */

/* 4: compute: every owned cell from its face neighbours (owned or halo) */
for (l = 0; l < nown; l++) {
    double s = 0.0;
    for (e = xadj[l]; e < xadj[l + 1]; e++) s += x[adj[e]];
    y[l] = 0.25 * s + 0.01;
}

/* 5: synchronize: a global residual, like a convergence check */
MPI_Allreduce(&res, &res_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
```

The three phases are timed separately. Since the `Allreduce` cannot finish until every process reaches it, its time is mostly **waiting for the slowest process**: it measures load imbalance. At the end, the result is checked against a serial run on one core: it matches for every method.

I ran it on two meshes, because the answer depends on where the data lives:
- **small**: 800 × 400 box, **232 768 cells**. Each subdomain is about 58 000 cells, a couple of MB: it fits in the caches.
- **big**: the 1600 × 800 mesh above, **931 032 cells**. Each subdomain is about 8 MB, and the four of them together do not fit in the 12 MB L3: the cores share the memory bandwidth.

![Time per iteration on the small mesh](/assets/projects/Data-Between-Cores/halo_small.svg)

![Time per iteration on the big mesh](/assets/projects/Data-Between-Cores/halo_big.svg)

| Mesh | Method | Compute | Exchange | Sync/wait | Iteration | Speedup |
|---|---|---|---|---|---|---|
| small | random | 0.591 | 0.912 | 0.152 | 1.639 | 0.46× |
| small | simple | 0.286 | 0.003 | 0.118 | 0.329 | 2.30× |
| small | scotch\* | 0.239 | 0.005 | 0.054 | 0.297 | 2.92× |
| small | metis\* | 0.213 | 0.005 | 0.043 | 0.253 | 3.00× |
| big | random | 2.846 | 5.979 | 0.599 | 9.374 | 0.36× |
| big | simple | 1.557 | 0.013 | 0.576 | 1.723 | 1.90× |
| big | scotch\* | 1.650 | 0.025 | 0.242 | 1.863 | 1.96× |
| big | metis\* | 1.536 | 0.024 | 0.217 | 1.743 | 2.01× |

*Times in milliseconds per iteration, median of 5 runs. Compute: slowest process. Exchange and sync/wait: average over processes. Iteration: wall time of the slowest process. Speedup: time on 1 core divided by time on 4 cores (1 core: 0.74 ms small, 3.30 ms big).*

### Random: four cores slower than one
With a random split on the big mesh, each process owns 233 000 cells but needs a halo of about **477 000 cells** from the others: twice its own size. Every iteration moves 15 MB between cores, and that alone takes longer than the whole serial iteration (6.0 ms against 3.3 ms; on the small mesh, 0.91 ms against 0.74 ms). Even the compute part is slower than it should be: neighbours are scattered across a huge halo array, so the sweep loses the locality that the BFS numbering gives to the other methods. Result: **0.36×**, four cores are almost three times *slower* than one.

### Simple: no communication, but someone waits
Simple exchanges almost nothing (16 KB per iteration, 0.003–0.013 ms). Its problem is the empty block: one core does nothing, and one has 37.5% more cells than the others. Everybody waits for that one in the `Allreduce`: the sync time is 2–3× higher than with the graph methods. On the small mesh this costs a lot (2.30× vs 3.00×).

### Graph methods: what is left is the hardware
With scotch\* and metis\* the exchange is about **5 µs per iteration** on the small mesh, about 2% of the time, and the load is balanced. On the small mesh they reach **3×** on 4 cores. The missing part is the remaining waiting (the cores do not run at exactly the same speed) and the fact that one core alone runs at a higher turbo clock than four together.

### The big mesh: bandwidth is shared
On the big mesh *every* method stops at about **2×**, and simple looks almost as good as the others. The reason is not communication, which is still around 0.02 ms: the sweep is memory-bound, and the four cores share the same memory bus. Running the graph methods on 2 and 4 cores shows it:

| Mesh | 2 cores | 4 cores |
|---|---|---|
| small (fits in cache) | 1.7× | 3.0× |
| big (memory-bound) | 1.4× | 2.0× |

*Speedup over 1 core, scotch\* and metis\*, median of 3 runs on 2 cores and 5 runs on 4 cores.*

On the big mesh each added core contributes less and less, because it competes for the same bandwidth. This probably also hides part of the imbalance of simple: its idle core leaves more bandwidth to the other three. It is the same wall as in the previous article, seen from the parallel side: **adding cores does not add memory bandwidth**.

<a name="realcfd"></a>
## 5. What this means for real CFD

A 4-core laptop is not a cluster, but the mechanisms are the same, and on a cluster they are usually worse:

- **Use a graph partitioner on real geometries.** In OpenFOAM, `scotch` (the usual default) or `metis` in `decomposeParDict`. `simple` and `hierarchical` are fine on box-shaped domains and a trap on everything else. `decomposePar` prints the number of processor faces and the cell imbalance: they are exactly the numbers in the tables above, and worth reading before launching a long run.
- **Imbalance is paid by everyone.** The run goes at the speed of the slowest process, at every iteration, at every global reduction. A few percent of extra cells on one processor is a few percent on the whole run.
- **Messages have a fixed cost.** With 0.3 µs on-chip and more between nodes, a solver that sends many small messages is limited by latency, not bandwidth. Fewer neighbours per subdomain means fewer messages.
- **Surface grows relative to volume.** Splitting a mesh into more and more subdomains (strong scaling) shrinks the cells per core faster than the halo. At some point exchange and waiting dominate: this is why every case has a sweet spot in cells per core, and why it is worth measuring it instead of guessing.
- **More cores do not mean more memory bandwidth.** On a memory-bound solver, filling every core of a node often gives much less than linear speedup. Sometimes using fewer cores per node, spread over more nodes, is faster.
- **Pin the processes** (`--bind-to core`, or the equivalent of your MPI and scheduler), so that a process and its data stay on the same core and its caches.

The common thread with the previous article: in scientific computing the cost of moving data shapes everything, from a single loop to a whole cluster. Inside a core it is cache lines and strides; between cores it is halos, messages and waiting. In both cases the numerics stay the same, and the time can change by an order of magnitude.

### Reproduce it
All the files are available: [pingpong.c](/assets/projects/Data-Between-Cores/pingpong.c), [halo.c](/assets/projects/Data-Between-Cores/halo.c), [decompose.c](/assets/projects/Data-Between-Cores/decompose.c), [Makefile](/assets/projects/Data-Between-Cores/Makefile), [plot_comm.py](/assets/projects/Data-Between-Cores/plot_comm.py) and the raw data [pp_cores.csv](/assets/projects/Data-Between-Cores/pp_cores.csv), [halo_runs.csv](/assets/projects/Data-Between-Cores/halo_runs.csv).

```bash
make all                                   # needs an MPI implementation (mpicc, mpirun)
mpirun -np 2 --cpu-set 0,1 --bind-to hwthread ./pingpong > pp_cores.csv
mpirun -np 4 --bind-to core ./halo_small   # ~3 s
mpirun -np 4 --bind-to core ./halo_big     # ~10 s
mkdir -p img && ./decompose 4 img          # partition statistics + one .ppm image per method
```

Check which logical CPUs belong to different physical cores (`lscpu -e`) before choosing the `--cpu-set` for the ping-pong.
