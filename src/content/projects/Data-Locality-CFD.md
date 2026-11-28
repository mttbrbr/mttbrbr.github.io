---
title: "Data Locality in CFD"
date: 2026-11-22
tech: ["C", "CFD", "Linux", "Hardware"]
excerpt: "Same equations, same compiler, same machine, up to 20× difference in runtime. A small 2D Navier-Stokes solver shows why scientific codes must be designed around how data sits in memory."
---

## Technical Index
1. [**The problem**](#problem) — Why FLOPs are not the bottleneck
2. [**The test case**](#testcase) — A tiny lid-driven cavity solver in C
3. [**Results**](#results) — Layout, loop order and grid size
4. [**Reading the assembly**](#asm) — What the compiler does with SoA and with AoS
5. [**The power-of-two trap**](#pow2) — When the grid size itself hurts you
6. [**What this means for real CFD**](#realcfd) — From a toy solver to OpenFOAM

<a name="problem"></a>
## 1. The problem

When we talk about the cost of a CFD simulation we usually count cells, iterations and floating-point operations. That mental model is incomplete. A modern core can do several double-precision multiply-adds per clock cycle, but a single access to main memory costs on the order of a hundred nanoseconds, which is hundreds of cycles in which the core does nothing. In between sit the caches: L1, L2 and L3, each one larger and slower than the previous one.

Data is not moved one number at a time: it travels in **cache lines** of 64 bytes, that is 8 doubles. Every time we touch one value, the hardware brings its 7 neighbours along "for free". Whether those neighbours are useful or not depends entirely on how *we* laid out the data and in which order we walk through it. This is **data locality**:

- **spatial locality**: if I use `a[k]`, I will soon use `a[k+1]`;
- **temporal locality**: if I use `a[k]` now, I will use it again soon, while it is still in cache.

Finite-volume and finite-difference solvers are almost always **memory-bound**: a stencil update reads several neighbours, does a handful of operations and writes one value. The arithmetic intensity (FLOPs per byte moved) is low, so the speed of the code is set by how efficiently we feed the core, not by how fast it computes.

<a name="testcase"></a>
## 2. The test case

To make this concrete I took a small 2D incompressible solver I had already written for compiler experiments: a **lid-driven cavity** on a uniform $N \times N$ grid, solved with a fractional step method:

1. predictor step for the velocity (convection + diffusion), giving $u^*, v^*$;
2. pressure Poisson equation $\nabla^2 p = \frac{\rho}{\Delta t}\nabla\cdot\mathbf{u}^*$, solved with 40 Jacobi sweeps;
3. projection, $\mathbf{u} = \mathbf{u}^* - \frac{\Delta t}{\rho}\nabla p$.

Every time step is 42 sweeps over the grid of a 5-point stencil. Six fields are stored: $u, v, u^*, v^*, p, p_{new}$. The time step is scaled with the grid so the scheme stays stable at every size.

The important point: **the numerics never change**. I compiled the same source in four variants, changing only *where the data lives* and *in which order it is visited*.

### Contiguous vs strided: the order of the loops
The grid is two-dimensional, but memory is a one-dimensional list of addresses. The solver stores it **row by row** (row-major): cell `(i, j)` is at position `i*N + j`. So the cells of one row are neighbours in memory, while two cells of the same column are $N$ positions apart. With $N = 4$:

```text
grid (i = row, j = column)       memory address of each cell

        j=0  j=1  j=2  j=3
i=0   [  0    1    2    3 ]      row 0 -> addresses  0  1  2  3
i=1   [  4    5    6    7 ]      row 1 -> addresses  4  5  6  7
i=2   [  8    9   10   11 ]      row 2 -> addresses  8  9 10 11
i=3   [ 12   13   14   15 ]      row 3 -> addresses 12 13 14 15

along a row    (j -> j+1): +1 address
down a column  (i -> i+1): +N addresses
```

The two loop orders visit exactly the same cells, only in a different sequence:

```c
/* contiguous: j in the inner loop, walk along each row      */
for (i = 1; i < N-1; i++)          /* visits 0, 1, 2, 3, 4, ... */
    for (j = 1; j < N-1; j++)      /* next cell = next address  */
        update(i, j);

/* strided: i in the inner loop, walk down each column       */
for (j = 1; j < N-1; j++)          /* visits 0, 4, 8, 12, 1, ...*/
    for (i = 1; i < N-1; i++)      /* next cell = N addresses ahead */
        update(i, j);
```

- **Contiguous**: consecutive iterations read consecutive addresses. One 64-byte cache line holds 8 doubles, so one load from memory serves the next 8 cells, and the hardware prefetcher can easily guess what comes next.
- **Strided**: consecutive iterations are $N$ doubles apart (the *stride*), which is $8N$ bytes, e.g. 16 KB for $N = 2048$. Every iteration lands on a different cache line, uses 8 bytes of it and moves on. The other 7 values in that line are only needed when the loop comes back to the next column, and by then the line has usually been evicted.

Same arithmetic, same result, different path through memory.

### Memory layout: SoA vs AoS
The second variable is how the six fields are grouped:


```c
/* Structure of Arrays (SoA): each field is contiguous */
double *u, *v, *us, *vs, *p, *pn;

/* Array of Structures (AoS): all fields of one cell are together, 48 bytes */
typedef struct { double u, v, us, vs, p, pn; } Cell;
Cell *C;
```

Combining the two choices gives the four variants (the last column is the distance in memory between two consecutively visited cells):

| Variant | Memory layout | Loop order | Step in memory |
|---|---|---|---|
| SoA, contiguous | one array per field | `i` outer, `j` inner | 8 B (next double) |
| AoS, contiguous | one struct per cell | `i` outer, `j` inner | 48 B (next struct) |
| SoA, strided | one array per field | `j` outer, `i` inner | $8N$ B (next row) |
| AoS, strided | one struct per cell | `j` outer, `i` inner | $48N$ B (next row) |

and the kernel, written once with macros so that both layouts compile from the very same lines (here the Jacobi sweep):

```c
FOR_INTERIOR {
    size_t c = IDX(i, j), up = c + N, dn = c - N, rt = c + 1, lt = c - 1;
    double div = (US(up) - US(dn)) * inv2dx + (VS(rt) - VS(lt)) * inv2dy;
    PN(c) = ((P(up) + P(dn)) * idx2 + (P(rt) + P(lt)) * idy2
             - (rho / dt) * div) / denom;
}
```

At the end of every run the solver prints a checksum of $|u|$ and $|p|$: it is **bit-for-bit identical** across all variants for every grid size, so the differences below are purely about memory.

### How the program is organised
`locality.c` is a single file of about 130 lines, and the whole experiment lives in three layers.

**1. Compile-time switches.** `-DLAYOUT_AOS`, `-DORDER_JI` and `-DPAD` choose the variant. The kernels never name an array directly: they use the macros `U(k)`, `V(k)`, `US(k)`, `VS(k)`, `P(k)`, `PN(k)` for field access and `FOR_INTERIOR` for the double loop. The preprocessor expands them differently per variant:

```c
/* SoA */  #define P(k)  p_[k]          /* p_ is its own array of doubles   */
/* AoS */  #define P(k)  C[k].p        /* p sits at byte offset 32 of a 48 B struct */

/* contiguous */ #define FOR_INTERIOR for (i = 1; i < N-1; ++i) for (j = 1; j < N-1; ++j)
/* strided    */ #define FOR_INTERIOR for (j = 1; j < N-1; ++j) for (i = 1; i < N-1; ++i)
```

So there is exactly one copy of the numerical code, and the four binaries cannot differ in anything except memory. This is also why the checksums match bit for bit.

**2. Setup.** The grid size `N` comes from the command line. The number of time steps is computed as $6\cdot10^8 / (42 \cdot N^2)$ so that every run does the same total work, and `dt` scales with the grid so the scheme stays stable. Arrays are allocated with `posix_memalign` on 64-byte boundaries, that is on a cache-line boundary, and then written once with zeros so that the pages really exist in RAM before the clock starts.

**3. The time loop.** Each step runs four stages, and only the clock around this loop is measured:

| Stage | Reads | Writes | Sweeps per step |
|---|---|---|---|
| Boundary conditions | none | edges of `u`, `v` | 1 (edges only) |
| Predictor | `u`, `v` | `us`, `vs` | 1 |
| Jacobi, pressure | `p` or `pn`, `us`, `vs` | `pn` or `p` | 40 |
| Projection | `us`, `vs`, `p` | `u`, `v` | 1 |

Two details matter for what follows. First, Jacobi is 40 of the 42 sweeps, so it is the kernel that decides the result and the one I disassemble below. Second, **no kernel touches all six fields**: the Jacobi sweep uses four of them and the predictor four, so in the AoS layout at least a third of every 48-byte struct is carried along for nothing. The Jacobi sweeps ping-pong between `p` and `pn`: an even number of iterations leaves the final pressure in `p`, and each half-sweep first copies the ghost values of its source array, then runs the stencil. At the end the program sums $|u|$ and $|p|$ over the grid (the checksum) and prints one CSV line with the layout, the loop order, the working set and the time per cell update.

### Test configuration
- CPU: Intel Core i7-1165G7 (4 cores, L1d 48 KB/core, L2 1.25 MB/core, L3 12 MB shared), single thread
- RAM: 16 GB
- Compiler: GCC 16.2, `-O3 -march=native`
- OS: CachyOS Linux
- Grid sizes from $64^2$ to $2896^2$, i.e. a working set from 0.2 MB to 384 MB. The number of time steps is chosen so that every run does the same total work (~$6\cdot10^8$ cell updates). Each point is the median of 3 runs.

The metric is **nanoseconds per cell update**: if the code scaled perfectly, it would be a flat line.

<a name="results"></a>
## 3. Results

![Time per cell update vs working set for the four variants](/assets/projects/Data-Locality-CFD/locality_sweep.svg)

| Grid | Working set (MB) | SoA contig. | AoS contig. | AoS strided | SoA strided |
|---|---|---|---|---|---|
| 64² | 0.2 | 1.38 | 6.56 | 3.72 | 3.06 |
| 128² | 0.8 | 1.17 | 6.64 | 4.00 | 5.43 |
| 256² | 3.0 | 1.18 | 6.61 | 4.73 | 12.49 |
| 512² | 12.0 | 1.43 | 8.32 | 13.80 | 17.73 |
| 724² | 24.0 | 5.60 | 9.39 | 11.96 | 10.12 |
| 1024² | 48.0 | 2.13 | 8.37 | 29.47 | 32.95 |
| 1448² | 96.0 | 5.26 | 7.96 | 12.56 | 10.03 |
| 2048² | 192.0 | 2.32 | 9.00 | 38.89 | 47.41 |
| 2896² | 383.9 | 2.08 | 8.29 | 10.49 | 14.96 |

*Grid: cells. Working set: memory occupied by the six fields, in MB. The four variant columns are the time to update one grid cell once (one stencil sweep), in **nanoseconds**: lower is better.*

### Loop order: up to 20× for swapping two `for`
The most striking result is also the simplest change. Walking the grid along the storage order (SoA, contiguous) costs about **1.2 ns** per update while the grid fits in cache, and only about **2 ns** once it lives in RAM: the hardware prefetcher sees a sequential stream and hides most of the latency.

Swapping the two loops makes every access jump $N$ doubles. Each cache line brought in is used for a single value and then, by the time the loop comes back to its neighbour, it has been evicted. At $256^2$ (3 MB, still inside L3) the strided version is already **10.6× slower**; at $2048^2$ it is **20× slower** (47.4 vs 2.3 ns).

Hardware counters (`perf stat`, $2048^2$, 3 time steps) show exactly where the time goes:

| Counter | SoA contig. | AoS contig. | SoA strided |
|---|---|---|---|
| IPC (instr./cycle) | 1.02 | 0.77 | 0.19 |
| L1 misses (× 10⁹) | 0.47 | 1.20 | 2.47 |
| L3 misses (× 10⁶) | 6.8 | 24.8 | 1176 |

*Whole run. IPC: instructions completed per clock cycle, higher is better. L1 misses: loads not found in the L1 data cache, in billions. L3 misses: loads not found in the last-level cache, i.e. trips to RAM, in millions. Fewer misses is better.*

The strided version misses the last-level cache **172 times more often**. With an IPC of 0.19 the core spends four cycles out of five waiting for memory. There is a second, hidden penalty too: with a stride the compiler can no longer vectorize the loop, so it also executes about 3× more instructions.

### Layout: AoS pays even when everything fits in cache
The AoS version is about **5.6× slower** than SoA even on tiny grids that fit in L2, and around **4×** on large ones. Two reasons:

1. **Wasted bandwidth.** The Jacobi sweep needs only `p`, `us`, `vs` and writes `pn`, but every cache line it loads also contains `u` and `v`, which it never uses. A large fraction of each 64-byte transfer is dead weight.
2. **Vectorization goes wrong.** With SoA the compiler fills 256-bit AVX2 registers with 4 consecutive values of the same field: one load, four useful numbers. With AoS consecutive values of `p` are 48 bytes apart, so GCC (as `-fopt-info-vec` reports) falls back to 128-bit vectors and has to gather the values with shuffles. That code is actually *slower* than not vectorizing at all: recompiling the AoS variant with `-fno-tree-vectorize` brings it from 6.4 to 2.8 ns at $256^2$, and from 7.8 to 5.0 ns at $2048^2$. Even then it stays 2.4× behind SoA.

This also explains a detail of the chart: at small sizes the *strided* AoS variant beats the *contiguous* one. The strided loop cannot be vectorized, so the compiler emits plain scalar code, and since all six fields of a cell share a cache line, one miss brings in everything the stencil needs for that cell. That is the one thing AoS is good at, and it is why "AoS vs SoA" has no universal answer: it depends on which fields each kernel actually reads. The lesson for the compiler, instead, is general: `-O3` is not a guarantee, and the data layout decides what the optimizer can do.

### The cache cliff
Look at the SoA strided curve: 3 ns at 0.2 MB, 5.4 ns at 0.8 MB, 12.5 ns at 3 MB, 17.7 ns at 12 MB. Every time the working set crosses a cache level, the cost per update jumps. The physics did not change, the algorithm did not change, only the size of the problem did. A code profiled on a small test case and then run on the production mesh can behave completely differently.

<a name="asm"></a>
## 4. Reading the assembly

The timings say that AoS is slower, and the `-fopt-info-vec` report says GCC vectorizes it differently. To see what that means I isolated the Jacobi stencil in a small function ([jacobi_kernel.c](/assets/projects/Data-Locality-CFD/jacobi_kernel.c)), one row of the grid, written once for SoA and once for AoS, and compiled both with the same flags as the benchmark (`-O3 -march=native`, GCC 16.2). `make asm` regenerates [jacobi_soa.s](/assets/projects/Data-Locality-CFD/jacobi_soa.s) and [jacobi_aos.s](/assets/projects/Data-Locality-CFD/jacobi_aos.s). The listings below are the inner loops, copied from the compiler output; the comments are mine. Register choices depend on the compiler and the CPU, the structure does not.

The stencil in C, for reference:

```c
div   = (US(up) - US(dn)) * inv2dx + (VS(rt) - VS(lt)) * inv2dy;
PN(c) = ((P(up) + P(dn)) * idx2 + (P(rt) + P(lt)) * idy2 - rho_dt * div) / denom;
```

It needs 8 loads, 13 floating-point operations (including one division) and 1 store per cell.

### SoA: one load per vector of useful numbers
GCC vectorizes the loop with 256-bit `ymm` registers, 4 doubles each, so one iteration updates **4 cells in 21 instructions**:

```asm
.L4:
    vmovupd  (%rdx,%rax), %ymm20        # us[up .. up+3]   4 consecutive doubles
    vmovupd  (%rcx,%rax), %ymm23        # vs[rt .. rt+3]
    vmovupd  0(%r13,%rax), %ymm0        # p[rt .. rt+3]
    vmovupd  (%r15,%rax), %ymm19        # p[up .. up+3]
    vsubpd   (%rsi,%rax), %ymm20, %ymm21   # us[up] - us[dn]   (load folded into the subtraction)
    vsubpd   (%rbx,%rax), %ymm23, %ymm24   # vs[rt] - vs[lt]
    vaddpd   (%r12,%rax), %ymm0,  %ymm0    # p[rt] + p[lt]
    vaddpd   (%r14,%rax), %ymm19, %ymm19   # p[up] + p[dn]
    vmulpd   %ymm11, %ymm21, %ymm22     # * inv2dx
    vmulpd   %ymm10, %ymm24, %ymm25     # * inv2dy
    vmulpd   %ymm12, %ymm0,  %ymm0      # * idy2
    vmulpd   %ymm17, %ymm19, %ymm19     # * idx2
    vaddpd   %ymm25, %ymm22, %ymm26     # div
    vaddpd   %ymm19, %ymm0,  %ymm0      # sum of the two p terms
    vmulpd   %ymm9,  %ymm26, %ymm27     # rho_dt * div
    vsubpd   %ymm27, %ymm0,  %ymm0
    vdivpd   %ymm8,  %ymm0,  %ymm0      # / denom
    vmovupd  %ymm0, (%r11,%rax)         # pn[c .. c+3]: one 32-byte store
    addq     $32, %rax                  # next 4 cells = next 32 bytes
    cmpq     %r10, %rax
    jne      .L4
```

Read it from the top: four `vmovupd` loads, four arithmetic instructions that read the other four operands straight from memory, then the arithmetic, then one store. Every one of the 8 memory accesses is 32 contiguous bytes, **all useful**, and the pointer advances by 32 bytes per iteration. Nothing else is needed, because the four neighbouring cells of `p` are already adjacent in memory.

### AoS: the same arithmetic, buried under data movement
For AoS, GCC still vectorizes, but only with 128-bit `xmm` registers (2 doubles), and each iteration advances by 96 bytes, that is 2 structs of 48 bytes. One iteration updates **2 cells in 34 instructions**:

```asm
.L4:
    vmovupd  16(%rdx,%rax), %xmm7       # 16 bytes load ...
    vmovupd  (%rdi,%rax),   %xmm15
    vmovhpd  64(%rdx,%rax), %xmm7, %xmm0   # ... high half replaced by the value 48 bytes further
    vmovupd  16(%rcx,%rax), %xmm7
    vmovhpd  48(%rdi,%rax), %xmm15, %xmm8
    vmovhpd  64(%rcx,%rax), %xmm7, %xmm6
    vmovupd  48(%r8,%rax),  %xmm7
    vaddpd   %xmm6, %xmm0, %xmm0
    vmovlpd  8(%r8,%rax),   %xmm7, %xmm6   # low half replaced, again with a value from another struct
    vmovupd  48(%rdi,%rax), %xmm7
    vmulpd   %xmm14, %xmm0, %xmm0
    vmovlpd  8(%rdi,%rax),  %xmm7, %xmm7
    vaddpd   %xmm7, %xmm6, %xmm6
    vmovupd  (%rcx,%rax),   %xmm7
    vmulpd   %xmm13, %xmm6, %xmm6
    vaddpd   %xmm6, %xmm0, %xmm0
    vmovhpd  48(%rcx,%rax), %xmm7, %xmm6
    vmovupd  (%rdx,%rax),   %xmm7
    vmovhpd  48(%rdx,%rax), %xmm7, %xmm7
    vsubpd   %xmm7, %xmm6, %xmm6
    vmovupd  (%r8,%rax),    %xmm7
    vmovhpd  48(%r8,%rax),  %xmm7, %xmm7
    vsubpd   %xmm8, %xmm7, %xmm7
    vmulpd   %xmm12, %xmm6, %xmm6
    vmulpd   %xmm11, %xmm7, %xmm7
    vaddpd   %xmm7, %xmm6, %xmm6
    vmulpd   %xmm10, %xmm6, %xmm6
    vsubpd   %xmm6, %xmm0, %xmm0
    vdivpd   %xmm9,  %xmm0, %xmm0
    vmovlpd  %xmm0, (%rsi,%rax)         # store cell 1 ...
    vmovhpd  %xmm0, 48(%rsi,%rax)       # ... and cell 2, 48 bytes later: two 8-byte stores
    addq     $96, %rax                  # 2 structs = 96 bytes
    cmpq     %r9, %rax
    jne      .L4
```

The pattern is `vmovupd` followed by `vmovhpd` or `vmovlpd`. The two values the vector needs belong to two different structs, so they are 48 bytes apart and cannot be loaded with one instruction. The compiler loads 16 bytes (of which 8 are thrown away) and overwrites one half of the register with the value from the neighbouring struct. The result is stored back the same way, as two separate 8-byte writes.

### Counting
| | SoA (`ymm`) | AoS (`xmm`) |
|---|---|---|
| Cells per iteration | 4 | 2 |
| Instructions per iteration | 21 | 34 |
| **Instructions per cell** | **5.3** | **17** |
| Arithmetic instructions | 13 | 13 |
| Loads (vector) | 4 + 4 folded into arithmetic | 8, each followed by a half-register merge |
| Stores | 1 of 32 bytes | 2 of 8 bytes |

The arithmetic is exactly the same, 13 instructions per vector. What changes is the width of the vector (4 lanes against 2) and the 18 data-movement instructions, which in SoA are 5. Per cell, AoS executes **3.2 times more instructions** for the same result, and every cache line it loads carries `u` and `v`, which the stencil never reads. This is the same 3× instruction count I saw with `perf` on the strided variant: the cause is different (there the compiler gives up on vectorizing, here it vectorizes badly), but the lesson is the same. The layout in memory decides which instructions the compiler is allowed to write.

Two last remarks. The `vdivpd` is still there: without `-ffast-math` the compiler is not allowed to replace `/ denom` with a multiplication by a precomputed reciprocal, because the result would differ in the last bit. And if you want AoS to be vectorized well, the way out is not a better compiler flag but a better layout: AoSoA, i.e. small structs of 4 or 8 values of the same field (`struct { double p[8]; }`), which keeps the fields of a block close together and still lets each load be a full vector.

<a name="pow2"></a>
## 5. The power-of-two trap

Two things in the chart are not explained by the cache sizes alone.

**The strided spikes at $N = 1024$ and $N = 2048$.** Caches are *set-associative*: an address can only be stored in a small group of slots (a "set"), chosen from some bits of the address. When the row length is a power of two, the cells `(i, j)`, `(i+1, j)`, `(i+2, j)`… are separated by exactly $8N$ bytes, and they all map to the **same set**. Only a few of them fit there, so lines are evicted long before they could be reused. At $N = 1448$ or $2896$, which are not powers of two, the same strided loop is 3–5× faster than at $2048$.

**The contiguous bumps at $N = 724$ and $N = 1448$.** These surprised me more. Here the loop order is right, and yet the cost jumps from ~2 to ~5.5 ns. The reason is the size of each array: $724^2 \cdot 8\ \text{B} \approx 4\ \text{MiB}$ and $1448^2 \cdot 8\ \text{B} \approx 16\ \text{MiB}$, just below a power of two. The allocator places large arrays on page boundaries, and printing the addresses shows that the six fields end up **exactly 4 194 304 bytes apart**. So `u[k]`, `v[k]`, `us[k]`, `vs[k]`, `p[k]`, `pn[k]` all collide in the same cache set, the same problem as before but between arrays instead of between rows.

The fix is old and simple: **padding**. Allocating the six fields in one pool, with a few extra cache lines between them, breaks the alignment:

```c
const size_t stride = cells + 8 * 5;          /* +5 cache lines per array */
double *pool = aligned_alloc_64(6 * stride * sizeof(double));
u = pool; v = pool + stride; us = pool + 2 * stride; /* ... */
```

| Grid (cells) | No padding (ns) | Padding (ns) |
|---|---|---|
| 724² | 5.60 | 2.02 |
| 1448² | 5.26 | 2.05 |

*SoA, contiguous variant; time per cell update in nanoseconds.*

A **2.6–2.8× speedup by adding 320 bytes** to each array, with identical results. This is why, for example, many structured-grid codes pad their leading dimension, and why "nice" grid sizes like 1024 or 2048 are often the worst possible choice.

<a name="realcfd"></a>
## 6. What this means for real CFD

A toy cavity solver is far from an industrial code, but the mechanisms are exactly the same, and in real solvers they are usually worse:

- **Unstructured meshes** (OpenFOAM, Fluent, SU2…) access neighbours through indirection: `phi[neighbour[face]]`. If cells that are close in space have distant indices in memory, every face loop becomes a strided, cache-missing loop like my "strided" variant. This is why mesh **renumbering** (e.g. Cuthill–McKee bandwidth reduction, `renumberMesh` in OpenFOAM) can noticeably speed up a case without changing a single number in the solution.
- **Field layout** matters in the design of the data structures: a vector field stored as $(x,y,z)$ triplets is a small AoS. It is a good choice for kernels that use all three components together and a bad one for kernels that only need one.
- **Problem size per core**: with domain decomposition, each process only works on its own sub-mesh. If that sub-mesh fits in the cache, the time per cell drops: this is one of the reasons why parallel CFD runs sometimes show *superlinear* speedup, and why the "cells per core" ratio is a hardware question as much as a numerical one.
- **Memory bandwidth is shared**: on a multi-core node the cores compete for the same RAM bandwidth. A memory-bound kernel that wastes half of every cache line, like the AoS Jacobi sweep, hits the bandwidth wall with fewer cores.
- **Benchmark on realistic sizes.** As the cache cliff shows, the performance measured on a small test mesh says little about the production run.

The broader lesson is that in scientific computing the algorithm is only half of the story. The same discretization, the same linear solver and the same compiler gave anything from 1.2 to 47 ns per cell update, a factor of about 40, depending only on how data is placed and visited. On a simulation that takes a day, that is the difference between a result tomorrow and a result next month. Thinking about memory is not premature optimization: in CFD it is part of the design.

### Reproduce it
All the files are available: [locality.c](/assets/projects/Data-Locality-CFD/locality.c), [Makefile](/assets/projects/Data-Locality-CFD/Makefile), [plot.py](/assets/projects/Data-Locality-CFD/plot.py) and the raw data [results.csv](/assets/projects/Data-Locality-CFD/results.csv).

```bash
make all
make bench > results.csv          # sweep all grid sizes and variants
python3 plot.py results.csv out.svg
perf stat -e cycles,instructions,L1-dcache-load-misses,LLC-load-misses ./bin_soa_ji 2048 3
```

Results will of course change on a different CPU (cache sizes, associativity and prefetchers differ), but the shape of the curves should not.
