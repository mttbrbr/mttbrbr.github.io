---
title: "Compiler Flags, Precision and Cache: What Really Speeds Up a CFD Solver"
date: 2026-11-22
tech: ["C", "OpenFOAM", "Linux", "Hardware"]
excerpt: "Native compilation, -ffast-math, single and mixed precision, big caches: which of them actually make a CFD solver faster? A small C model of OpenFOAM's linear solver answers with numbers, and explains why -march=native gave me 2% and mixed precision 27%."
---

## Technical Index
1. [**The starting point**](#start) — Two numbers from a real OpenFOAM run
2. [**A tiny OpenFOAM**](#model) — The LDU matrix and the conjugate gradient
3. [**Compiler flags**](#flags) — Why native compilation does so little
4. [**The price of -ffast-math**](#fastmath) — Two things it is allowed to break
5. [**Single, double and mixed precision**](#precision) — Speed per iteration
6. [**Precision and convergence**](#convergence) — When single precision lies
7. [**Cache, CCDs and asymmetric processors**](#cache) — What a bigger L3 changes
8. [**What this means in practice**](#practice) — Where the real gains are

<a name="start"></a>
## 1. The starting point

On my workstation (Ryzen 9 7900X3D, 12 cores on two CCDs, one of them with 64 MB of extra stacked L3 cache, DDR5-6000) I ran the [HPC motorbike case](https://develop.openfoam.com/committees/hpc/-/tree/develop/incompressible/simpleFoam/HPC_motorbike) of OpenFOAM v2312, 17 million cells, `simpleFoam`, 500 iterations in parallel, with three different builds:

| OpenFOAM build | Time (s) | Speedup |
|---|---|---|
| Standard package (generic x86-64) | 8 826 | 1.00× |
| Compiled with `-march=znver4` | 8 656 | 1.02× |
| Mixed precision (SPDP) | 6 399 | 1.38× |

*Wall time for 500 iterations, same case, same parallel setup.*

Recompiling the whole of OpenFOAM for my exact processor gained **2%**. Switching to mixed precision gained **38%** (27% less time). The common advice usually says the opposite: compile natively, enable AVX, and precision is a matter of accuracy, not speed.

Running every combination on the real case takes days. Instead, as in the previous articles, I wrote a small C program that reproduces the part of OpenFOAM where the time goes, and measured it with every option. All the measurements below come from that program on a laptop (Intel Core i7-1165G7, 4 cores, L2 1.25 MB per core, L3 12 MB, GCC 16.2 and Clang, CachyOS Linux), median of 3 runs.

<a name="model"></a>
## 2. A tiny OpenFOAM

In an incompressible OpenFOAM run most of the time goes into the linear solvers, especially the pressure equation. [ldu_cg.c](/assets/projects/OpenFOAM-socket-Optimization/ldu_cg.c) reproduces exactly that, in about 300 lines:

- **Mesh:** a cube of $n \times n \times n$ hexahedra, stored as if it were unstructured: a list of internal faces, each with an *owner* and a *neighbour* cell.
- **Matrix:** OpenFOAM's **LDU** format: one diagonal coefficient per cell and one off-diagonal coefficient per face. It is a discretized Laplacian, the same kind of matrix as a pressure equation, with non-uniform coefficients.
- **Solver:** preconditioned conjugate gradient (PCG) with a diagonal preconditioner. The matrix-vector product is OpenFOAM's `Amul`, a loop over faces:

```c
for (c = 0; c < nCells; c++) y[c] = diag[c] * x[c];
for (f = 0; f < nFaces; f++) {
    y[own[f]] += upper[f] * x[nei[f]];
    y[nei[f]] += upper[f] * x[own[f]];
}
```

- **Field operations:** a second kernel represents the rest of an iteration: interpolate a vector field to the faces, compute the face flux, accumulate its divergence in the cells and update the field. It is a simplified version of `fvc::interpolate` and `fvc::div`.

The metric is the time to process **one cell in one iteration**, in nanoseconds.

<a name="flags"></a>
## 3. Compiler flags

The flags I compared:

- `-O0`: no optimization, used for debugging.
- `-O2`: the standard level for distributed binaries, a compromise between speed and binary size. With `-march=x86-64 -mtune=generic` it produces code that runs on any 64-bit x86 processor, like a distribution package.
- `-O3`: more aggressive optimization (more unrolling and vectorization); it is the default of OpenFOAM's own build.
- `-march=native`: use every instruction the processor has (AVX2, AVX-512 on this CPU), at the price of a binary that only runs on similar machines. `-mprefer-vector-width=512` asks GCC to actually use the 512-bit registers.
- `-ffast-math`: allow the compiler to ignore some IEEE floating-point rules (section 4).

![Time per cell per iteration of the solver kernel for each set of flags](/assets/projects/OpenFOAM-socket-Optimization/flags.svg)

| Flags | 24³ cells, in cache (ns) | 96³ cells, in RAM (ns) |
|---|---|---|
| gcc -O0 | 47.0 | 48.5 |
| gcc -O2, generic x86-64 | 14.0 | 19.1 |
| gcc -O3 | 14.0 | 18.9 |
| gcc -O3 -march=native | 13.1 | 18.8 |
| … + 512-bit vectors | 13.3 | 19.7 |
| … + -ffast-math | **9.4** | **18.0** |
| clang -O3 | 14.5 | 18.9 |
| clang -O3 -march=native | 13.6 | 19.8 |
| … + -ffast-math | 9.8 | 18.6 |

*Solver kernel in double precision: time per cell per CG iteration in nanoseconds, lower is better. 24³ cells take 1.3 MB (they fit in cache), 96³ cells take 84 MB (they live in RAM).*

What the numbers say:

- **Optimizing at all is what matters**: `-O0` is 3.4 times slower than everything else.
- **`-march=native` changes almost nothing**: 0–7%, and 512-bit vectors are even slightly slower. This is my 2% on the real case. The reason is in the face loop above: every face writes to two cells chosen by `own[f]` and `nei[f]`. The compiler cannot know that two faces never write the same cell at the same time, so it cannot process several faces in one vector instruction. The most important loop of the solver simply **cannot be vectorized**, whatever the instruction set.
- **`-ffast-math` gives 28% in cache, 4% in RAM.** The gain comes from the dot products of the conjugate gradient. Under strict IEEE rules a sum must be done in the written order, one term after the other; with `-ffast-math` the compiler may reorder it and add four or eight terms at a time in vector registers. When the data comes from RAM, waiting for memory hides most of this gain.
- **GCC and Clang are equivalent** here, within a few percent.

<a name="fastmath"></a>
## 4. The price of -ffast-math

`-ffast-math` is a bundle of options that let the compiler assume that floating-point arithmetic behaves like real-number algebra: operations can be reordered, $a/b$ can become $a \cdot (1/b)$, and NaN and infinity never occur. Most of the time the result changes only in the last digits. But [fastmath_traps.c](/assets/projects/OpenFOAM-socket-Optimization/fastmath_traps.c) shows two things it is allowed to break, with the same source compiled twice:

| Test | gcc -O3 | gcc -O3 -ffast-math |
|---|---|---|
| Kahan (compensated) sum, error | 0 | 2.2e-3, same as the plain sum |
| Residual computed as 0/0 | NaN | **1.0** |
| Divergence check `isnan(r) \|\| isinf(r)` | detects it | **does not** |

*Clang with `-ffast-math` behaves the same way.*

- **Compensated summation disappears.** The Kahan algorithm keeps a correction term that is algebraically zero; under real-number algebra the compiler is entitled to remove it, and it does.
- **NaN detection disappears.** With `-ffast-math` the compiler assumes NaN never happens: it simplifies $x/x$ to 1 and turns `isnan()` into "false". A solver that checks its residual this way to stop a diverging run would keep going with garbage.

None of this means `-ffast-math` is forbidden, but it has to be tested on the real code, not assumed safe.

<a name="precision"></a>
## 5. Single, double and mixed precision

OpenFOAM can be compiled in three precisions (`WM_PRECISION_OPTION`):

- **DP:** everything in double precision (64-bit, about 16 significant digits). The default.
- **SP:** everything in single precision (32-bit, about 7 digits).
- **SPDP:** fields and matrix coefficients stored in single precision, but the linear solvers work in double (OpenFOAM's `solveScalar`).

The model reproduces the three options. Halving the size of a number halves the bytes moved from memory, so the first question is how much faster each iteration becomes.

![Time per cell per iteration vs mesh size, for the three precisions](/assets/projects/OpenFOAM-socket-Optimization/precision_size.svg)

On a single core the gain is moderate: 1.2–1.4 times for SP. The surprise is SPDP: in the **solver** it gains nothing, because the solver's work vectors stay in double; in the **field operations** it gains exactly as much as SP, because there everything is in single precision.

A single core is not how CFD runs, though. The real case ran in parallel, with many cores sharing the same memory bus. To reproduce that, I ran 4 copies of the program at once, one per core, each with its own 96³ mesh:

![Speedup of single and mixed precision over double, with 1 and 4 processes](/assets/projects/OpenFOAM-socket-Optimization/precision_speedup.svg)

| Kernel | Processes | DP (ns) | SP (ns) | SPDP (ns) | SP speedup | SPDP speedup |
|---|---|---|---|---|---|---|
| solver | 1 | 19.1 | 13.6 | 18.9 | 1.41× | 1.01× |
| solver | 4 | 41.8 | 23.2 | 41.3 | 1.80× | 1.01× |
| fields | 1 | 19.0 | 13.7 | 13.7 | 1.39× | 1.38× |
| fields | 4 | 55.6 | 32.5 | 34.7 | 1.71× | 1.60× |

*Time per cell per iteration with 96³ cells per process, in nanoseconds.*

With 4 processes, double precision becomes more than **twice as slow per cell**: the cores are now waiting for the shared memory bus, and that is when moving half the bytes pays off most, up to 1.8 times faster.

This also explains the 38% of the real case. With SPDP, OpenFOAM's field operations, matrix assembly, interpolations and gradients all run in single precision and get the full gain, while the linear solvers keep their double-precision vectors and gain little. The total lands in between, depending on how the time is split. The model does not predict the exact number, but it shows where it comes from.

<a name="convergence"></a>
## 6. Precision and convergence

Speed is only half of the question. The other half is whether the solver still converges to the right answer. The model solves the same system to a relative residual of $10^{-12}$ and, at every step, also computes the **true residual**, with the matrix rebuilt in double precision.

![Convergence of the conjugate gradient in double, single and mixed precision](/assets/projects/OpenFOAM-socket-Optimization/convergence.svg)

| Cells | Precision | CG iterations | True residual reached | Error vs DP solution |
|---|---|---|---|---|
| 24³ | DP | 139 | 8.5e-13 | – |
| 24³ | SPDP | 139 | 2.0e-6 | 9.3e-8 |
| 24³ | SP | 152 | 1.9e-5 | 2.6e-7 |
| 96³ | DP | 497 | 1.7e-12 | – |
| 96³ | SPDP | 496 | 3.4e-5 | 1.4e-7 |
| 96³ | SP | **762** | **5.9e-4** | 6.3e-7 |

*True residual: relative residual of the final solution with the double-precision matrix. Error: relative L2 difference from the double-precision solution.*

- **Single precision lies.** The residual that the SP conjugate gradient computes for itself (dashed line) keeps falling and reaches $10^{-12}$, so the solver reports convergence. The true residual stopped at $6 \cdot 10^{-4}$ hundreds of iterations earlier. On top of that, SP needed **53% more iterations** on the 96³ mesh, which eats most of its speed advantage per iteration.
- **The floor rises with the mesh.** Finer meshes make the pressure matrix worse conditioned, and the limited precision shows more. A residual floor that does not matter on a coarse mesh can stall a fine one.
- **Mixed precision converges like double.** SPDP takes the same number of iterations as DP. Its residual floor (from $2 \cdot 10^{-6}$ on 24³ to $3 \cdot 10^{-5}$ on 96³) comes from storing the matrix in single precision, which perturbs the problem slightly. It is below the relative tolerances normally used for the pressure equation inside each outer iteration (`relTol` 0.01–0.1), but it would stop a solver asked to converge much further.
- **The solution error is tiny.** In all cases the solution differs from the double-precision one by less than $10^{-6}$, far below the discretization error of any real mesh. In this well-conditioned test, the risk of single precision is not a wrong answer but a solver that does not really converge. In cases that are harder to converge, or that compute small differences between large numbers (like the pressure at low Mach number, see the [article on density-based solvers](/projects/pressure-vs-density-based/)), double precision is what keeps the solution meaningful.

<a name="cache"></a>
## 7. Cache, CCDs and asymmetric processors

The 7900X3D is an asymmetric processor: its two CCDs have 6 cores each, but one has 96 MB of L3 cache (32 MB plus 64 MB of stacked "3D V-Cache") and the other 32 MB with higher clock speeds. Intel's hybrid processors are asymmetric in a different way, with performance and efficiency cores. I cannot reproduce this hardware on a laptop with 4 identical cores, so this section is not a measurement of the 7900X3D, but the model shows the two concepts behind it.

**What a bigger L3 changes.** The plot of section 5 is also a plot against the size of the problem. The solver needs about 95 bytes per cell in double precision: 32³ cells (3 MB) fit in the 12 MB L3 of the laptop, 64³ (25 MB) do not. Between the two, the time per cell grows by about 20%, and at 96³ it is 30–50% higher than in cache. A larger L3 moves this transition to larger meshes: that is all the V-Cache does. Two consequences:

- The benefit depends on the **cells per CCD**, not on the size of the case. A real OpenFOAM solver uses much more memory than this model, of the order of 1–2 KB per cell for an incompressible case, so 96 MB of L3 holds tens of thousands of cells. For the 17-million-cell motorbike on one workstation, every CCD works on millions of cells and runs mostly from RAM, whatever its cache. The extra cache helps when the mesh is split across many CCDs or nodes, so that each piece gets small.
- Cache-related effects can be large and surprising. At exactly 128³ cells the double-precision solver was **about 75% slower** than at 127³ or 129³: every vector then occupies exactly 16 MiB, the vectors land on the same cache sets and evict each other. This is the power-of-two trap described in the [data locality article](/projects/data-locality-cfd/).

**Why asymmetry is a problem for MPI.** In a parallel run every process waits for the slowest one at every iteration (see the [article on moving data between cores](/projects/data-between-cores/)). If some processes run on the V-Cache CCD and others on the high-frequency CCD, or some on performance cores and others on efficiency cores, they run at different speeds even with a perfectly balanced mesh, and the whole run goes at the speed of the slowest. Two practical options: pin all the processes to identical cores (for example `mpirun --bind-to core` with an explicit CPU list, or `taskset`), or decompose the mesh with weights proportional to the speed of each core.

<a name="practice"></a>
## 8. What this means in practice

- **A generic optimized build is already almost as fast as a native one.** The face loops at the heart of a finite-volume solver cannot be vectorized, so `-march=native` and AVX-512 bring a few percent. Compile natively if it is easy, but do not expect miracles.
- **Never run unoptimized code.** `-O0` (or a debug build) is 3–4 times slower.
- **`-ffast-math` gives a real gain in cache-resident code, but it can silently break NaN checks and compensated sums.** If you use it, test the whole code, including what happens when a run diverges.
- **Mixed precision (SPDP) is the biggest lever on a memory-bound multi-core run:** fields in single precision move half the bytes, while the linear solvers keep double-precision vectors and converge like DP. On the real case it gained 38%; in the model, the field operations ran 1.6–1.7 times faster on 4 cores.
- **Pure single precision is risky:** more iterations, a residual floor that rises with the mesh, and a solver that reports convergence it has not reached.
- **Cache helps when the cells per core are few.** A big L3 changes the picture only if each piece of the decomposed mesh is small enough; on asymmetric processors, pin the processes to identical cores or balance the decomposition accordingly.

The common thread of this series: in a memory-bound code, the way to go faster is rarely a better instruction and almost always fewer bytes moved.

### Reproduce it
All the files are available: [ldu_cg.c](/assets/projects/OpenFOAM-socket-Optimization/ldu_cg.c), [fastmath_traps.c](/assets/projects/OpenFOAM-socket-Optimization/fastmath_traps.c), [bench.sh](/assets/projects/OpenFOAM-socket-Optimization/bench.sh), [plot_cp.py](/assets/projects/OpenFOAM-socket-Optimization/plot_cp.py) and the data [results.csv](/assets/projects/OpenFOAM-socket-Optimization/results.csv), [solve_dp_96.csv](/assets/projects/OpenFOAM-socket-Optimization/solve_dp_96.csv), [solve_sp_96.csv](/assets/projects/OpenFOAM-socket-Optimization/solve_sp_96.csv), [solve_spdp_96.csv](/assets/projects/OpenFOAM-socket-Optimization/solve_spdp_96.csv).

```bash
./bench.sh > results.csv                 # builds all variants and runs every timing test, ~8 min
./bin/dp 96 solve > solve_dp_96.csv      # convergence test (also sp, spdp)
gcc -O3 fastmath_traps.c -o traps_safe -lm && ./traps_safe
gcc -O3 -ffast-math fastmath_traps.c -o traps_fast -lm && ./traps_fast
mkdir -p figures && python3 plot_cp.py . figures
```
