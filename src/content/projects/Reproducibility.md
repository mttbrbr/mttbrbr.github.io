---
title: "Same Case, More Cores, Different Answer"
date: 2026-11-22
tech: ["C", "MPI", "CFD", "Turbulence"]
excerpt: "Run the same turbulent simulation on 1, 2, 4 and 8 processes and you get four different flows. Nothing is wrong: the last bit of a parallel sum changes, and chaos does the rest. A small 2D turbulence solver shows how fast, why it only happens in turbulent flows, and what can and cannot be compared."
---

## Technical Index
1. [**Addition is not associative**](#associative) — In floating point, order matters
2. [**The experiment**](#experiment) — 2D turbulence with a conjugate-gradient pressure solve
3. [**From the last bit to a different flow**](#divergence) — How the difference grows
4. [**What stays the same**](#statistics) — Statistics, not snapshots
5. [**What this means in practice**](#practice) — Comparing, debugging and validating

<a name="associative"></a>
## 1. Addition is not associative

In exact arithmetic $(a + b) + c = a + (b + c)$. In floating point it is not: every addition rounds the result to 53 significant bits, and different orders round differently. The difference is tiny, about one part in $10^{16}$, but it is real.

A parallel CFD solver adds numbers in an order that depends on the number of processes. A dot product in the conjugate gradient, or a residual norm, is computed as one partial sum per process, then the partial sums are added by `MPI_Allreduce`. Change the number of processes and you change the order. [mpi_dot.c](/assets/projects/Reproducibility/mpi_dot.c) computes the same dot product of a million numbers with 1 to 8 MPI processes:

| Processes | Result |
|---|---|
| 1 | 110.622246360**81555** |
| 2 | 110.622246360**80362** |
| 3 | 110.622246360**81555** |
| 4 | 110.622246360**79418** |
| 5 | 110.622246360**81794** |
| 6 | 110.622246360**77531** |
| 7 | 110.622246360**81328** |
| 8 | 110.622246360**80896** |

*The same sum, computed on 1 to 8 processes on one machine. Bold: the digits that change.*

Mathematically identical, numerically different from the 13th significant digit on (with 3 processes it happens to round exactly like with 1). In a laminar steady flow this does not matter at all. In a turbulent flow, it changes everything.

<a name="experiment"></a>
## 2. The experiment

[kolmogorov.c](/assets/projects/Reproducibility/kolmogorov.c) solves the 2D incompressible Navier–Stokes equations in vorticity–streamfunction form on a periodic box with 64 × 64 points, driven by a sinusoidal body force (*Kolmogorov flow*). At every stage of the time integration it solves a Poisson equation for the streamfunction with the conjugate gradient, exactly like a CFD code solves its pressure equation.

The same case is run five times, side by side, from the same initial condition:

- with the dot products of the conjugate gradient computed as **1, 2, 4 and 8 partial sums**, one per "process", added at the end like `MPI_Allreduce` does;
- once more with 1 process, but with **one single value** of the initial condition changed by 1 ulp (the smallest possible change of a double, about one part in $10^{16}$).

Everything else, the code, the compiler, the machine and the time step, is identical. Two flow regimes:

- **turbulent**: low viscosity, strong forcing, a chaotic flow of interacting vortices;
- **laminar**: high viscosity, weak forcing, a steady flow.

<a name="divergence"></a>
## 3. From the last bit to a different flow

The measure is the difference between each run and the 1-process run, as a relative norm of the vorticity field.

![Difference from the 1-process run vs time](/assets/projects/Reproducibility/divergence.svg)

| Time | 2 processes | 4 processes | 8 processes | 1 ulp changed | laminar, 8 processes |
|---|---|---|---|---|---|
| 5 | 2·10⁻¹³ | 5·10⁻¹³ | 4·10⁻¹³ | 4·10⁻¹³ | 1·10⁻¹⁵ |
| 20 | 2·10⁻⁸ | 4·10⁻⁸ | 2·10⁻⁸ | 4·10⁻⁸ | 2·10⁻¹⁵ |
| 40 | 5·10⁻⁵ | 3·10⁻⁴ | 3·10⁻⁴ | 2·10⁻⁴ | 2·10⁻¹⁴ |
| 50 | 0.010 | 0.033 | 0.042 | 0.028 | – |
| 60 | 0.96 | 1.30 | 1.37 | 1.22 | – |

*Relative L2 difference of the vorticity from the 1-process run.*

- **In the turbulent case the difference grows exponentially**: about a factor of 10 every 2.5 time units, from the last digit at the start to 100% at $t \approx 60$. After that, the runs are completely different flows, as different as two independent realizations can be.
- **It does not matter where the first difference comes from.** The runs with a different number of processes and the run with a single bit changed in the initial condition follow the same curve. The number of processes is just one of countless ways to perturb a chaotic system at the level of round-off.
- **In the laminar case nothing happens.** The difference stays at $10^{-15}$–$10^{-14}$ for the whole run. A stable flow damps perturbations; a chaotic one amplifies them. The parallel sum is the same in both; only the physics differs.

At the end, the three vorticity fields have nothing in common except their general appearance:

![Vorticity at t = 80 with 1 process, 8 processes and one bit changed](/assets/projects/Reproducibility/vorticity_t80.png)

*Vorticity at t = 80. Orange: counter-clockwise, blue: clockwise. Same code, same case, same machine.*

<a name="statistics"></a>
## 4. What stays the same

If the instantaneous flow is different, is the result useless? No: what a turbulent simulation can predict is its **statistics**, not its snapshots. The energy of the 1- and 8-process runs is identical until $t \approx 55$, then the two curves separate but keep fluctuating around the same level:

![Kinetic energy vs time, 1 and 8 processes](/assets/projects/Reproducibility/energy.svg)

| Run | Mean kinetic energy | Mean enstrophy |
|---|---|---|
| 1 process | 0.912 | 9.47 |
| 2 processes | 1.013 | 10.52 |
| 4 processes | 0.994 | 10.35 |
| 8 processes | 1.042 | 10.89 |
| 1 process, 1 ulp changed | 0.966 | 10.08 |
| **all five** | **0.985 ± 0.050** | **10.26 ± 0.53** |

*Time averages over $40 \le t \le 80$. Last row: mean and standard deviation of the five runs.*

The averages differ by about ±5%, and the run with one bit changed falls in the middle of the others. These differences are not caused by the number of processes: they are the **statistical uncertainty** of an average taken over 40 time units of a flow whose energy fluctuates slowly. With a single run you would see only one of these values and could easily mistake it for "the" answer. A longer averaging window would shrink the spread; changing the number of processes does not shift the mean, it only produces another sample.

<a name="practice"></a>
## 5. What this means in practice

- **Bitwise differences between runs on different core counts are normal**, for LES, DES, URANS of unsteady separated flows, and any chaotic case. Most MPI implementations and most CFD codes do not guarantee bitwise reproducibility across decompositions, and achieving it (fixed reduction orders, compensated or exact summation) costs performance.
- **Compare statistics, not snapshots.** For turbulent cases, compare time-averaged fields, spectra and integral quantities with their uncertainty, never instantaneous fields. Two runs that "look different" at the same instant may be in perfect agreement.
- **Estimate the sampling uncertainty.** The spread of the averages here is the error bar of a 40-time-unit average. In a real case, run long enough that the average stops changing, or split the run into blocks and look at the scatter of the block averages.
- **Steady and stable cases should reproduce closely.** For a converged RANS or laminar steady case, results on different core counts should agree to the solver tolerance. If they do not, that is a bug worth investigating (a decomposition-dependent boundary condition, an uninitialized value, a race condition), not chaos.
- **For debugging, fix the decomposition.** When you need bitwise-identical runs to track a problem, keep the number of processes, the decomposition, the compiler and the MPI library the same.

This closes the loop with the [article on moving data between cores](/projects/data-between-cores/): splitting a simulation across processes changes how fast it runs, and for chaotic flows it also changes which of the many possible flows you get. Both are fine, as long as you know it and compare the right things.

### Reproduce it
All the files are available: [kolmogorov.c](/assets/projects/Reproducibility/kolmogorov.c), [mpi_dot.c](/assets/projects/Reproducibility/mpi_dot.c), [plot_rp.py](/assets/projects/Reproducibility/plot_rp.py) with [svgchart.py](/assets/projects/Reproducibility/svgchart.py), and the data [divergence (turbulent)](/assets/projects/Reproducibility/divergence_nu0.002_F1.csv), [divergence (laminar)](/assets/projects/Reproducibility/divergence_nu0.05_F0.2.csv), [mpi_dot.csv](/assets/projects/Reproducibility/mpi_dot.csv).

```bash
gcc -O3 -march=native kolmogorov.c -o kolmogorov -lm
mkdir -p out
./kolmogorov out 0.002 80 1.0     # turbulent, ~8 min
./kolmogorov out 0.05 40 0.2      # laminar, ~4 min
mpicc -O2 mpi_dot.c -o mpi_dot -lm
for p in 1 2 3 4 5 6 7 8; do mpirun --oversubscribe -np $p ./mpi_dot; done
mkdir -p figures && python3 plot_rp.py out figures
```
