---
title: "Why Multigrid Always Wins: Solving the Pressure Equation"
date: 2026-11-22
tech: ["C", "CFD", "Numerical Methods", "Linear Algebra"]
excerpt: "Gauss-Seidel, conjugate gradient and multigrid on the same pressure equation, from 225 to a million unknowns. The first two get slower and slower as the mesh grows; multigrid needs 8 cycles every time. This is why GAMG is the default for pressure in OpenFOAM."
---

## Technical Index
1. [**The equation that eats the time**](#problem) — The pressure Poisson equation
2. [**Four solvers**](#solvers) — Gauss-Seidel, CG, multigrid, CG with multigrid
3. [**Results**](#results) — From 225 to a million unknowns
4. [**Why multigrid works**](#why) — Smooth errors and coarse grids
5. [**What this means in practice**](#practice) — GAMG, PCG and `fvSolution`

<a name="problem"></a>
## 1. The equation that eats the time

In every article of this series the same thing kept coming back: in an incompressible solver, most of the time goes into the **pressure equation**. Whatever the algorithm, SIMPLE, PISO or a coupled solver, at some point it has to solve a Poisson equation

$$
-\nabla^2 p = f
$$

on the whole mesh, once or several times per iteration. It is elliptic: the value in every cell depends on every other cell of the domain, so information has to travel across the whole mesh before the solution is right.

How fast that happens depends on the linear solver. In OpenFOAM's `fvSolution` you pick it for every equation: `smoothSolver`, `PCG`, `GAMG`… This article measures the difference on the simplest possible case, the Poisson equation on a square with $N \times N$ unknowns, discretized with the standard 5-point stencil.

[mg.c](/assets/projects/Multigrid/mg.c) solves it with four methods, on grids from $15^2 = 225$ to $1023^2 \approx 10^6$ unknowns. All methods stop when the residual has dropped by a factor $10^8$, and all reach the same solution. The right-hand side is chosen so that the exact solution is a smooth polynomial that does not favour any method.

<a name="solvers"></a>
## 2. Four solvers

- **Gauss–Seidel.** Sweep through the grid and replace each value with the one that satisfies its own equation, given the current neighbours. It is the simplest iterative method and the building block of OpenFOAM's `smoothSolver`.
- **Conjugate gradient (PCG).** A Krylov method: at every iteration it builds the best possible correction in a growing subspace, using one matrix–vector product and two dot products. Here it is preconditioned with the diagonal (Jacobi); OpenFOAM's `PCG` usually uses DIC, an incomplete Cholesky factorization, which is better but behaves the same way as the mesh grows.
- **Multigrid (V-cycle).** Do a few Gauss–Seidel sweeps, compute the residual, transfer it to a grid with half the points in each direction, solve the correction there (recursively, on coarser and coarser grids), bring the correction back and smooth again. Here the coarse grids are built geometrically: 2 sweeps before and 2 after on every level.
- **CG preconditioned by multigrid.** A conjugate gradient where each preconditioning step is one multigrid V-cycle. This is what `PCG` with `preconditioner GAMG` does in OpenFOAM.

<a name="results"></a>
## 3. Results

![Iterations to converge vs number of unknowns](/assets/projects/Multigrid/iterations.svg)

| Unknowns | Gauss–Seidel | PCG (Jacobi) | Multigrid | PCG + multigrid |
|---|---|---|---|---|
| 225 | 480 | 38 | 8 | 6 |
| 3 969 | 7 600 | 164 | 8 | 6 |
| 16 129 | 30 390 | 333 | 8 | 6 |
| 65 025 | 121 530 | 677 | 8 | 6 |
| 261 121 | – | 1 368 | 8 | 6 |
| 1 046 529 | – | 2 760 | 8 | 6 |

*Iterations (for multigrid: V-cycles) to reduce the residual by $10^8$. Gauss–Seidel was not run beyond 65 025 unknowns: it would have needed hours.*

The pattern is the whole story:

- **Gauss–Seidel** needs **4 times more iterations** every time the grid is refined (twice the points in each direction). Its iteration count grows like $N^2$, proportional to the number of unknowns.
- **PCG** needs **2 times more** at each refinement: its count grows like $N$, like the number of cells across the domain.
- **Multigrid** needs **8 V-cycles on every grid**, from 225 unknowns to a million. Its iteration count does not depend on the mesh.

Iterations do not cost the same, so here is the time:

![Time to solve vs number of unknowns](/assets/projects/Multigrid/time.svg)

| Unknowns | Gauss–Seidel | PCG (Jacobi) | Multigrid | PCG + multigrid |
|---|---|---|---|---|
| 16 129 | 3.3 s | 0.033 s | 0.005 s | 0.004 s |
| 65 025 | 56 s | 0.30 s | 0.025 s | 0.019 s |
| 261 121 | – | 3.1 s | 0.11 s | 0.083 s |
| 1 046 529 | – | 31 s | **0.41 s** | **0.34 s** |

*Wall time on one core (Intel i7-1165G7).*

One V-cycle costs about as much as a handful of Gauss–Seidel sweeps, because every coarse grid has a quarter of the points of the one above: the whole hierarchy adds only a third to the cost of the finest level. So the time of multigrid grows **linearly with the number of unknowns**: 4 times more unknowns, 4 times the time. The others grow much faster:

- with a million unknowns, multigrid is **77 times faster** than PCG;
- with 65 000 unknowns, it is **2 250 times faster** than Gauss–Seidel;
- every time the mesh is refined once more, the gap grows: by 2 against PCG and by 4 against Gauss–Seidel.

Using multigrid as a preconditioner for CG is a little faster still (6 iterations instead of 8): the Krylov method cleans up what a single V-cycle leaves behind.

The convergence histories show the different behaviours. Gauss–Seidel's residual barely moves for thousands of iterations and then falls; PCG has a long plateau and then drops; multigrid reduces the residual by about a factor of 10 at every cycle, from the first one:

![Convergence histories on 16 129 unknowns](/assets/projects/Multigrid/history.svg)

<a name="why"></a>
## 4. Why multigrid works

The error of an approximate solution can be seen as a sum of waves of different wavelengths, from the size of a cell to the size of the domain. Gauss–Seidel treats them very differently:

- **Short waves** (errors that oscillate from cell to cell) are damped very quickly: each cell is updated with its neighbours, and an oscillation between neighbours is exactly what the update averages out. After a couple of sweeps, the error is *smooth*.
- **Long waves** (errors that vary slowly across the domain) are damped extremely slowly: each sweep moves information only one cell, so a correction that must cross $N$ cells needs a number of sweeps that grows like $N^2$. This is the long plateau in the history above.

The idea of multigrid is that a wave that is long on a fine grid is short on a coarse one. After a few sweeps the remaining error is smooth, so it can be represented on a grid with half the points, where it oscillates twice as fast and Gauss–Seidel is effective again. Repeat on coarser and coarser grids, and every wavelength meets a grid on which it is "short". Each V-cycle reduces all the components of the error by roughly the same factor, here about 10, and that factor does not depend on the size of the grid. Hence a constant number of cycles.

The conjugate gradient sits in between: it does not move information one cell at a time, but its convergence still depends on the ratio between the longest and the shortest wave the grid can represent (the condition number of the matrix), which grows like $N^2$. The number of iterations grows like its square root, $N$.

<a name="practice"></a>
## 5. What this means in practice

OpenFOAM's `GAMG` is an *algebraic* multigrid: instead of a geometric hierarchy of grids, it builds the coarse levels by agglomerating neighbouring cells of the unstructured mesh, using only the matrix. On real 3D meshes it is not as perfect as the geometric version here, and the iteration count grows slowly with the mesh, but the principle and the advantage are the same.

- **Use `GAMG` for the pressure equation**, either as a solver or as a preconditioner of `PCG`. It is the default in most tutorials for a reason: on large meshes it is often an order of magnitude faster than `PCG` with DIC.
- **The advantage grows with the mesh.** On a few thousand cells all methods are fast enough. On millions of cells, a solver whose iterations grow with the mesh becomes the bottleneck of the whole simulation.
- **Momentum and turbulence equations are different.** They are dominated by convection and by the diagonal, and with relaxation or small time steps their matrix is well conditioned: a few sweeps of `smoothSolver` with Gauss–Seidel are usually enough. Multigrid pays off on elliptic equations like pressure.
- **GAMG has its own knobs.** `smoother` (`GaussSeidel`, `DIC`, `DICGaussSeidel`), `nPreSweeps` and `nPostSweeps`, `nCellsInCoarsestLevel`, `agglomerator`. They correspond directly to the pieces of the V-cycle above: the smoother, the sweeps on every level and the size of the coarsest grid.
- **In parallel, the coarse levels are the weak point.** The coarsest grids have very few cells per process, so the communication of the halo exchange ([article on moving data between cores](/projects/data-between-cores/)) dominates there. This is why GAMG offers options to gather the coarsest levels onto fewer processes.

The same theme as the rest of the series: the equations do not change, the answer does not change, and the cost changes by orders of magnitude. Here not because of how data sits in memory, but because of how information travels through the mesh.

### Reproduce it
All the files are available: [mg.c](/assets/projects/Multigrid/mg.c), [plot_mg.py](/assets/projects/Multigrid/plot_mg.py) with [svgchart.py](/assets/projects/Multigrid/svgchart.py), and the data [scaling.csv](/assets/projects/Multigrid/scaling.csv), [history.csv](/assets/projects/Multigrid/history.csv).

```bash
gcc -O3 -march=native mg.c -o mg -lm
mkdir -p out && ./mg out                  # ~1.5 min, most of it Gauss-Seidel on 255 x 255
mkdir -p figures && python3 plot_mg.py out figures
```
