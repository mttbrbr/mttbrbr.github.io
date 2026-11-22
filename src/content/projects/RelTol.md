---
title: "How Accurately Should You Solve Each Linear System? relTol, Measured"
date: 2026-11-22
tech: ["C", "CFD", "OpenFOAM", "Linear Algebra"]
excerpt: "Every outer iteration of a CFD solver solves linear systems that the next iteration will throw away. Solving them too accurately wastes time, too loosely breaks convergence. A sweep of the tolerances of SIMPLE and a coupled solver shows where the sweet spot is: solving the pressure to 1e-6 was 10 times slower and converged in exactly the same number of iterations."
---

## Technical Index
1. [**Inner and outer iterations**](#inner) — Two loops, two tolerances
2. [**The experiment**](#experiment) — The lid-driven cavity, again
3. [**SIMPLE**](#simple) — Pressure and momentum tolerances
4. [**The coupled solver**](#coupled) — One tolerance, a clear minimum
5. [**What this means in practice**](#practice) — Reading `fvSolution`

<a name="inner"></a>
## 1. Inner and outer iterations

A steady CFD solver runs two nested loops. The **outer** loop handles the non-linearity and the pressure–velocity coupling: assemble the equations with the current solution, solve, update, repeat. Inside each outer iteration, the **inner** loop is an iterative linear solver (Gauss–Seidel, conjugate gradient, multigrid…) that solves each linear system to some tolerance.

The inner systems are only approximations: their coefficients come from a solution that is itself not converged yet, and the next outer iteration will rebuild them. Solving them to machine precision is pointless. Solving them too loosely can make the outer loop slower, or stop it from converging. In OpenFOAM, `fvSolution` controls this with two numbers for each equation:

- `tolerance`: an absolute residual at which the linear solver stops;
- `relTol`: a relative tolerance, the fraction by which the residual has to drop *within this solve*. With `relTol 0.1`, each solve stops as soon as the residual has dropped by 10.

This article measures what these numbers do, on the solvers from the [segregated vs coupled article](/projects/openfoam-coupled-solver/).

<a name="experiment"></a>
## 2. The experiment

[reltol.c](/assets/projects/RelTol/reltol.c) solves the steady lid-driven cavity at Re = 100 on 64 × 64 cells, with exactly the same discretization as before, until the full non-linear residual has dropped by $10^8$.

- **SIMPLE**, with under-relaxation 0.9 for momentum and 0.2 for pressure (the fastest combination found earlier). Momentum solved with Gauss–Seidel, pressure correction with Jacobi-preconditioned conjugate gradient. I swept the momentum `relTol` (0.5, 0.1, 0.01) and the pressure `relTol` (from 0.9 down to $10^{-6}$).
- **Coupled**: the whole system solved with FGMRES and a block preconditioner, sweeping its `relTol` from 0.9 down to $10^{-6}$.

Every combination converges to the same solution; only the cost changes.

<a name="simple"></a>
## 3. SIMPLE

![SIMPLE: outer iterations and wall time vs the pressure and momentum relTol](/assets/projects/RelTol/simple_reltol.svg)

| Momentum relTol | Pressure relTol | Outer iterations | Time (s) | Pressure CG iterations |
|---|---|---|---|---|
| 0.5 | 1e-6 | 871 | 15.4 | 307 574 |
| 0.5 | 0.9 | 873 | 1.47 | 17 803 |
| 0.1 | 1e-6 | 530 | 9.6 | 186 973 |
| 0.1 | 0.05 | 530 | 3.9 | 68 774 |
| 0.1 | 0.5 | 530 | **1.61** | 19 155 |
| 0.1 | 0.9 | 801 | **1.35** | 10 245 |
| 0.01 | 0.5 | 487 | 1.63 | 11 787 |
| 0.01 | 0.9 | did not converge | – | – |

*64 × 64 cells, Re = 100. Full data in [simple_reltol_N64.csv](/assets/projects/RelTol/simple_reltol_N64.csv).*

Two very different behaviours:

- **The pressure tolerance does not change the number of outer iterations.** With momentum `relTol` 0.1, SIMPLE needs 530 outer iterations whether each pressure correction is solved to 0.5 or to $10^{-6}$. The only thing a tighter pressure tolerance changes is the cost: from 1.6 to 9.6 seconds, **6 times slower for nothing** (10 times with momentum 0.5). The reason is the pressure under-relaxation: only 20% of each pressure correction is applied, and the next outer iteration recomputes the whole correction from scratch anyway. Precision in a quantity that is mostly thrown away is wasted.
- **The momentum tolerance does change the outer iterations.** Solving momentum to 0.5, 0.1 or 0.01 gives 871, 530 or 487 outer iterations. A poorly solved momentum equation gives a poor intermediate velocity, and the pressure correction has to make up for it.

The extremes have a limit, though. With the loosest pressure tolerance (0.9), the outer iterations start to grow (801 instead of 530), and combined with a tightly solved momentum it **did not converge at all**. A very inaccurate pressure correction no longer drives the velocity towards continuity. The fastest settings are close to that edge, which is a good reason not to sit exactly on it.

<a name="coupled"></a>
## 4. The coupled solver

![Coupled solver: outer iterations and wall time vs relTol](/assets/projects/RelTol/coupled_reltol.svg)

| FGMRES relTol | Outer iterations | Time (s) | Matrix–vector products |
|---|---|---|---|
| 0.9 | 147 | 2.21 | 794 |
| 0.5 | 36 | **1.62** | 619 |
| 0.1 | 16 | 2.10 | 802 |
| 0.01 | 11 | 2.29 | 876 |
| 0.001 | 11 | 3.38 | 1 288 |
| 1e-6 | 11 | 7.50 | 2 620 |

The coupled solver behaves differently, because its linear system contains the whole pressure–velocity coupling:

- Below $10^{-2}$, solving more accurately brings **nothing**: 11 outer iterations whatever the tolerance, at up to 3.3 times the cost.
- Above $10^{-2}$, the outer iterations grow, but each one becomes much cheaper. The minimum time is at `relTol 0.5`: 36 outer iterations, each solved very loosely.
- At 0.9 the linear solves barely do anything and the outer loop needs 147 iterations: too loose.

The total time varies by less than a factor of 5 over six orders of magnitude of tolerance, so the coupled solver is much less sensitive than SIMPLE. But the pattern is the same: there is no reason to solve an inner system much more accurately than the outer loop needs.

<a name="practice"></a>
## 5. What this means in practice

- **In steady SIMPLE runs, `relTol` of 0.05–0.1 for pressure and 0.1 for momentum are sensible defaults,** and the measurements support them: tighter pressure solves only cost time. Typical tutorial settings (`p` with GAMG, `relTol 0.1` or `0.05`; `U` with smoothSolver, `relTol 0.1`) are in this range.
- **Pressure is where the time goes**, so it is where loose tolerances save most. With momentum `relTol` 0.1, tightening the pressure tolerance from 0.5 to $10^{-6}$ multiplied the conjugate-gradient iterations by 10 (from 19 000 to 187 000), while the Gauss–Seidel sweeps for momentum stayed the same.
- **Do not go to the extreme.** The very loosest settings saved a little more time in some cases and broke convergence in another. A small safety margin costs little.
- **Transient runs are different.** In PISO or PIMPLE, the final corrector of each time step should leave continuity well satisfied, because there is no next outer iteration to fix it before time moves on. That is why tutorials use a `pFinal` entry with `relTol 0`: loose solves for the intermediate correctors, a full solve at the end of the step.
- **The absolute `tolerance` matters near convergence.** Once the residual falls below it, the solver does no iterations at all. If it is too large, the run stops improving; if it is far below what the outer loop can reach, it does nothing.
- **Coupled solvers are less sensitive but still have a minimum.** A loose inner tolerance with a few more outer iterations is often the fastest path.

The broader point, as in the rest of the series: accuracy is only useful where it survives. A solution computed to $10^{-6}$ and then multiplied by 0.2 and overwritten at the next iteration is time spent on nothing.

### Reproduce it
All the files are available: [reltol.c](/assets/projects/RelTol/reltol.c) (needs LAPACK), [plot_rt.py](/assets/projects/RelTol/plot_rt.py) with [svgchart.py](/assets/projects/RelTol/svgchart.py), and the data [simple_reltol_N64.csv](/assets/projects/RelTol/simple_reltol_N64.csv), [coupled_reltol_N64.csv](/assets/projects/RelTol/coupled_reltol_N64.csv).

```bash
gcc -O3 -march=native reltol.c -o reltol -llapack -lm
mkdir -p out && ./reltol out 64               # ~3 min
mkdir -p figures && python3 plot_rt.py out figures 64
```
