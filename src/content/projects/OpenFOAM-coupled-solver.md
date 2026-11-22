---
title: "Segregated vs Coupled: Solving Pressure and Velocity Together"
date: 2026-11-22
tech: ["C", "CFD", "Numerical Methods", "Linear Algebra"]
excerpt: "SIMPLE solves velocity and pressure one after the other; a coupled solver solves them at once. Same equations, same mesh, same answer: measured on a lid-driven cavity, the difference is 12 outer iterations against thousands, paid with memory and a harder linear system."
---

## Technical Index
1. [**The pressure–velocity problem**](#problem) — Why incompressible flow has no equation for pressure
2. [**Two strategies**](#strategies) — Segregated (SIMPLE) and coupled
3. [**The experiment**](#experiment) — A lid-driven cavity, solved both ways
4. [**Outer iterations**](#iterations) — 12 against thousands
5. [**Under-relaxation**](#relaxation) — The knobs SIMPLE needs and coupled does not
6. [**Cost with a direct solver**](#direct) — A first, misleading answer
7. [**Cost with iterative solvers**](#iterative) — The real trade-off
8. [**What this means in practice**](#practice) — When coupling pays off

Every pressure-based CFD code has to answer the same question: how do you solve for velocity and pressure when the two are tied together but only one of them has its own equation? Most codes answer with a **segregated** algorithm (SIMPLE, PISO, PIMPLE): solve one, then the other, and iterate. Some answer with a **coupled** solver: put everything in one linear system and solve it at once. Fluent and STAR-CCM+ offer both; OpenFOAM's standard solvers are segregated.

In this article I compare the two on the same problem, with the same discretization, and measure what each one really costs.

<a name="problem"></a>
## 1. The pressure–velocity problem

The steady incompressible Navier–Stokes equations are

$$
\nabla\cdot(\mathbf{u}\otimes\mathbf{u}) - \nabla\cdot(\nu\nabla\mathbf{u}) + \nabla p = 0, \qquad \nabla\cdot\mathbf{u} = 0
$$

The momentum equation gives the velocity, if you know the pressure. But there is no equation that gives the pressure: continuity does not contain $p$ at all. It is a *constraint* on the velocity, and the pressure is whatever it must be for the velocity to satisfy it (as discussed in the [article on low-Mach flows](/projects/pressure-vs-density-based/), in the incompressible limit pressure is a Lagrange multiplier).

After discretization, with velocity unknowns $\mathbf{u}$ and pressure unknowns $\mathbf{p}$, the linearized system has the block structure

$$
\begin{bmatrix} \mathbf{A} & \mathbf{G} \\ \mathbf{D} & \mathbf{0} \end{bmatrix}
\begin{bmatrix} \mathbf{u} \\ \mathbf{p} \end{bmatrix} =
\begin{bmatrix} \mathbf{b} \\ \mathbf{0} \end{bmatrix}
$$

where $\mathbf{A}$ is the momentum (convection–diffusion) matrix, $\mathbf{G}$ the discrete pressure gradient and $\mathbf{D}$ the discrete divergence. The zero block on the diagonal is the mathematical form of "continuity has no pressure in it". A matrix like this is called a **saddle-point** system: it is indefinite, and the standard iterative solvers and preconditioners used for momentum or pressure alone do not work on it directly.

There is one more difficulty: the momentum equation is non-linear, because $\mathbf{A}$ depends on $\mathbf{u}$ through the convective fluxes. Every method therefore needs an **outer iteration**, where the matrix is rebuilt from the latest velocity and a linear problem is solved.

<a name="strategies"></a>
## 2. Two strategies

### Segregated: SIMPLE
SIMPLE (Patankar and Spalding, 1972) never forms the block system. At every outer iteration it:

1. solves the momentum equations with the pressure of the previous iteration, obtaining a velocity $\mathbf{u}^*$ that does not satisfy continuity;
2. builds a **pressure-correction** equation by approximating $\mathbf{A}$ with its diagonal: a velocity change responds to a pressure change as $\mathbf{u}' \approx -\text{diag}(\mathbf{A})^{-1}\mathbf{G}\,\mathbf{p}'$, and imposing continuity on $\mathbf{u}^* + \mathbf{u}'$ gives a Poisson-like equation for $\mathbf{p}'$;
3. corrects velocity and pressure.

Each step is cheap, because it only involves scalar systems the size of one field. But the diagonal approximation neglects the influence of neighbouring cells, so the correction overshoots. SIMPLE only converges with **under-relaxation**: momentum is relaxed with a factor $\alpha_u$, the pressure correction with $\alpha_p$, both below 1. And since velocity is always computed with an out-of-date pressure, information about the coupling crawls across the domain one outer iteration at a time.

### Coupled
A coupled solver assembles the full block system above at every outer iteration and solves it for $\mathbf{u}$ and $\mathbf{p}$ together. The pressure–velocity coupling is exact within each linear solve, so nothing lags and no pressure relaxation is needed. The price is the linear system: it contains all the unknowns at once and it is a saddle-point matrix.

The standard way to solve it iteratively is a Krylov method (GMRES) with a **block preconditioner** that splits velocity and pressure. Eliminating the velocity from the block system gives the pressure **Schur complement**

$$
\mathbf{S} = -\mathbf{D}\,\mathbf{A}^{-1}\mathbf{G} \;\approx\; -\mathbf{D}\,\text{diag}(\mathbf{A})^{-1}\mathbf{G}
$$

The exact $\mathbf{A}^{-1}$ is a dense matrix and cannot be formed, so it is approximated with the inverse of the diagonal: exactly the same approximation that SIMPLE uses to build its pressure-correction equation (in PETSc this is the field-split preconditioner with the `selfp` Schur option). This is the key point of the whole comparison: **SIMPLE and the coupled solver use the same approximation**. SIMPLE applies it as an outer iteration, accepting the error and fixing it with many more iterations; the coupled solver applies it inside a Krylov method, which corrects the error at every linear solve.

<a name="experiment"></a>
## 3. The experiment

[cavity.c](/assets/projects/OpenFOAM-coupled-solver/cavity.c) solves the classic **lid-driven cavity**: a square box of fluid, closed on all sides, with the top wall sliding at constant speed. The whole program is about 600 lines of C.

- **Discretization:** finite volumes on a staggered (MAC) grid, with velocity on the faces and pressure at the cell centres, so no Rhie–Chow interpolation is needed. First-order upwind convection, central diffusion. Identical for both methods.
- **Non-linearity:** treated the same way in both, with Picard iteration (convective fluxes from the previous outer iteration).
- **Convergence:** both methods stop when the **same** full non-linear residual (momentum plus continuity) has dropped by a factor $10^8$.
- **Grids:** 16², 32², 64² and 128² cells; Reynolds numbers 100, 400 and 1000.

The two methods converge to the same discrete solution. The minimum horizontal velocity on the vertical centreline agrees to five digits in every case:

| Cells | $u_{min}$, SIMPLE | $u_{min}$, coupled |
|---|---|---|
| 32² | −0.18867 | −0.18867 |
| 64² | −0.20047 | −0.20047 |
| 128² | −0.20686 | −0.20686 |

*Re = 100. Reference value (Ghia et al., 1982, 128² cells, higher-order scheme): −0.2109.*

The gap from the reference shrinks as the grid is refined and comes from the numerical diffusion of first-order upwind (see the [article on numerical diffusion](/projects/numerical-diffusion/)), not from the solvers. Every difference below is about *how* the coupling is solved, not *what* is solved.

<a name="iterations"></a>
## 4. Outer iterations

The first, cleanest comparison is the number of outer (non-linear) iterations.

![Outer iterations vs grid size for SIMPLE and the coupled solver](/assets/projects/OpenFOAM-coupled-solver/outer_iterations.svg)

| Cells | SIMPLE 0.7 / 0.3 | SIMPLE 0.9 / 0.2 | Coupled |
|---|---|---|---|
| 32² | 557 | 156 | 11 |
| 64² | 2 000 | 530 | 11 |
| 128² | 7 279 | 1 907 | 12 |

*Outer iterations to reduce the residual by $10^8$, Re = 100. SIMPLE 0.7 / 0.3: common default under-relaxation ($\alpha_u$ / $\alpha_p$). 0.9 / 0.2: the fastest combination found on 32² cells.*

- The coupled solver needs **11–12 outer iterations on every grid**. Its convergence is limited only by the non-linearity, which does not depend on the mesh.
- SIMPLE needs hundreds to thousands, and about **3.5 times more each time the grid is refined**. A finer grid means more cells for the pressure information to cross, and SIMPLE moves it one lagged iteration at a time. Even with the best under-relaxation, at 128² it is still 160 times more outer iterations.

The convergence histories show the two behaviours: the coupled residual falls steeply from the start, while SIMPLE first crawls for hundreds of iterations.

![Convergence history on 64 x 64 cells](/assets/projects/OpenFOAM-coupled-solver/history.svg)

<a name="relaxation"></a>
## 5. Under-relaxation

SIMPLE has two tuning parameters, and they matter a lot. I ran it on the 32² grid with every combination of $\alpha_u$ and $\alpha_p$ in a small table:

![SIMPLE outer iterations as a function of the under-relaxation factors](/assets/projects/OpenFOAM-coupled-solver/relaxation_map.svg)

- The number of iterations ranges from **142 to 1 178**, a factor of 8, just by changing two numbers.
- About half of the combinations **diverge**.
- The fastest settings sit right next to the divergence boundary, along the classic rule of thumb $\alpha_p \approx 1 - \alpha_u$. Finding them means trial and error, and the best values change with the mesh and the flow.

The coupled solver has no such parameter. It converged without any under-relaxation at all three Reynolds numbers tested:

| Re | SIMPLE 0.7 / 0.3 | Coupled |
|---|---|---|
| 100 | 511 | 11 |
| 400 | 467 | 18 |
| 1000 | 423 | 22 |

*Outer iterations on 32² cells, direct linear solvers.*

The coupled iterations grow with the Reynolds number because convection becomes stronger and the Picard linearization (fluxes from the previous iteration) becomes a worse approximation: that is a non-linearity effect, not a pressure–velocity coupling effect. A Newton linearization would reduce it further.

<a name="direct"></a>
## 6. Cost with a direct solver

Fewer iterations do not automatically mean less time: each coupled iteration solves a much larger and harder system. The first version of the program solved every linear system exactly, with a banded LU factorization (LAPACK), for both methods:

| Cells | Method | Outer iterations | Time (s) | Memory (MB) |
|---|---|---|---|---|
| 16² | SIMPLE 0.7 / 0.3 | 151 | 0.05 | 0.1 |
| 16² | Coupled | 12 | 0.01 | 0.8 |
| 32² | SIMPLE 0.7 / 0.3 | 511 | 1.03 | 0.8 |
| 32² | Coupled | 11 | 0.41 | 6.9 |
| 64² | SIMPLE 0.7 / 0.3 | 1 829 | 67.3 | 6.3 |
| 64² | Coupled | 11 | 5.96 | 55.8 |

*Re = 100, one core (i7-1165G7). Memory: storage of the LU factors of the largest system, in MB.*

The coupled solver looks 11 times faster, but this comparison is unfair to SIMPLE. Its strength is that each iteration can be *cheap and approximate*: real codes never solve its momentum and pressure equations exactly at every outer iteration, because the next iteration changes everything anyway. With an exact direct solve, SIMPLE pays the full price thousands of times. The memory column, on the other hand, is already telling: the coupled LU takes **9 times more memory**, and it grows much faster with the grid, because the bandwidth of the coupled matrix is three times larger.

<a name="iterative"></a>
## 7. Cost with iterative solvers

So I wrote the version that matches what production codes do:

- **SIMPLE:** momentum solved with Gauss–Seidel until its residual drops by 10× (relative tolerance 0.1); pressure correction with Jacobi-preconditioned conjugate gradient to a relative tolerance of 0.05. These are typical values in an OpenFOAM `fvSolution`.
- **Coupled:** the full system solved with flexible GMRES to a relative tolerance of $10^{-2}$, with the block preconditioner of section 2: two Gauss–Seidel sweeps for the velocity, the Schur complement approximated with $\text{diag}(\mathbf{A})^{-1}$ and solved with conjugate gradient. Inexact linear solves, the same outer iteration.

![Wall time and memory vs grid size](/assets/projects/OpenFOAM-coupled-solver/time_memory.svg)

| Cells | Method | Outer iterations | Time (s) | Memory (MB) |
|---|---|---|---|---|
| 32² | SIMPLE 0.7 / 0.3 | 557 | 0.53 | 0.2 |
| 32² | SIMPLE 0.9 / 0.2 | 156 | **0.14** | 0.2 |
| 32² | Coupled | 11 | 0.16 | 1.7 |
| 64² | SIMPLE 0.7 / 0.3 | 2 000 | 14.2 | 0.9 |
| 64² | SIMPLE 0.9 / 0.2 | 530 | 4.0 | 0.9 |
| 64² | Coupled | 11 | **2.4** | 6.8 |
| 128² | SIMPLE 0.7 / 0.3 | 7 279 | 410 | 3.4 |
| 128² | SIMPLE 0.9 / 0.2 | 1 907 | 106 | 3.4 |
| 128² | Coupled | 12 | **37.5** | 26.9 |

*Re = 100, one core. Memory: storage of the linear solvers (matrix coefficients, work vectors and, for the coupled solver, the Krylov basis), in MB.*

Now the picture is a real trade-off:

- **On a small grid a well-tuned SIMPLE is as fast or faster**: 0.14 s against 0.16 s on 32². Each of its iterations costs almost nothing.
- **As the grid grows, the coupled solver wins**: 1.7 times faster than the best SIMPLE on 64², **2.8 times faster** on 128², and 11 times faster than SIMPLE with default under-relaxation. The break-even point moves in favour of coupling as the mesh grows, because SIMPLE's iteration count grows with the mesh and the coupled one does not.
- **The memory stays about 8 times larger** for the coupled solver, now mostly because FGMRES keeps its Krylov basis: up to 60 vectors, each as long as all the unknowns together.

Where does the time go? Almost entirely into the pressure. On 128² cells, SIMPLE with default factors performed **1.94 million** conjugate-gradient iterations on the pressure correction in total, the coupled solver **164 000** on its Schur complement. The coupled solver solves the pressure a few times, approximately, inside a method that corrects the approximation; SIMPLE solves it thousands of times and corrects the approximation with more outer iterations.

Two caveats. The pressure solves here use a simple Jacobi-preconditioned conjugate gradient. With algebraic multigrid (GAMG in OpenFOAM, BoomerAMG or GAMG in PETSc), both methods would become much faster and the ratio between them could change. And this is a 2D laminar problem: in 3D, with turbulence and distorted meshes, the coupled system is larger and harder to precondition, while turbulence equations are usually still solved separately in both approaches.

<a name="practice"></a>
## 8. What this means in practice

- **Coupled solvers converge in far fewer outer iterations, and their count barely depends on the mesh.** This is the reason commercial codes offer them for steady problems: Fluent's coupled pressure-based scheme and STAR-CCM+'s coupled flow solver typically allow much larger pseudo-time steps (Courant numbers) than their segregated counterparts. In the OpenFOAM world, block-coupled solvers exist in foam-extend (`pUCoupledFoam`) and in several research codes built on PETSc.
- **They pay with memory.** A factor of 8–9 in this test. On a large 3D mesh, memory per node can be what decides between the two.
- **They remove the under-relaxation guesswork.** SIMPLE's speed depended on two factors by a factor of 8, and half of the combinations diverged. Variants like SIMPLEC (available in OpenFOAM with `consistent yes`) relax this by improving the same diagonal approximation, allowing larger factors.
- **Segregated is often the better choice for transients.** In a time-accurate simulation with a small time step (PISO, PIMPLE), each time step starts from a solution that is already close, and a few segregated corrections are enough. The advantage of coupling is largest when many outer iterations are needed: steady problems, large meshes, strongly coupled physics.
- **Small problems do not need it.** On the smallest grid a well-tuned SIMPLE was just as fast.
- **The two are less different than they look.** Both rely on the same approximation of the pressure Schur complement. The difference is whether that approximation is used as the algorithm (SIMPLE) or as a preconditioner inside a Krylov method (coupled).

The general lesson is the same as in the previous articles: the equations do not change, the result does not change, and the cost can change by an order of magnitude depending on how the solution is organized. Counting iterations is not enough; time, memory and robustness have to be measured together.

### Reproduce it
All the files are available: [cavity.c](/assets/projects/OpenFOAM-coupled-solver/cavity.c) (needs LAPACK), [plot_sc.py](/assets/projects/OpenFOAM-coupled-solver/plot_sc.py) and the data: [iterative.csv](/assets/projects/OpenFOAM-coupled-solver/iterative.csv), [sweep.csv](/assets/projects/OpenFOAM-coupled-solver/sweep.csv) (direct solvers), [relaxation.csv](/assets/projects/OpenFOAM-coupled-solver/relaxation.csv), [history_iter.csv](/assets/projects/OpenFOAM-coupled-solver/history_iter.csv).

```bash
gcc -O3 -march=native cavity.c -o cavity -llapack -lm
mkdir -p out && ./cavity out        # everything, ~12 min (most of it: SIMPLE on 128 x 128)
./cavity out iter                   # only the iterative-solver comparison, ~10 min
mkdir -p figures && python3 plot_sc.py out figures
```
