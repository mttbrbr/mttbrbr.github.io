---
title: "Is My Mesh Fine Enough? Richardson Extrapolation and the GCI"
date: 2026-11-22
tech: ["C", "CFD", "Verification", "Numerical Methods"]
excerpt: "Three meshes and a formula tell you the order your solver really achieves, an estimate of the mesh-independent answer and an error band. Applied to a lid-driven cavity on five meshes from 16² to 256², including what happens when the meshes are too coarse for the method to work."
---

## Technical Index
1. [**The question**](#question) — How far is the answer from the mesh-independent one?
2. [**The method**](#method) — Observed order, Richardson extrapolation, GCI
3. [**The experiment**](#experiment) — Five meshes of the lid-driven cavity
4. [**Results**](#results) — When the GCI works and when it does not
5. [**What this means in practice**](#practice) — Running a mesh study that means something

<a name="question"></a>
## 1. The question

Every CFD result carries a discretization error: the difference between the solution on the mesh you used and the solution of the equations on an infinitely fine mesh. It cannot be measured directly, because the exact solution is unknown. The common practice, "I refined the mesh and the result did not change much", is a good instinct but not a number: how much is "not much", and how far is the result still from the limit?

The standard answer is a **grid convergence study** with **Richardson extrapolation** and the **Grid Convergence Index** (GCI), proposed by Roache in 1994 and written down as a procedure by Celik et al. in 2008 (it is the editorial policy of the ASME *Journal of Fluids Engineering*). It needs three meshes and a few lines of arithmetic.

<a name="method"></a>
## 2. The method

Take a quantity of interest $f$ (a velocity at a point, a drag coefficient, a pressure drop) computed on three meshes with cell sizes $h_1 < h_2 < h_3$ and a constant refinement ratio $r = h_2/h_1 = h_3/h_2$. If the meshes are fine enough, the discretization error behaves like $C h^p$, where $p$ is the order of the scheme. Three values are enough to find the three unknowns:

**Observed order:**
$$
p = \frac{\ln\left[(f_3 - f_2)/(f_2 - f_1)\right]}{\ln r}
$$

**Richardson extrapolation**, the estimate of the mesh-independent value:
$$
f_{ext} = f_1 + \frac{f_1 - f_2}{r^p - 1}
$$

**Grid Convergence Index** of the finest mesh, an error band with a safety factor $F_s = 1.25$:
$$
\text{GCI}_{fine} = F_s\,\frac{|(f_1 - f_2)/f_1|}{r^p - 1}
$$

The observed order is the most important number of the three. If it is close to the formal order of the scheme, the meshes are in the *asymptotic range* and the extrapolation is trustworthy. If it is not, the method is telling you that its own assumptions do not hold.

<a name="experiment"></a>
## 3. The experiment

[gci.c](/assets/projects/Grid-Convergence/gci.c) solves the steady lid-driven cavity at Re = 100 with the coupled solver from the [segregated vs coupled article](/projects/openfoam-coupled-solver/) (staggered grid, first-order upwind convection), on five meshes: 16², 32², 64², 128² and 256² cells, so $r = 2$. On each mesh it records three quantities:

- $u$ at the centre of the cavity,
- $v$ at the centre of the cavity,
- the total kinetic energy of the flow.

[gci_analysis.py](/assets/projects/Grid-Convergence/gci_analysis.py) applies the formulas above to every triplet of consecutive meshes.

The scheme is first-order upwind on purpose: its formal order is 1, so the expected observed order is known, and its [numerical diffusion](/projects/numerical-diffusion/) produces errors large enough to see.

<a name="results"></a>
## 4. Results

| Cells | u(0.5, 0.5) | v(0.5, 0.5) | Kinetic energy |
|---|---|---|---|
| 16² | −0.16684 | 0.04071 | 0.024555 |
| 32² | −0.18656 | 0.04624 | 0.029746 |
| 64² | −0.19726 | 0.05104 | 0.032210 |
| 128² | −0.20300 | 0.05406 | 0.033365 |
| 256² | −0.20601 | 0.05574 | 0.033917 |

The values keep changing, by 1.5% from 128² to 256² for $u$: the solution is clearly not mesh-independent yet. How far is it?

![u at the centre vs cell size, with the Richardson extrapolation](/assets/projects/Grid-Convergence/gci.svg)

| Quantity | Meshes | Observed order p | Extrapolated value | GCI (fine mesh) |
|---|---|---|---|---|
| u(0.5, 0.5) | 16 / 32 / 64 | 0.88 | −0.20995 | 8.0% |
| u(0.5, 0.5) | 32 / 64 / 128 | 0.90 | −0.20964 | 4.1% |
| u(0.5, 0.5) | 64 / 128 / 256 | 0.93 | −0.20935 | 2.0% |
| v(0.5, 0.5) | 16 / 32 / 64 | **0.21** | 0.0822 | **76%** |
| v(0.5, 0.5) | 32 / 64 / 128 | 0.67 | 0.0592 | 11.9% |
| v(0.5, 0.5) | 64 / 128 / 256 | 0.85 | 0.0578 | 4.7% |
| kinetic energy | 16 / 32 / 64 | 1.08 | 0.03443 | 8.6% |
| kinetic energy | 32 / 64 / 128 | 1.09 | 0.03439 | 3.8% |
| kinetic energy | 64 / 128 / 256 | 1.07 | 0.03442 | 1.9% |

*Reference values from Ghia et al. (1982), 129 × 129 points, second-order scheme: u = −0.20581, v = 0.05454.*

What the numbers say:

- **For $u$ and the kinetic energy the method works as it should.** The observed order is 0.9–1.1 on every triplet, close to the formal order 1 of upwind. The extrapolated value is almost the same whichever triplet is used (−0.2099, −0.2096, −0.2094 for $u$; 0.0344 for the energy), and the GCI halves at each refinement, as it should for a first-order method. On the finest mesh, $u = -0.2060$ with a GCI of 2%: the mesh-independent value is estimated to lie within about ±0.004.
- **For $v$ on the coarse meshes it does not.** The first triplet gives an observed order of 0.21, an extrapolation of 0.082 (about 40% above the other estimates) and a GCI of 76%. The coarsest mesh, 16 × 16, is simply too coarse for $v$ at the centre: its error does not follow $C h^p$ yet. Only on the finer triplets does the order approach 1 and the extrapolation settle around 0.058. The lesson is that the observed order is a diagnostic: when it is far from the formal order, the GCI is not a reliable error band.
- **The comparison with the reference is instructive.** The extrapolated $u$ is −0.2094; the classic reference value of Ghia et al. is −0.2058, 1.7% away, and lies inside the GCI band of the finest mesh. The reference is itself a numerical solution on a 129 × 129 grid with its own discretization error, so a disagreement of this size says as much about the reference as about the extrapolation. Even "benchmark" numbers come with an error bar.

There is also a practical cost: the 256² mesh took **14 minutes** with this solver, against 35 seconds for 128² and 2 seconds for 64². A first-order scheme needs a very fine mesh to get a small error; the same study with a second-order scheme would reach a 2% band on a much coarser mesh.

<a name="practice"></a>
## 5. What this means in practice

- **Use at least three meshes with a constant refinement ratio.** In 3D a ratio of 2 means 8 times more cells per mesh, so ratios between 1.3 and 1.5 are common; the procedure works with any $r > 1.3$, as long as it is applied consistently (refine the whole mesh, not just some regions).
- **Look at the observed order before the GCI.** If it is far from the formal order of your schemes, or changes a lot between quantities or triplets, the meshes are not in the asymptotic range, or something else is going on (limiters switching, a turbulence model reacting to $y^+$, non-converged iterations). In that case the GCI is just a number.
- **Converge the iterations first.** An iteration error larger than the discretization error makes the differences between meshes meaningless. Here every mesh was converged by $10^8$.
- **Report the quantities that matter.** The method is applied quantity by quantity: a mesh can be converged for the drag and not for a local velocity. Integral quantities (forces, energy, flow rates) usually converge more regularly than point values.
- **Remember what it does not cover.** The GCI estimates the discretization error only. It says nothing about modelling errors (turbulence model, boundary conditions, geometry), and it cannot catch an *inconsistent* discretization: as the [non-orthogonality article](/projects/non-orthogonality/) showed, an uncorrected scheme can converge beautifully to the wrong answer.

The theme of the series once more: a result without an error estimate is a guess. Three meshes and a formula turn "looks converged" into "$-0.206 \pm 0.004$", and occasionally they tell you that you do not know as much as you thought.

### Reproduce it
All the files are available: [gci.c](/assets/projects/Grid-Convergence/gci.c) (needs LAPACK), [gci_analysis.py](/assets/projects/Grid-Convergence/gci_analysis.py) with [svgchart.py](/assets/projects/Grid-Convergence/svgchart.py), and the data [grids.csv](/assets/projects/Grid-Convergence/grids.csv), [gci.csv](/assets/projects/Grid-Convergence/gci.csv).

```bash
gcc -O3 -march=native gci.c -o gci -llapack -lm
mkdir -p out && ./gci out 256          # ~15 min, almost all of it on 256 x 256
./gci out 128                           # ~40 s, without the finest mesh
mkdir -p figures && python3 gci_analysis.py out figures
```
