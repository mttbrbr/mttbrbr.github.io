---
title: "Gauss linear or leastSquares? Gradient Schemes on Imperfect Meshes"
date: 2026-11-22
tech: ["C", "CFD", "Meshing", "OpenFOAM"]
excerpt: "Two lines in every fvSchemes file decide how cell gradients are computed. On a regular mesh the choice does not matter. On a skewed or irregular one, the Gauss gradient keeps a 12-42% error at any resolution while least squares converges. Measured on five kinds of mesh, plus what a limiter costs."
---

## Technical Index
1. [**Why gradients matter**](#why) — They are inside almost every term
2. [**Two ways to compute a gradient**](#schemes) — Gauss and least squares
3. [**The experiment**](#experiment) — Five meshes, one exact field
4. [**Results**](#results) — Where each scheme fails
5. [**The cost of a limiter**](#limiter) — `cellLimited` on a smooth field
6. [**What this means in practice**](#practice) — Choosing `gradSchemes`

<a name="why"></a>
## 1. Why gradients matter

A finite-volume solver stores one value per cell, but it needs gradients almost everywhere: the pressure gradient in the momentum equation, the reconstruction of face values in second-order convection schemes (`linearUpwind`, MUSCL-type schemes), the explicit part of the [non-orthogonal correction](/projects/non-orthogonality/), the production term of turbulence models, wall distances, limiters. An error in the cell gradient propagates into all of them.

In OpenFOAM the choice is the `gradSchemes` entry of `fvSchemes`: usually `Gauss linear`, sometimes `leastSquares`, possibly with a limiter. This article compares them on meshes that are regular, skewed, irregular and stretched.

<a name="schemes"></a>
## 2. Two ways to compute a gradient

- **Gauss (Green–Gauss).** The divergence theorem turns the volume integral of the gradient into a sum over faces: $\nabla\phi_P = \frac{1}{V}\sum_f \phi_f\,\mathbf{S}_f$. The face values $\phi_f$ come from interpolation, here linear with OpenFOAM's geometric weights. Its accuracy depends on how well $\phi_f$ represents the value at the face centre, which is exactly what [skewness](/projects/skewness/) breaks.
- **Least squares.** Fit a linear function through the cell and its neighbours: find the gradient $\mathbf{g}$ that minimizes $\sum_f w_f\,(\phi_N - \phi_P - \mathbf{d}_f\cdot\mathbf{g})^2$, where $\mathbf{d}_f$ is the vector between the cell centres and $w_f = 1/|\mathbf{d}_f|^2$ (OpenFOAM's weights). It uses only cell centres, not face values. For any linear field it is exact, on any mesh.

<a name="experiment"></a>
## 3. The experiment

[grad.c](/assets/projects/Gradient-Schemes/grad.c) computes both gradients on periodic quadrilateral meshes, for a smooth field whose exact gradient is known, on 16² to 256² cells. Five meshes:

- **uniform**: squares;
- **skewed**: rows of parallelograms leaning alternately left and right (the mesh of the skewness article, skewness 0.3);
- **random**: every node moved randomly by up to 25% of the cell size, a stand-in for an irregular unstructured mesh;
- **stretched**: cells 50 times longer than they are thick, as in a boundary layer;
- **stretched and curved**: the same thin cells following a wavy surface, like a boundary layer over a bump.

<a name="results"></a>
## 4. Results

![Gradient error vs mesh size for Gauss linear and leastSquares](/assets/projects/Gradient-Schemes/gradients.svg)

*On the uniform and stretched meshes the two schemes give identical errors (in the left panel the uniform curve is hidden under the curved one).*

| Mesh | Gauss linear, 32² | Gauss linear, 256² | leastSquares, 32² | leastSquares, 256² |
|---|---|---|---|---|
| uniform | 0.64% | 0.010% | 0.64% | 0.010% |
| stretched (50:1) | 0.64% | 0.010% | 0.64% | 0.010% |
| stretched + curved | 0.64% | 0.010% | 0.74% | 0.012% |
| skewed | **42%** | **42%** | 0.64% | 0.010% |
| random | **12%** | **12%** | 1.8% | 0.21% |

*Relative L2 error of the gradient.*

- **On regular meshes the choice does not matter.** On squares and on straight stretched cells the two schemes are identical and second order.
- **Gauss linear does not converge on skewed or irregular meshes.** Its error stays at 42% on the skewed mesh and 12% on the random one, whatever the resolution: the same inconsistency as in the skewness article. On an irregular mesh every face is a little skewed, in random directions, and the errors do not cancel.
- **Least squares converges everywhere**: second order on the skewed mesh (it does not see skewness at all), first order on the random one (the error halves at each refinement), still 60 times more accurate than Gauss on the finest random mesh.
- **On thin curved cells least squares is slightly worse** in this test (about 20%), but still second order. The literature reports much larger least-squares errors on highly stretched cells around strongly curved walls (Mavriplis, 2003), where neighbours along the wall are far away and far from the tangent plane; this smooth, mildly curved test does not reproduce that worst case.

<a name="limiter"></a>
## 5. The cost of a limiter

Gradient limiters (`cellLimited`, `cellMDLimited`, `faceLimited`) scale the gradient down so that values extrapolated to the faces do not exceed the range of the neighbouring cells. They prevent new extrema and keep bounded quantities bounded, which is essential near shocks, interfaces and steep fronts. But they cannot tell a steep front from an ordinary smooth maximum.

![Least squares with and without cellLimited](/assets/projects/Gradient-Schemes/limiter.svg)

| Mesh | leastSquares, 256² | with cellLimited 1, 256² |
|---|---|---|
| uniform | 0.010% | 0.040% |
| random | 0.21% | 0.22% |

On the smooth field the limiter makes the gradient error **4 times larger** on the regular mesh: at every maximum and minimum of the field it flattens the gradient, because any non-zero gradient would extrapolate beyond the neighbours. The order of convergence does not change, so the penalty is a constant factor; on the irregular mesh, where the error is dominated by the mesh, the limiter barely adds to it.

<a name="practice"></a>
## 6. What this means in practice

- **On good hexahedral meshes, `Gauss linear` is fine.** It is cheap, robust and as accurate as least squares.
- **On tetrahedral, polyhedral or skewed meshes, consider `leastSquares`** (or `pointCellsLeastSquares`, which uses all point neighbours) for the gradients that matter, typically the velocity gradient feeding `linearUpwind` and the turbulence production. Alternatively, `Gauss skewCorrected linear` removes the skewness error at a higher cost.
- **Watch thin cells on curved walls.** Least squares can lose accuracy in high-aspect-ratio layers on strongly curved surfaces; this is why some codes switch to Gauss gradients in boundary layers.
- **Use limiters where they are needed, not everywhere.** `cellLimited Gauss linear 1` on velocity is a common robustness choice, but it costs accuracy at every smooth extremum. `cellMDLimited` (limiting each direction separately) and coefficients below 1 (for example 0.5) are gentler compromises.
- **A gradient error is not visible in the residuals.** The solver converges happily with an inconsistent gradient; only a comparison with a known solution, or a careful mesh study, reveals it.

### Reproduce it
All the files are available: [grad.c](/assets/projects/Gradient-Schemes/grad.c), [plot_gr.py](/assets/projects/Gradient-Schemes/plot_gr.py) with [svgchart.py](/assets/projects/Gradient-Schemes/svgchart.py), and the data [gradients.csv](/assets/projects/Gradient-Schemes/gradients.csv).

```bash
gcc -O3 -march=native grad.c -o grad -lm
mkdir -p out && ./grad out                  # < 1 s
mkdir -p figures && python3 plot_gr.py out figures
```
