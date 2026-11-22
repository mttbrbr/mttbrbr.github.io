---
title: "Mesh Non-Orthogonality: What checkMesh Is Warning You About"
date: 2026-11-22
tech: ["C", "CFD", "OpenFOAM", "Meshing"]
excerpt: "On a non-orthogonal mesh, an uncorrected Laplacian does not just lose accuracy: it converges to the wrong answer, and refining the mesh does not help. Measured on a sheared mesh with an exact solution, with OpenFOAM's orthogonal, uncorrected and corrected schemes."
---

## Technical Index
1. [**What non-orthogonality is**](#what) — The angle checkMesh reports
2. [**The correction**](#correction) — Splitting the face gradient in two
3. [**The experiment**](#experiment) — A sheared mesh with an exact solution
4. [**Error vs angle**](#angle) — How much the correction matters
5. [**Error vs mesh size**](#grid) — The error that refinement does not fix
6. [**How many correctors**](#correctors) — `nNonOrthogonalCorrectors`
7. [**What this means in practice**](#practice) — Schemes, limits and meshing

<a name="what"></a>
## 1. What non-orthogonality is

Every finite-volume solver computes fluxes through faces, and diffusive fluxes (viscous stresses, heat conduction, and above all the pressure equation) need the gradient of a variable normal to the face. The simplest estimate uses the two cell centres on either side:

$$
\nabla u \cdot \mathbf{S}_f \approx |\mathbf{S}_f|\,\frac{u_N - u_P}{|\mathbf{d}|}
$$

where $\mathbf{S}_f$ is the face area vector (normal to the face, as long as the face is large) and $\mathbf{d}$ is the vector from the centre of cell $P$ to the centre of cell $N$. This is exact for a linear field only if $\mathbf{d}$ is parallel to $\mathbf{S}_f$. The angle between them is the **non-orthogonality** of the face, the number `checkMesh` reports as *Mesh non-orthogonality Max* and flags above 70°.

On a hexahedral mesh aligned with the geometry it is close to zero. On tetrahedral and polyhedral meshes, around refinement interfaces, and in distorted cells near complex surfaces, it is typically 30–70°, sometimes more.

<a name="correction"></a>
## 2. The correction

The standard fix (as implemented in OpenFOAM, following Jasak's thesis, 1996) splits the face vector in two parts:

$$
\mathbf{S}_f = \boldsymbol{\Delta}_f + \mathbf{k}_f, \qquad \boldsymbol{\Delta}_f = \mathbf{d}\,\frac{\mathbf{S}_f\cdot\mathbf{S}_f}{\mathbf{d}\cdot\mathbf{S}_f}
$$

$\boldsymbol{\Delta}$ is parallel to $\mathbf{d}$, so its part of the flux can be computed exactly from the two cell values and enters the matrix **implicitly**. The remainder $\mathbf{k}$ needs the full gradient on the face, interpolated from the cell gradients, and is added **explicitly** to the right-hand side:

$$
\nabla u \cdot \mathbf{S}_f \approx \underbrace{|\boldsymbol{\Delta}_f|\,\frac{u_N - u_P}{|\mathbf{d}|}}_{\text{implicit}} + \underbrace{\mathbf{k}_f \cdot (\nabla u)_f}_{\text{explicit}}
$$

Because the explicit part uses the current solution, the equation has to be solved again with the updated correction: this is what `nNonOrthogonalCorrectors` controls. In OpenFOAM the choice appears in the `snGradSchemes` and `laplacianSchemes` of `fvSchemes`:

- `orthogonal`: ignore the angle and use $|\mathbf{S}|/|\mathbf{d}|$;
- `uncorrected`: keep only the implicit $\boldsymbol{\Delta}$ part;
- `corrected`: implicit part plus the explicit correction;
- `limited corrected ψ`: the explicit correction, limited to a fraction of the implicit part (ψ = 0 is `uncorrected`, ψ = 1 is `corrected`).

<a name="experiment"></a>
## 3. The experiment

[nonortho.c](/assets/projects/Non-Orthogonality/nonortho.c) solves the Poisson equation $-\nabla^2 u = f$ on a mesh of parallelograms, obtained by shearing a square grid:

$$
x = \xi + \tan\theta\;\eta, \qquad y = \eta
$$

Every face of this mesh has exactly the same non-orthogonality, the angle $\theta$, so a single number describes the mesh. The domain is periodic, so there are no boundary conditions to blur the picture, and the source $f$ is chosen so that the exact solution is known: $u = \sin(2\pi\xi)\cos(2\pi\eta)$. The error is the relative L2 difference between the computed and the exact solution at the cell centres. Each case takes a fraction of a second.

<a name="angle"></a>
## 4. Error vs angle

![Relative error vs non-orthogonality on 64 x 64 cells](/assets/projects/Non-Orthogonality/error_angle.svg)

| Angle | orthogonal | uncorrected | 1 corr. | 2 corr. | converged |
|---|---|---|---|---|---|
| 0° | 0.08% | 0.08% | 0.08% | 0.08% | 0.08% |
| 10° | 18% | 17% | 2.9% | 0.57% | 0.08% |
| 30° | 59% | 50% | 24% | 12% | 0.16% |
| 50° | 122% | 70% | 48% | 34% | 0.36% |
| 70° | 256% | 58% | 33% | 19% | 0.21% |

*Relative L2 error on 64 × 64 cells. The last three columns use the `corrected` scheme: with 1 non-orthogonal corrector, with 2, and with the explicit correction iterated until it stops changing.*

- At **0°** all schemes are identical and the error is 0.08%: the pure discretization error of a second-order scheme on this mesh.
- Already at **10°**, ignoring the angle gives an error of **17–18%**, more than 200 times larger.
- **`orthogonal` is the worst**: it gets the size of the coefficient wrong as well as the direction, and its error grows without limit with the angle.
- **`uncorrected`** is better but still wrong by 50–70% in the middle range. It drops exactly the term that couples the two directions of the mesh (a mixed derivative $\partial^2 u/\partial\xi\,\partial\eta$), and how much that matters depends on the solution: with this solution the error peaks around 50–60°.
- **`corrected`**, once the correction has converged, stays within 0.08–0.36% at every angle. Non-orthogonality costs a little accuracy, not correctness.

<a name="grid"></a>
## 5. Error vs mesh size

The most important result is what happens when the mesh is refined.

![Relative error vs mesh size](/assets/projects/Non-Orthogonality/error_grid.svg)

| Cells per side | 0°, any scheme | 30°, uncorrected | 60°, uncorrected | 60°, corrected |
|---|---|---|---|---|
| 16 | 1.3% | 50.1% | 70.2% | 5.2% |
| 32 | 0.32% | 49.7% | 69.5% | 1.4% |
| 64 | 0.080% | 49.5% | 69.3% | 0.35% |
| 128 | 0.020% | 49.5% | 69.3% | 0.088% |
| 256 | 0.0050% | 49.5% | 69.3% | 0.022% |

*Relative L2 error. The non-orthogonality angle is the same on every mesh.*

With the correction, the error falls by 4 every time the mesh is refined, exactly as on an orthogonal mesh: second-order convergence. **Without it, the error does not fall at all.** On 256 × 256 cells the uncorrected solution is still 50% or 70% wrong, exactly as on 16 × 16.

This is the difference between an *inaccurate* scheme and an *inconsistent* one. An inaccurate scheme gets better with resolution. An inconsistent one converges, very precisely, to the solution of a different equation: here a Laplacian with the cross-derivative term missing. A mesh refinement study would show a perfectly converged result, and it would be wrong.

The mesh in this test is uniformly skewed, which is the worst case. On a real mesh the non-orthogonality varies from face to face and the high values are often local, so the error is smaller and concentrated. But it does not go away with refinement if the refined mesh keeps the same shape of cells, which is exactly what happens when a tetrahedral mesh is refined.

<a name="correctors"></a>
## 6. How many correctors

The explicit correction uses the gradient of the current solution, so it has to be iterated. How fast does it converge?

![Error vs number of non-orthogonal correctors](/assets/projects/Non-Orthogonality/error_correctors.svg)

| Angle | Correctors to converge | Error after 1 corrector | Error converged |
|---|---|---|---|
| 10° | 14 | 2.9% | 0.085% |
| 30° | 33 | 24% | 0.16% |
| 60° | 63 | 48% | 0.35% |
| 80° | 21 | 11% | 0.10% |

*Correctors needed until the solution changes by less than $10^{-10}$, 64 × 64 cells.*

Here each corrector is a full solve of the equation, and one or two of them are far from enough: at 30° the error is still 24% after one. It would seem to say that `nNonOrthogonalCorrectors 1` is useless. In a real solver it is not, for one reason: **the outer iterations also update the correction.**

- In a **steady** SIMPLE run, the pressure equation is solved again at every outer iteration with the latest solution, so the explicit correction is updated hundreds or thousands of times before the run converges. At convergence it is fully converged too, whatever `nNonOrthogonalCorrectors` says. The extra correctors mainly help stability and the rate of convergence.
- In a **transient** PISO or PIMPLE run with only a few corrections per time step, there is no such luxury: the correction has to converge within the time step, or the error stays in the solution. Here `nNonOrthogonalCorrectors` directly affects accuracy, and meshes with high non-orthogonality need more of them.

<a name="practice"></a>
## 7. What this means in practice

- **Use `corrected` (or `limited corrected`) on any mesh that is not orthogonal.** `orthogonal` and `uncorrected` are only safe on meshes that are, in fact, orthogonal. Their error does not go away with refinement.
- **High angles make the correction unstable**, because the explicit part grows like $\tan\theta$ and can dominate the implicit one. This is why `limited corrected 0.33` or `0.5` is a common choice when non-orthogonality exceeds about 70°: it trades some accuracy for robustness, between the two extremes measured here.
- **Steady and transient runs need different settings.** In steady SIMPLE, 0–1 correctors are usually enough because the outer iterations update the correction anyway. In transient runs on poor meshes, use 1–3 and check that the results do not change with more.
- **The real fix is the mesh.** Every corrector costs another pressure solve, the most expensive step of the iteration ([see the multigrid article](/projects/multigrid/)). Aligned cells, smooth size transitions and good prism layers remove the problem instead of patching it.
- **A converged mesh study is not proof of correctness.** As the 60° results show, an inconsistent discretization converges beautifully to the wrong answer.

As in the rest of this series: the equations are fine, the solver is fine, and the answer can still be wrong, here by 70%, because of one geometric detail of the discretization.

### Reproduce it
All the files are available: [nonortho.c](/assets/projects/Non-Orthogonality/nonortho.c), [plot_no.py](/assets/projects/Non-Orthogonality/plot_no.py) with [svgchart.py](/assets/projects/Non-Orthogonality/svgchart.py), and the data [angle.csv](/assets/projects/Non-Orthogonality/angle.csv), [grid.csv](/assets/projects/Non-Orthogonality/grid.csv), [correctors.csv](/assets/projects/Non-Orthogonality/correctors.csv).

```bash
gcc -O3 -march=native nonortho.c -o nonortho -lm
mkdir -p out && ./nonortho out          # ~10 s
mkdir -p figures && python3 plot_no.py out figures
```
