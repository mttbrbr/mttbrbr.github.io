---
title: "Mesh Skewness: The Interpolation Error That Refinement Does Not Fix"
date: 2026-11-22
tech: ["C", "CFD", "Meshing", "OpenFOAM"]
excerpt: "When a face centre is not where the line between two cell centres crosses the face, linear interpolation evaluates the field at the wrong point. On a skewed mesh the Gauss gradient keeps an error of 14% at skewness 0.1 whatever the resolution; the skew correction brings it back to second order."
---

## Technical Index
1. [**What skewness is**](#what) — Not the same as non-orthogonality
2. [**The experiment**](#experiment) — A mesh with pure skewness
3. [**Results**](#results) — An error that does not converge
4. [**Why**](#why) — Interpolating at the wrong point
5. [**What this means in practice**](#practice) — `skewCorrected`, least squares and the mesh

<a name="what"></a>
## 1. What skewness is

In the [non-orthogonality article](/projects/non-orthogonality/) the problem was the *angle* between the face normal and the line joining two cell centres. Skewness is a different defect, often confused with it: it is about *where* that line crosses the face.

A finite-volume code needs values on the faces: to compute convective fluxes, and to compute cell gradients with the Gauss theorem, $\nabla\phi_P = \frac{1}{V}\sum_f \phi_f \mathbf{S}_f$. With linear interpolation, the face value is a weighted average of the two cell values, and it is the value of the field at the point where the line joining the two centres crosses the face. If that point coincides with the face centre, the integral over the face is second-order accurate. If it does not, the value is taken at the wrong point. The distance between the two points, relative to the distance between the cell centres, is the **skewness** of the face.

Skewness typically appears where cells of very different shape or size meet: at refinement interfaces, in tetrahedral and polyhedral meshes, and in distorted cells near complex surfaces.

<a name="experiment"></a>
## 2. The experiment

[skew.c](/assets/projects/Skewness/skew.c) builds a periodic 2D mesh of parallelograms whose rows lean alternately to the left and to the right. On this mesh every face is **orthogonal** (no non-orthogonality at all), but on the horizontal faces the face centre is shifted sideways from the line joining the cell centres, by a fraction $a$ of the cell size: pure skewness. (OpenFOAM's `checkMesh` uses a differently normalized skewness measure; here $a$ is simply the displacement divided by the distance between cell centres.)

On this mesh the program evaluates, for the smooth field $\phi = \sin(2\pi x)\sin(2\pi y)$, two operators every solver computes all the time:

- the **Gauss gradient** of $\phi$;
- the **convective term** $\nabla\cdot(\mathbf{U}\phi)$ with $\mathbf{U} = (1, 0.5)$;

with face values from plain **linear** interpolation, or with the **skew correction** $\phi_f = \phi_{lin} + (\mathbf{x}_f - \mathbf{x}_i)\cdot(\nabla\phi)_f$, which moves the interpolated value from the crossing point $\mathbf{x}_i$ to the face centre $\mathbf{x}_f$ (OpenFOAM's `skewCorrected linear`). Since the exact answers are known, the errors can be measured directly, on meshes from 16² to 512² cells.

<a name="results"></a>
## 3. Results

![Gauss gradient error vs mesh size, with and without skew correction](/assets/projects/Skewness/skewness.svg)

| Skew | Grad. lin. 32² | Grad. lin. 512² | Grad. corr. 512² | Conv. lin. 512² |
|---|---|---|---|---|
| 0 | 0.64% | 0.0025% | 0.0025% | 0.0025% |
| 0.05 | 7.1% | 7.1% | 0.0025% | 4.5% |
| 0.1 | 14% | 14% | 0.0026% | 8.9% |
| 0.2 | 28% | 28% | 0.0027% | 18% |
| 0.4 | 56% | 57% | 0.0033% | 36% |

*Relative L2 errors of the Gauss gradient (Grad.) and of the convective term (Conv.), with linear (lin.) or skew-corrected (corr.) face values. The convective term with skew correction behaves like the corrected gradient (0.0025–0.0028% on 512²).*

- On the regular mesh the Gauss gradient is second order: the error drops by 4 at each refinement, to 0.0025% on 512².
- With linear interpolation on the skewed mesh, the error is about **1.4 times the skewness** (14% at 0.1, 57% at 0.4) and **does not decrease at all with refinement**: it is the same on 32² and on 512² cells. The convective term behaves the same way, with an error of about 0.9 times the skewness.
- With the skew correction, both operators are second order again, with errors within 30% of the regular mesh at every skewness.

As with non-orthogonality, this is not a loss of accuracy but a loss of **consistency**: the discrete operator converges, very precisely, to the wrong one.

<a name="why"></a>
## 4. Why

On a horizontal face the interpolated value is the field at the crossing point, $a h$ away from the face centre, so it carries an error of about $a h\,\partial\phi/\partial x$. The error is of the order of the cell size, which sounds harmless. But the Gauss gradient divides the sum of face values times face areas by the cell volume: an error of order $h$ in the face values becomes an error of order $h \cdot h / h^2 = 1$ in the gradient. The faces on opposite sides of the cell are skewed in opposite directions here, so their errors add up instead of cancelling.

This mesh is a worst case on purpose. On real meshes the skewness changes from face to face, and part of the error cancels between neighbouring faces; the error is then often first order rather than constant. But wherever skewness is high and systematic, as at refinement interfaces or along layers of distorted cells, the same mechanism is at work.

<a name="practice"></a>
## 5. What this means in practice

- **`Gauss linear` is not skew-corrected.** The default gradient and interpolation schemes in most OpenFOAM tutorials use plain linear interpolation. On meshes with significant skewness, the gradients that feed the pressure equation, the limiters and the non-orthogonal correction carry the error measured here.
- **The fix exists.** `skewCorrected linear` as interpolation scheme (and in `Gauss skewCorrected linear` for gradients) adds the correction term, at the cost of computing a gradient first. It is more expensive and can be less robust on very poor meshes, which is why it is not the default.
- **Least-squares gradients are immune by construction.** The `leastSquares` gradient uses only cell centres and their distances, not face values, so face skewness does not enter it at all. This is the subject of the [gradient schemes article](/projects/gradient-schemes/).
- **The best fix is still the mesh.** Watch the skewness reported by `checkMesh`, especially at refinement interfaces and in layers. A smooth mesh does not need a correction.
- **A converged mesh study is not enough.** As with non-orthogonality, an inconsistent discretization converges with refinement, just to the wrong answer.

### Reproduce it
All the files are available: [skew.c](/assets/projects/Skewness/skew.c), [plot_sk.py](/assets/projects/Skewness/plot_sk.py) with [svgchart.py](/assets/projects/Skewness/svgchart.py), and the data [skewness.csv](/assets/projects/Skewness/skewness.csv).

```bash
gcc -O3 -march=native skew.c -o skew -lm
mkdir -p out && ./skew out                  # < 1 s
mkdir -p figures && python3 plot_sk.py out figures
```
