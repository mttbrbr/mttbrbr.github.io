---
title: "Mesh Expansion Ratio: When Cell-Size Jumps Matter"
date: 2026-11-22
tech: ["C", "CFD", "Meshing", "OpenFOAM"]
excerpt: "A 2:1 jump in cell size wrecks the local truncation error, and the solution of a diffusion problem does not care. A wave crossing the same jump is partly reflected as grid-scale noise, a hundred times more than across a gradual transition. Two 1D experiments on what the expansion-ratio rules of thumb are really about."
---

## Technical Index
1. [**The rule of thumb**](#rule) — Keep neighbouring cells within 20%
2. [**The experiment**](#experiment) — Abrupt jumps and gradual transitions
3. [**Diffusion: the jump that does not matter**](#diffusion) — Truncation error vs solution error
4. [**Waves: the jump that does**](#waves) — Reflections at the transition
5. [**What this means in practice**](#practice) — Refinement levels, layers and LES

<a name="rule"></a>
## 1. The rule of thumb

Every meshing guide has some version of it: neighbouring cells should not differ in size by more than about 20–30%. Prism layers grow with an `expansionRatio` of 1.1–1.3; `blockMesh` grading should be gentle; refinement regions in `snappyHexMesh` create 2:1 jumps at every level, buffered by `nCellsBetweenLevels`.

The reason usually given is accuracy: on a non-uniform mesh the standard second-order formulas lose an order where the cell size changes. That is true, but the consequences are less simple than the rule suggests. This article measures them in the cleanest possible setting, one dimension, with two very different kinds of problem.

<a name="experiment"></a>
## 2. The experiment

[expansion.c](/assets/projects/Expansion-Ratio/expansion.c) builds 1D cell-centred finite-volume meshes that go from a fine region (cell size $h$) to a coarse one (cell size $2h$ or $4h$) in one of two ways:

- **abruptly**: one jump, like the interface between two refinement levels;
- **gradually**: with a constant expansion ratio $r = 1.1$ or $1.2$ between neighbouring cells, until the coarse size is reached.

On these meshes it solves:

- **A. a diffusion problem**, $-u'' = f$ with a known exact solution, the 1D analogue of a pressure or heat-conduction equation;
- **B. a wave problem**, $u_t + u_x = 0$ with central fluxes (linear interpolation with the proper geometric weights, as OpenFOAM's `linear` scheme), where a wave packet crosses the transition from fine to coarse.

<a name="diffusion"></a>
## 3. Diffusion: the jump that does not matter

Two different errors can be measured. The **truncation error** is what you get by plugging the exact solution into the discrete equations: how badly the scheme approximates the equation, cell by cell. The **solution error** is the difference between the computed and the exact solution, which is what you actually care about.

![Truncation error and solution error vs cell size](/assets/projects/Expansion-Ratio/diffusion.svg)

| Mesh | Trunc., h = 1/40 | Trunc., h = 1/1280 | Error, h = 1/1280 |
|---|---|---|---|
| uniform | 0.20 | 0.0002 | 1.6·10⁻⁶ |
| jump 2:1 | 14 | 13 | 6.5·10⁻⁶ |
| graded 1.2 to 2× | 11 | 2.6 | 6.6·10⁻⁶ |
| graded 1.1 to 2× | 10 | 1.3 | 6.6·10⁻⁶ |
| jump 4:1 | 47 | 39 | 2.7·10⁻⁵ |
| graded 1.2 to 4× | 4.2 | 2.6 | 2.8·10⁻⁵ |

*Trunc.: maximum local truncation error. Error: L2 error of the solution. h: fine cell size.*

The truncation error at a transition **does not decrease at all** when the mesh is refined: on the 2:1 jump it stays at 13 whether there are 28 or 900 cells. The scheme is locally inconsistent there. Gradual transitions reduce it by 5 to 15 times, but it does not vanish either.

And yet the solution error decreases by 4 at every refinement on *every* mesh, exactly like on the uniform one: second-order convergence. With the same fine cell size, the jump and the graded transitions give the same error, about 6.6·10⁻⁶; the meshes that go to 4× are four times less accurate simply because their coarse region is coarser.

This is a known property of cell-centred finite volumes for diffusion, called **supraconvergence**: the local errors at the transition are large but they cancel out in the solution, because they have the form of a difference of fluxes that telescopes across neighbouring cells. For elliptic problems like the pressure equation, a cell-size jump costs much less accuracy than the truncation error suggests.

<a name="waves"></a>
## 4. Waves: the jump that does

Waves are a different story. A wave packet travels from the fine region into the coarse one. On a perfect scheme it would pass the transition untouched; with central fluxes on a non-uniform mesh, part of it is **reflected backwards** as spurious waves.

![Fraction of the wave energy reflected by the transition](/assets/projects/Expansion-Ratio/waves.svg)

| Transition | Reflected, 16 cells/λ | Reflected, 32 cells/λ |
|---|---|---|
| uniform mesh | 10⁻¹⁴ | 10⁻¹⁷ |
| jump 2:1 | 1.4·10⁻³ | 1.0·10⁻⁵ |
| graded 1.2 to 2× | 2.1·10⁻⁴ | 6.2·10⁻⁷ |
| graded 1.1 to 2× | 6.3·10⁻⁵ | 3.1·10⁻⁸ |
| jump 4:1 | wave destroyed | 5.3·10⁻⁴ |
| graded 1.2 to 4× | wave destroyed | 5.4·10⁻⁵ |

*Fraction of the wave energy reflected. Cells per wavelength (cells/λ) counted on the fine side; the coarse side has half (2×) or a quarter (4×) of them.*

- **A gradual transition reflects much less than a jump**: 7 to 22 times less at 16 cells per wavelength, 16 to 330 times less at 32. The smoother the change in cell size, the smaller the reflection.
- **What is reflected is grid-scale noise.** The waves left behind are zig-zags two cells long, the parasitic mode that central schemes support and do not damp:

![Waves left behind the transition](/assets/projects/Expansion-Ratio/reflections.svg)

- **The coarse side must still resolve the wave.** With 16 cells per wavelength on the fine side and a 4:1 jump, only 4 remain on the coarse side: the wave cannot propagate there at all and is destroyed, whatever the transition looks like. No grading fixes a mesh that is too coarse where the wave is going.
- **Energy is not conserved across the transition.** Even when the wave passes, 3–18% of its energy is lost, depending on the resolution. With geometric interpolation weights, the central scheme is no longer exactly energy-conserving on a non-uniform mesh, a known issue for LES codes that rely on non-dissipative schemes.

<a name="practice"></a>
## 5. What this means in practice

- **For steady RANS and diffusion-dominated quantities, refinement jumps are rarely the main error.** The solution converges at the expected order across them. The expansion-ratio rule matters more for robustness and for the quality of the other mesh metrics (skewness, non-orthogonality) that sudden jumps usually bring along.
- **For anything that propagates, it matters a lot.** Acoustics, LES, vortex shedding, wakes crossing refinement interfaces: every abrupt jump acts like a partial mirror and a source of grid-scale noise. Here a gradual transition reduces the reflection by one to two orders of magnitude.
- **Buffer the refinement levels.** In `snappyHexMesh`, `nCellsBetweenLevels` (typically 2–4) keeps several cells of each size between levels, so that each 2:1 jump is isolated and the overall transition is gradual. Avoid placing level jumps in the region where the flow features you care about travel.
- **Resolve the flow on both sides.** The worst case is not the jump itself but a coarse region that cannot carry what arrives from the fine one. Check how many cells per wavelength (or per vortex diameter) the coarse side has.
- **Upwind-biased schemes hide the problem by damping it.** The reflected zig-zags are exactly the modes that [numerical diffusion](/projects/numerical-diffusion/) removes first. That makes a dissipative scheme look cleaner at interfaces, and it is one more reason why LES with central schemes needs smooth meshes.

The common thread of the series: the same mesh defect can be harmless for one equation and important for another. A rule of thumb is a summary; knowing what it protects against tells you when you can bend it.

### Reproduce it
All the files are available: [expansion.c](/assets/projects/Expansion-Ratio/expansion.c), [plot_ex.py](/assets/projects/Expansion-Ratio/plot_ex.py) with [svgchart.py](/assets/projects/Expansion-Ratio/svgchart.py), and the data [diffusion.csv](/assets/projects/Expansion-Ratio/diffusion.csv), [waves.csv](/assets/projects/Expansion-Ratio/waves.csv).

```bash
gcc -O3 -march=native expansion.c -o expansion -lm
mkdir -p out && ./expansion out            # < 1 s
mkdir -p figures && python3 plot_ex.py out figures
```
