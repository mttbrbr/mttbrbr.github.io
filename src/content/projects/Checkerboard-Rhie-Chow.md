---
title: "Checkerboards, Rhie-Chow and Spurious Currents"
date: 2026-11-22
tech: ["C", "CFD", "Numerical Methods", "OpenFOAM"]
excerpt: "Storing velocity and pressure in the same place makes a solver simple, and lets a zig-zag pressure field pass completely unnoticed. A small collocated solver shows the checkerboard, how Rhie-Chow interpolation removes it, and why interFoam needs one more trick to avoid spurious currents."
---

## Technical Index
1. [**Where to store the pressure**](#where) — Staggered and collocated grids
2. [**The blind spot of central differences**](#blind) — A pressure that the momentum equation cannot see
3. [**The experiment**](#experiment) — A collocated solver with three face interpolations
4. [**A fluid at rest**](#hydrostatic) — The test that exposes everything
5. [**Refining the mesh**](#refine) — Which errors go away and how fast
6. [**What this means in practice**](#practice) — OpenFOAM, `p_rgh` and VOF

<a name="where"></a>
## 1. Where to store the pressure

A finite-volume solver for incompressible flow has to decide where velocity and pressure live. There are two classic answers.

- **Staggered grid** (Harlow and Welch, 1965): pressure at the cell centres, each velocity component on the faces it crosses. This is what I used in the [segregated vs coupled article](/projects/openfoam-coupled-solver/). It is very robust, but it gets awkward on unstructured meshes, where faces have arbitrary orientations and every component would need its own control volumes.
- **Collocated grid**: velocity and pressure both at the cell centres. Every variable shares the same mesh, the same control volumes and the same data structures. This is what OpenFOAM, Fluent, STAR-CCM+ and practically every unstructured code do.

The collocated grid is much simpler to implement on arbitrary meshes, but it comes with a famous problem: the **checkerboard**.

<a name="blind"></a>
## 2. The blind spot of central differences

On a collocated grid, the pressure force on cell $P$ is computed from the pressures on its faces, and with linear interpolation those are averages of the neighbouring cells. In one direction:

$$
\left(\frac{\partial p}{\partial x}\right)_P \approx \frac{p_e - p_w}{\Delta x} = \frac{\frac{p_P + p_E}{2} - \frac{p_W + p_P}{2}}{\Delta x} = \frac{p_E - p_W}{2\,\Delta x}
$$

The pressure of the cell itself, $p_P$, cancels out. The gradient in cell $P$ only sees the cells two steps away. Now take a pressure field that alternates from cell to cell, $p = +1, -1, +1, -1, \ldots$: every $p_E - p_W$ is zero, so its gradient is **exactly zero everywhere**. The momentum equation does not feel it at all. In 2D the alternating pattern looks like a chessboard, hence the name.

The continuity equation has the same blind spot. If the velocity on a face is the average of the two cells, the mass balance of cell $P$ depends only on $u_E$ and $u_W$, never on $u_P$. The discrete equations split into two (in 3D, eight) families of cells that never talk to each other. Each family gets its own pressure level, and any zig-zag between them is a valid solution.

The classic cure is **Rhie–Chow interpolation** (Rhie and Chow, 1983). The face velocity is not the plain average of the cells, but the average corrected by the difference between the *compact* pressure gradient across the face and the average of the cell gradients:

$$
u_f = \overline{u}_f - \overline{D}_f\left[\frac{p_N - p_P}{\Delta x} - \overline{\left(\frac{\partial p}{\partial x}\right)}_f\right]
$$

where $D = \text{volume}/a_P$ comes from the momentum equation. For a smooth pressure the bracket is tiny, of the order of $\Delta x^2$, and nothing changes. For a zig-zag pressure the compact gradient $(p_N - p_P)/\Delta x$ is large while the average of the cell gradients is zero: the zig-zag now drives mass through the faces, the continuity equation sees it, and the pressure correction removes it. In OpenFOAM this is exactly what the `HbyA`, `rAU` and `phiHbyA` construction in `pEqn.H` does.

<a name="experiment"></a>
## 3. The experiment

[collocated.c](/assets/projects/Checkerboard-Rhie-Chow/collocated.c) is a small collocated SIMPLE solver on a square cavity: first-order upwind convection, Gauss–Seidel for momentum, conjugate gradient for the pressure correction. The only thing that changes is how the face velocity (the mass flux) is built:

1. **linear**: the average of the two cells;
2. **Rhie–Chow**: the average corrected as above;
3. **Rhie–Chow + balanced force**: Rhie–Chow, plus the body force treated on the faces in the same way as the pressure gradient (explained in the next section).

First, the classic lid-driven cavity at Re = 100 on 32 × 32 cells. All three variants converge to practically the same velocity: the minimum of the horizontal velocity on the centreline is −0.1877 with linear interpolation and −0.1878 with Rhie–Chow, against −0.1887 on a staggered grid of the same size. The linear version only shows faint stripes in the pressure near the corners, because the walls couple the families of cells a little. A smooth, lid-driven flow is simply not the case where the problem shows up.

The case where it shows up is a fluid at rest under a force that changes abruptly, which is exactly what happens at the interface between two fluids, or across a sharp temperature jump in a buoyant flow.

<a name="hydrostatic"></a>
## 4. A fluid at rest

The second test removes the moving lid and applies a vertical body force that jumps at mid-height:

$$
f_y = +g \ \text{ for } y < 0.5, \qquad f_y = -g \ \text{ for } y > 0.5, \qquad g = 10
$$

This force can be balanced by pressure alone. The exact solution is a fluid **perfectly at rest**, with a pressure shaped like a tent: rising linearly up to $y = 0.5$, then falling. Any velocity the solver produces is an error, and so is any deviation from the tent.

![Pressure along the vertical centreline, and its deviation from the exact solution](/assets/projects/Checkerboard-Rhie-Chow/pressure_profile.svg)

![Deviation from the exact pressure over the whole cavity](/assets/projects/Checkerboard-Rhie-Chow/perr_maps.png)

*Deviation of the computed pressure from the exact one, 32 × 32 cells. Orange: too high, blue: too low, light grey: exact.*

| Variant | Maximum velocity | Maximum pressure error |
|---|---|---|
| linear interpolation | 2·10⁻⁵ | **0.16**, alternating |
| Rhie–Chow | **0.125** | 0.08, at the jump |
| Rhie–Chow + balanced force | 2·10⁻⁵ | 2·10⁻⁶ |

*32 × 32 cells. Velocities around 10⁻⁵ are the level of the iterative convergence, i.e. zero. The pressure range is about ±2.5.*

Three different outcomes:

- **Linear interpolation keeps the fluid at rest, but the pressure is a checkerboard.** The odd and even rows of cells solve two separate problems, each correct on its own but shifted with respect to the other: the error alternates row by row over the whole domain. The velocity is right only because the momentum equation cannot see the zig-zag in the first place.
- **Rhie–Chow fixes the pressure, but the fluid starts to move.** The pressure is smooth, but a spurious flow appears, up to 0.125 in velocity, concentrated around the jump. The reason is subtle: the cell momentum equations balance the cell-averaged force, while the Rhie–Chow face velocity is driven by the *compact* pressure gradient across each face. At the face where the force jumps, the two no longer agree, and the mismatch pushes mass across the face.

![Spurious velocity with standard Rhie-Chow](/assets/projects/Checkerboard-Rhie-Chow/spurious_map.png)

*Magnitude of the velocity with standard Rhie–Chow in the hydrostatic test: the fluid should be at rest. Darker = faster.*

- **The balanced force fixes both.** If the body force is evaluated on the faces and enters the face velocity exactly like the pressure gradient, and the cell force is reconstructed from the face forces in the same way as the cell pressure gradient, then pressure and force cancel identically at every face. The fluid stays at rest to machine precision, and the pressure is exact.

<a name="refine"></a>
## 5. Refining the mesh

Do these errors disappear with a finer mesh?

| Cells | Linear: pressure error | Rhie–Chow: maximum velocity | Balanced: maximum velocity |
|---|---|---|---|
| 16² | 0.31 | 0.50 | 1·10⁻⁵ |
| 32² | 0.16 | 0.125 | 2·10⁻⁵ |
| 64² | 0.078 | 0.031 | 6·10⁻⁵ |
| 128² | 0.039 | 0.0079 | 4·10⁻⁵ |

*Same hydrostatic test, g = 10.*

- The checkerboard of the linear interpolation halves at each refinement: it is proportional to the cell size.
- The spurious velocity of standard Rhie–Chow drops by 4 at each refinement: it is proportional to the cell size squared. Rhie–Chow is consistent, its error goes away with the mesh.
- The balanced version is exact on every mesh.

So refinement helps, but on a real mesh "small" may not be small enough. With a density ratio of 1000, as between water and air, the force jump at an interface is a thousand times larger than in this test, and surface tension adds an even sharper force concentrated exactly on the interface. Spurious currents of that size can be comparable to the real flow.

<a name="practice"></a>
## 6. What this means in practice

- **Every collocated solver needs Rhie–Chow, or an equivalent.** In OpenFOAM it is built into the pressure equation of every solver (`HbyA`, `phiHbyA`, `rAU`). It is not optional and it is not a setting: without it, the pressure would not be determined.
- **Rhie–Chow has side effects.** Depending on how it is written, the converged solution can depend on the under-relaxation factor or on the time step. Modern implementations, including OpenFOAM's `ddtCorr`, add corrections to remove these dependencies.
- **Body forces must be balanced with the pressure gradient.** This is why OpenFOAM's buoyant and multiphase solvers do not solve for $p$ directly, but for $p_{rgh} = p - \rho\,\mathbf{g}\cdot\mathbf{h}$, and add the gravity term on the faces (`ghf * fvc::snGrad(rho)`) inside the flux, exactly like the "balanced force" variant here. The same idea for surface tension is the *balanced-force* algorithm of Francois et al. (2006).
- **Spurious currents are a symptom.** In a VOF simulation of a droplet or a bubble at rest, any velocity is spurious; its magnitude is a standard check of how well a solver balances pressure, gravity and surface tension.
- **When a pressure field looks noisy cell by cell, think of decoupling first.** A pressure that zig-zags while the velocity looks fine is the classic signature, whether from a missing correction, an inconsistent boundary condition or a force that is not balanced.

The general lesson is the one that runs through this series: a discretization is not just "more or less accurate". It can have blind spots, and the fix is not more resolution but the right consistency between terms that the continuous equations treat as one.

### Reproduce it
All the files are available: [collocated.c](/assets/projects/Checkerboard-Rhie-Chow/collocated.c), [plot_cb.py](/assets/projects/Checkerboard-Rhie-Chow/plot_cb.py) with [svgchart.py](/assets/projects/Checkerboard-Rhie-Chow/svgchart.py), and the data [grid_hydrostatic.csv](/assets/projects/Checkerboard-Rhie-Chow/grid_hydrostatic.csv) and the pressure fields [linear](/assets/projects/Checkerboard-Rhie-Chow/p_linear_N32_a0_g10.csv), [Rhie–Chow](/assets/projects/Checkerboard-Rhie-Chow/p_rhie-chow_N32_a0_g10.csv), [balanced](/assets/projects/Checkerboard-Rhie-Chow/p_balanced_N32_a0_g10.csv).

```bash
gcc -O3 -march=native collocated.c -o collocated -lm
mkdir -p out
./collocated out 32 0 0 1      # lid-driven cavity, Re = 100
./collocated out 32 0 10 0     # fluid at rest with a jump in the body force
mkdir -p figures && python3 plot_cb.py out figures
```
