---
title: "When the Scheme Becomes the Viscosity: the Taylor-Green Vortex"
date: 2026-11-22
tech: ["C", "CFD", "Turbulence", "LES"]
excerpt: "The Taylor-Green vortex at Re = 1600 solved with central, linearUpwind and upwind convection. Upwind gets the peak dissipation almost right, and 89% of it comes from the scheme, not from the viscosity. Measuring how much of a turbulent flow's dissipation is physics and how much is numerics."
---

## Technical Index
1. [**The test case**](#case) — A vortex that becomes turbulent and decays
2. [**The experiment**](#experiment) — One solver, three convection schemes
3. [**No viscosity at all**](#inviscid) — Which schemes conserve energy
4. [**Re = 1600**](#results) — Physical and numerical dissipation
5. [**The right answer for the wrong reason**](#wrong) — Why upwind looks good
6. [**What this means in practice**](#practice) — LES, schemes and resolution

<a name="case"></a>
## 1. The test case

The **Taylor–Green vortex** is the standard benchmark for schemes used in turbulence simulations. It starts from a smooth, periodic arrangement of vortices,

$$
u = \sin x \cos y \cos z, \qquad v = -\cos x \sin y \cos z, \qquad w = 0
$$

in a box of side $2\pi$. The vortices stretch and break up, the flow becomes turbulent around $t \approx 8$–$9$, and then it decays. At a Reynolds number of 1600 the reference solution is known from high-resolution spectral DNS: the dissipation rate $-dE/dt$ peaks at about **0.0127** around $t \approx 9$ (Brachet et al., 1983; van Rees et al., 2011).

The question of this article: when a coarser simulation reproduces that peak, how much of its dissipation comes from the viscosity of the fluid, and how much from the numerical scheme?

<a name="experiment"></a>
## 2. The experiment

[tgv.c](/assets/projects/Taylor-Green/tgv.c) solves the 3D incompressible Navier–Stokes equations with the kind of discretization a finite-volume code uses on a Cartesian mesh: staggered grid, second-order viscous term, exact pressure projection, third-order Runge–Kutta in time. Only the face values of the convective term change:

- **central**: the average of the two neighbouring values, second order, the scheme used for LES;
- **linearUpwind**: second-order upwind, a common default in RANS (`linearUpwind grad(U)` in OpenFOAM);
- **upwind**: first order.

Meshes of 32³, 64³ and 128³ cells, up to $t = 20$.

The key diagnostic is that kinetic energy can be removed by only two things. The **physical dissipation** is $\nu\langle|\nabla\mathbf{u}|^2\rangle$, computed from the velocity gradients. The **total dissipation** is the actual rate of energy loss, $-dE/dt$. With an exact scheme the two would coincide; their difference is the **numerical dissipation** of the scheme.

<a name="inviscid"></a>
## 3. No viscosity at all

First the simplest check: remove the viscosity ($Re = \infty$). The exact equations then conserve kinetic energy exactly, and so should a good scheme.

| Scheme (32³, from t = 0 to 2) | Energy at t = 0 | Energy at t = 2 | Lost |
|---|---|---|---|
| central | 0.12500 | 0.12500 | 0.00006% |
| linearUpwind | 0.12500 | 0.12329 | 1.4% |
| upwind | 0.12500 | 0.09943 | **20%** |

The central scheme conserves energy to within the time-integration error: on a staggered grid, its convective term transfers energy between scales but never creates or destroys it. The upwind schemes destroy energy even when there is no viscosity, and even when the flow is still smooth: upwind loses a fifth of the energy in two time units.

<a name="results"></a>
## 4. Re = 1600

![Total dissipation vs time on 64³ and 128³](/assets/projects/Taylor-Green/dissipation_total.svg)

| Mesh | Scheme | Peak of −dE/dt | at t | Share of the dissipation that is numerical |
|---|---|---|---|---|
| 64³ | central | 0.0119 | 8.3 | 0% |
| 64³ | linearUpwind | 0.0118 | 8.1 | 70% |
| 64³ | upwind | 0.0126 | 5.0 | **89%** |
| 128³ | central | 0.0136 | 8.6 | 0% |
| 128³ | linearUpwind | 0.0121 | 8.2 | 46% |
| 128³ | upwind | 0.0125 | 5.3 | **83%** |

*Reference DNS: peak 0.0127 at t ≈ 9. Share: integral of (total − physical) over the integral of the total dissipation, 0 ≤ t ≤ 20 (32³: 0%, 84% and 94%).*

![Physical and total dissipation on 64³](/assets/projects/Taylor-Green/physical_vs_numerical.svg)

*For the central scheme the total and physical curves coincide, so only one is visible.*

- **Central dissipates only through viscosity.** The physical and total curves coincide. On 64³ the peak is a little low and early (0.0119 at 8.3): the mesh cannot resolve the smallest scales, so the velocity gradients, and with them the viscous dissipation, are underestimated. On 128³ it slightly overshoots (0.0136 at 8.6), typical of a non-dissipative scheme that is still not fully resolving the flow. With refinement it converges to the reference.
- **linearUpwind gets the timing right, but most of its dissipation is numerical**: 70% on 64³, still 46% on 128³. The share drops with refinement, as it should for a second-order scheme, but slowly.
- **upwind is dominated by numerics**: 83–94% of the energy loss comes from the scheme at every resolution tested, and the dissipation starts at $t = 0$, when the flow is still a smooth, laminar vortex. The physical dissipation (the velocity gradients) is ten times too small, because the scheme smears the small scales before they can form.

![Kinetic energy vs time](/assets/projects/Taylor-Green/energy.svg)

<a name="wrong"></a>
## 5. The right answer for the wrong reason

The most instructive number in the table is the peak dissipation of upwind on 64³: **0.0126**, within 1% of the reference. Looked at alone, it would suggest an excellent result. It is not. It comes four time units too early, 89% of it is produced by the scheme, and the turbulent small scales the real flow develops are not there at all.

This is the general trap of dissipative schemes in turbulent flows: an under-resolved turbulent flow needs to dissipate energy somewhere, and numerical dissipation provides a sink. The total energy budget can then look plausible while the physics that should produce it, the cascade of energy to small eddies, is replaced by the scheme. The upwind solution behaves like the flow of a much more viscous fluid, the same conclusion as the [numerical diffusion article](/projects/numerical-diffusion/), now in a turbulent flow.

There is an approach that exploits this on purpose: **implicit LES** uses the dissipation of a carefully designed scheme as the subgrid model. The difference is that the scheme is chosen and tuned so that its dissipation acts mainly at the smallest resolved scales, much like [hyperviscosity](/projects/hyperviscosity/); plain upwind acts on every scale.

<a name="practice"></a>
## 6. What this means in practice

- **LES and DNS need non-dissipative (or very weakly dissipative) convection schemes.** In OpenFOAM, `linear` for the velocity in LES, often blended (`LUST`, or `filteredLinear`) for robustness; never `upwind`, and `linearUpwind` only with care. Otherwise the subgrid model competes with a numerical dissipation that can be larger than its own.
- **Check the energy budget.** Comparing the physical dissipation with the actual energy loss is a direct way to see how much of a simulation is controlled by the scheme. A large gap means the results reflect the discretization more than the model.
- **A matching integral quantity is not validation.** Upwind matched the peak dissipation within 1% for the wrong reason. Look at when it happens, at spectra and at the flow structures, not only at one number.
- **Central schemes need resolution and good meshes.** They conserve energy, so they do not hide under-resolution; they show it as noise or overshoot (as on 128³ here) and can become unstable on poor meshes. Smooth meshes (see the [expansion-ratio article](/projects/expansion-ratio/)) and enough cells are the price.
- **For RANS the situation is different.** In a RANS simulation the turbulence model provides a large eddy viscosity, and second-order upwind schemes are generally adequate; the trap described here concerns simulations that try to resolve turbulence.

### Reproduce it
All the files are available: [tgv.c](/assets/projects/Taylor-Green/tgv.c) (needs FFTW 3), [plot_tgv.py](/assets/projects/Taylor-Green/plot_tgv.py) with [svgchart.py](/assets/projects/Taylor-Green/svgchart.py), and the energy histories `tgv_N*_Re1600.csv`.

```bash
gcc -O3 -march=native tgv.c -o tgv -lfftw3 -lm
mkdir -p out
./tgv out 64 0 1600 20     # central, 64^3, ~1 min (1 = upwind, 2 = linearUpwind)
./tgv out 128 0 1600 20    # central, 128^3, ~12 min
./tgv out 32 1 0 2         # inviscid check
mkdir -p figures && python3 plot_tgv.py out figures
```
