---
title: "Turbulence in Two Dimensions: Vortices That Grow Instead of Breaking Up"
date: 2026-11-22
tech: ["C", "DNS", "Turbulence", "FFTW"]
excerpt: "In 3D turbulence large eddies break into smaller ones and energy is dissipated at the smallest scales. In 2D the opposite happens: vortices merge and grow, energy is almost conserved and only enstrophy cascades to small scales. A pseudo-spectral DNS of decaying 2D turbulence shows the difference, and why the textbook spectrum is not what you get."
---

## Technical Index
1. [**Why two dimensions are different**](#why) — No vortex stretching
2. [**Two invariants, two cascades**](#cascades) — Energy and enstrophy
3. [**The experiment**](#experiment) — Decaying turbulence, pseudo-spectral DNS
4. [**Vortices that merge**](#vortices) — From noise to a few large vortices
5. [**Energy survives, enstrophy does not**](#budget) — The opposite of 3D
6. [**The spectrum**](#spectrum) — Steeper than the theory
7. [**Where 2D turbulence matters**](#practice) — Atmosphere, oceans, thin layers

<a name="why"></a>
## 1. Why two dimensions are different

In the [hyperviscosity article](/projects/hyperviscosity/) the turbulence was three-dimensional: large eddies break into smaller ones, the energy cascades down to the Kolmogorov scale, and there it is dissipated at a rate that does not depend on the viscosity. The engine of that cascade is **vortex stretching**: a vortex tube pulled along its axis becomes thinner and spins faster, transferring energy to smaller scales.

In two dimensions vortex stretching is impossible: the vorticity is a scalar perpendicular to the plane, and it is simply carried around by the flow. Without the viscosity, the vorticity of each fluid particle is conserved. This single difference changes everything about how 2D turbulence behaves.

<a name="cascades"></a>
## 2. Two invariants, two cascades

Without viscosity, 2D flows conserve two quadratic quantities: the kinetic **energy** $E = \frac{1}{2}\langle|\mathbf{u}|^2\rangle$ and the **enstrophy** $Z = \frac{1}{2}\langle\omega^2\rangle$, the mean squared vorticity. In Fourier space, the enstrophy of a mode is its energy multiplied by $k^2$.

Kraichnan (1967) and Batchelor (1969) showed that conserving both forces them to go in opposite directions:

- the **enstrophy cascades to small scales**, where the viscosity eventually dissipates it, with a predicted spectrum $E(k) \propto k^{-3}$;
- the **energy cascades to large scales** (the *inverse cascade*), where it accumulates in larger and larger structures, with a $k^{-5/3}$ spectrum when the flow is forced at an intermediate scale.

Reproducing the inverse cascade cleanly requires a forced simulation with a large separation between the forcing scale and the box, typically grids of 1024² or more with random forcing; my attempts on 256² did not give a clean $k^{-5/3}$ range, so this article focuses on what a moderate-resolution DNS shows reliably: **decaying** 2D turbulence.

<a name="experiment"></a>
## 3. The experiment

[turb2d.c](/assets/projects/Turbulence-2D/turb2d.c) solves the 2D Navier–Stokes equations in vorticity form,

$$
\frac{\partial\omega}{\partial t} + \mathbf{u}\cdot\nabla\omega = \nu\nabla^2\omega
$$

on a periodic box of side $2\pi$ with a pseudo-spectral method (FFTW), 2/3-rule dealiasing and fourth-order Runge–Kutta with an integrating factor, the same approach as the 3D code of the hyperviscosity article. The initial condition is a random field with its energy concentrated around wavenumber 10, and the flow is left to decay freely. Two runs:

- $\nu = 5\cdot10^{-4}$ on 256² points;
- $\nu = 1.25\cdot10^{-4}$ on 512² points, four times less viscous.

<a name="vortices"></a>
## 4. Vortices that merge

![Vorticity at t = 0, 5, 10 and 30](/assets/projects/Turbulence-2D/vorticity.png)

*Vorticity of the 512² run. Orange: counter-clockwise, blue: clockwise.*

The random initial field organizes itself within a few time units into **coherent vortices**: round, isolated structures that survive for a long time, surrounded by thin filaments of vorticity. Vortices of the same sign that come close **merge** into a larger one, and the filaments thrown off during each merger are stretched and dissipated. By $t = 30$ only a handful of large vortices remain.

The mean wavenumber of the energy spectrum shows the growth of the dominant scale:

| Time | 5 | 10 | 20 | 30 |
|---|---|---|---|---|
| Energy-weighted mean wavenumber | 4.5 | 3.2 | 2.3 | 1.8 |

*512² run.*

The typical size of the energy-carrying structures grows by a factor of 2.5 in 25 time units: energy moves to large scales, the signature of the inverse transfer. In 3D the opposite would happen.

<a name="budget"></a>
## 5. Energy survives, enstrophy does not

![Energy and enstrophy vs time](/assets/projects/Turbulence-2D/energy_enstrophy.svg)

| Run | Energy left at t = 30 | Enstrophy left at t = 30 |
|---|---|---|
| ν = 5·10⁻⁴, 256² | 57% | 2.5% |
| ν = 1.25·10⁻⁴, 512² | 84% | 4.6% |

- **Enstrophy is destroyed quickly** in both runs: more than 95% of it is gone in 30 time units, almost all in the first five, while the vortices form and the filaments are stretched to the dissipation scale.
- **Energy is largely kept**, and the less viscous the flow, the more of it survives: dividing the viscosity by 4 reduces the energy loss from 43% to 16%. As $\nu \to 0$, the energy of decaying 2D turbulence tends to be conserved.

This is the opposite of 3D. In 3D turbulence the energy dissipation rate stays finite when the viscosity goes to zero (the *dissipative anomaly*): the cascade simply makes smaller scales until the viscosity can act. In 2D the energy goes to large scales, where viscosity cannot reach it; only the enstrophy reaches small scales and is dissipated.

<a name="spectrum"></a>
## 6. The spectrum

![Energy spectra at different times](/assets/projects/Turbulence-2D/spectra.svg)

| Time | Slope of E(k), 15 ≤ k ≤ 80 |
|---|---|
| 5 | −4.1 |
| 10 | −4.2 |
| 20 | −4.2 |
| 30 | −3.9 |

*Least-squares fit of the log-log spectrum, 512² run.*

The theory predicts $k^{-3}$ in the enstrophy range; the simulation gives slopes around **−4**, steeper than the theory. This is not a numerical error but a well-known property of *decaying* 2D turbulence (McWilliams, 1984; Benzi et al., 1988). The coherent vortices contain sharp edges and concentrated vorticity, and they do not take part in the cascade: they are islands that the enstrophy flux goes around. The Kraichnan–Batchelor spectrum assumes a self-similar cascade with no such structures; with them, the spectrum steepens, typically to between $k^{-4}$ and $k^{-5}$. The textbook $k^{-3}$ appears more cleanly in forced simulations, where the forcing continuously regenerates small-scale vorticity.

It is a useful reminder that turbulence theories describe idealized states. A simulation that "does not match the $-3$" is not necessarily wrong; it may be showing physics the theory leaves out.

<a name="practice"></a>
## 7. Where 2D turbulence matters

Real flows are three-dimensional, but many are close to 2D at large scales:

- **The atmosphere and the oceans** are thin compared with their horizontal extent and strongly influenced by rotation and stratification, which suppress vertical motion. Large-scale weather systems and ocean eddies behave much more like 2D turbulence than like 3D: energy moves to large scales, which is one reason why cyclones and jets are so persistent.
- **Thin layers and soap films**, flows between closely spaced plates, and strongly magnetized plasmas are used in the laboratory to study 2D turbulence directly.
- **For CFD**, the main lesson is a warning: a 2D simulation of a flow that is turbulent in reality (a 2D "LES" of a wake, for example) does not have the right physics. Without vortex stretching, energy does not cascade to small scales, vortices become too strong and too long-lived, and the results can be qualitatively wrong. 2D is fine for laminar and RANS problems; it is not a cheap substitute for 3D turbulence.
- **Numerically**, 2D turbulence is a good test for schemes: an energy- and enstrophy-conserving discretization (such as the Arakawa Jacobian used in the [reproducibility article](/projects/reproducibility/)) is needed to avoid spurious accumulation or loss at the grid scale.

### Reproduce it
All the files are available: [turb2d.c](/assets/projects/Turbulence-2D/turb2d.c) (needs FFTW 3), [plot_2d.py](/assets/projects/Turbulence-2D/plot_2d.py) with [svgchart.py](/assets/projects/Turbulence-2D/svgchart.py), and the data [history_N256_p1.csv](/assets/projects/Turbulence-2D/history_N256_p1.csv), [history_N512_p1.csv](/assets/projects/Turbulence-2D/history_N512_p1.csv), [spectra_N512_p1.csv](/assets/projects/Turbulence-2D/spectra_N512_p1.csv).

```bash
gcc -O3 -march=native turb2d.c -o turb2d -lfftw3 -lm
mkdir -p snap out
./turb2d snap 512 1 1.25e-4 30      # ~9 min, with vorticity snapshots
./turb2d out 256 1 5e-4 30          # ~1 min
mkdir -p figures && python3 plot_2d.py snap out figures
```
