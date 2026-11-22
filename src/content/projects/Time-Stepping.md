---
title: "Stable Is Not Accurate: Time Steps and Crank-Nicolson"
date: 2026-11-22
tech: ["C", "CFD", "Numerical Methods", "OpenFOAM"]
excerpt: "Implicit schemes let you take huge time steps without blowing up. But a stable answer can still lose 80% of a wave in one period, and Crank-Nicolson can turn a step upside down. Euler, backward, Crank-Nicolson and RK4 measured with the spatial error removed."
---

## Technical Index
1. [**Stability and accuracy**](#question) — Two different questions
2. [**The experiment**](#experiment) — Time error only
3. [**Carrying a wave**](#advection) — Damping and phase lag vs Courant number
4. [**Diffusing a step**](#diffusion) — Why Crank-Nicolson oscillates
5. [**What this means in practice**](#practice) — `ddtSchemes`, Courant numbers and `CrankNicolson 0.9`

<a name="question"></a>
## 1. Stability and accuracy

In a transient simulation the time step is usually chosen by one number, the **Courant number** $C = u\,\Delta t / \Delta x$: how many cells the flow crosses in one step. Explicit schemes become unstable above a limit of order 1, so with an explicit scheme the question answers itself. Implicit schemes remove that limit: OpenFOAM's `Euler`, `backward` and `CrankNicolson` will happily run at Courant 5, 10 or 50 without blowing up.

But stability only says that the solution stays bounded. It says nothing about whether it is *right*. This article measures how wrong a stable solution can be, with the schemes available in OpenFOAM's `ddtSchemes`.

<a name="experiment"></a>
## 2. The experiment

[timestep.c](/assets/projects/Time-Stepping/timestep.c) solves two model equations on a periodic domain with 256 points. To measure the time error alone, space is treated exactly: the solution is split into Fourier modes and each mode is advanced in time on its own, so the spatial derivative is computed without error and every difference from the exact solution comes from the time scheme.

The schemes:

- **Euler explicit** and **RK4**: explicit, first and fourth order.
- **Euler implicit** (`Euler` in OpenFOAM): implicit, first order.
- **BDF2** (`backward` in OpenFOAM): implicit, second order, uses the two previous time levels.
- **Crank–Nicolson** ($\theta = 0.5$, `CrankNicolson 1`): implicit, second order, the average of explicit and implicit Euler.
- **θ = 0.55**: Crank–Nicolson pushed slightly towards implicit Euler. OpenFOAM's `CrankNicolson ψ` blends the two in a similar way, and ψ = 0.9 is the common choice; this θ plays the same role.

The whole program runs in 0.15 seconds.

<a name="advection"></a>
## 3. Carrying a wave

The first test transports a Gaussian bump, $u_t + a\,u_x = 0$, for exactly one period: after the last step it should be back where it started, with the same height.

![Peak height and phase lag after one period vs Courant number](/assets/projects/Time-Stepping/advection.svg)

| Courant | RK4 | Euler implicit | BDF2 | Crank–Nicolson | θ = 0.55 |
|---|---|---|---|---|---|
| 0.1 | 1.000 | 0.897 | 1.000 | 1.000 | 0.988 |
| 0.5 | 1.000 | 0.671 | 0.994 | 1.000 | 0.944 |
| 1 | unstable | 0.540 | 0.969 | 0.998 | 0.896 |
| 2 | unstable | 0.413 | 0.868 | 0.986 | 0.817 |
| 5 | unstable | 0.277 | 0.612 | 0.860 | 0.653 |
| 10 | unstable | 0.203 | 0.425 | 0.675 | 0.508 |

*Height of the peak after one period (exact: 1). Explicit Euler is unstable at every Courant number for pure advection, so it is not in the table.*

Three observations:

- **Explicit schemes are not an option here.** Explicit Euler amplifies every wave of a pure advection problem, at any time step. RK4 is very accurate but stops being stable around Courant 0.9. (In a real finite-volume code with upwind convection, explicit Euler becomes stable for $C \le 1$, because the [numerical diffusion](/projects/numerical-diffusion/) of upwind damps the amplification.)
- **Implicit Euler is stable and very diffusive.** Even at Courant 0.1, it loses 10% of the bump in one period, and at Courant 1 it loses half of it. Its error is the time-domain twin of upwind's numerical diffusion: a first-order scheme that damps everything.
- **Second-order schemes are excellent at small Courant numbers and degrade quickly above 1.** Crank–Nicolson is the best of all: it does not damp a single wave mode at all. What it loses at large time steps is phase: the waves travel at the wrong speed and the bump spreads out (dispersion). BDF2 both damps and delays.

The profiles at Courant 2 show the different kinds of error:

![The Gaussian after one period at Courant number 2](/assets/projects/Time-Stepping/advection_profiles.svg)

Implicit Euler has flattened the bump to 41% of its height. BDF2 lags behind and leaves a negative dip in front of the bump. Crank–Nicolson keeps almost all the height but slightly behind its exact position. None of them blew up, and at Courant 2 none of them is right.

<a name="diffusion"></a>
## 4. Diffusing a step

The second test is diffusion, $u_t = \nu\,u_{xx}$, of a step, with the usual three-point Laplacian of a finite-volume code and a large time step: $\nu\,\Delta t/\Delta x^2 = 20$, a "diffusive Courant number" far above the explicit limit of 0.5. Physically the step should only get smoother: the exact solution is monotone and its total variation (the sum of the jumps between neighbouring points) stays at 2.

![The edge of the step after 1 and after 20 time steps](/assets/projects/Time-Stepping/diffusion.svg)

| After | Euler implicit | BDF2 | Crank–Nicolson | θ = 0.55 |
|---|---|---|---|---|
| 1 step | 0.06 / 2.00 | 0.06 / 2.00 | **0.38 / 4.75** | 0.31 / 4.19 |
| 2 steps | 0.03 / 2.00 | 0.02 / 2.00 | 0.33 / 3.71 | 0.22 / 2.90 |
| 5 steps | 0.01 / 2.00 | 0.003 / 2.00 | 0.24 / 3.80 | 0.09 / 2.58 |
| 20 steps | 0.003 / 1.90 | 0.0001 / 1.91 | 0.08 / 2.60 | 0.001 / 1.90 |

*Maximum error / total variation. Exact total variation: 2.00 (1.905 after 20 steps). A value above the exact one means the profile is no longer monotone.*

- **Implicit Euler and BDF2 stay monotone.** The step is smoothed, a bit too much by Euler, almost exactly by BDF2.
- **Crank–Nicolson turns the edge of the step upside down.** Next to the jump, the points just below the edge jump up to 0.84 and the points just above drop to 0.16: the profile zig-zags, and the total variation more than doubles. After 20 steps the oscillation is still there.

The reason is in how each scheme treats the shortest waves on the grid. Diffusion should kill them almost instantly. Implicit Euler multiplies them by $1/(1 + z)$, which goes to zero when the time step is large. Crank–Nicolson multiplies them by $(1 - z/2)/(1 + z/2)$, which goes to **−1**: they are not damped at all, they just flip sign at every step. The same property that made Crank–Nicolson perfect for carrying a wave makes it unable to get rid of grid-scale noise. In the jargon, it is A-stable but not L-stable.

A small step towards implicit Euler is enough to fix this: with θ = 0.55 the shortest waves are multiplied by about −0.8 per step instead of −1, so the oscillation is still created but dies out within a few steps, and after 20 steps the solution is as accurate as BDF2.

<a name="practice"></a>
## 5. What this means in practice

- **An implicit scheme removes the stability limit, not the accuracy limit.** For transient results that matter (vortex shedding, wakes, frequencies), the Courant number has to stay small, typically below 1 in the regions of interest, whatever the scheme.
- **`Euler` is for robustness, not accuracy.** It is the time-domain equivalent of first-order upwind. It is fine for pseudo-transient runs towards a steady state, where the time history does not matter, and wrong for anything where it does.
- **`backward` (BDF2) is a solid default for accurate transients:** second order, monotone on diffusion, but it damps and delays waves above Courant 1.
- **`CrankNicolson 1` is the most accurate for waves and the least robust:** it never damps the shortest modes, so any grid-scale noise, from a sharp initial condition, a poor mesh region or a sudden change in boundary conditions, survives and can grow through the non-linear terms. That is why `CrankNicolson 0.9` (slightly off-centred towards Euler) is the usual recommendation: almost all the accuracy, with some damping of the noise.
- **Stable is not the same as converged in time.** As with mesh refinement, the only way to know is to halve the time step and check that the result does not change.

As in the rest of the series: a scheme that runs without errors is not necessarily solving the problem you think. Here it is the time step that can quietly turn a wave into a smear or a step into a zig-zag.

### Reproduce it
All the files are available: [timestep.c](/assets/projects/Time-Stepping/timestep.c), [plot_ts.py](/assets/projects/Time-Stepping/plot_ts.py) with [svgchart.py](/assets/projects/Time-Stepping/svgchart.py), and the data [advection.csv](/assets/projects/Time-Stepping/advection.csv), [diffusion.csv](/assets/projects/Time-Stepping/diffusion.csv).

```bash
gcc -O3 -march=native timestep.c -o timestep -lm
mkdir -p out && ./timestep out              # 0.15 s
mkdir -p figures && python3 plot_ts.py out figures
```
