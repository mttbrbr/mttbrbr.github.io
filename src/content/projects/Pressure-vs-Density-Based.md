---
title: "Pressure-Based vs Density-Based at Low Mach"
date: 2026-11-22
tech: ["C", "CFD", "Numerical Methods", "Compressible Flow"]
excerpt: "Why a density-based solver cannot simply be used for incompressible flow: stiffness, dissipation that scales with the speed of sound and round-off, measured on the 1D Euler equations, and what the all-speed (preconditioning) treatment changes."
---

## Technical Index
1. [**Two ways to solve the same equations**](#families) — What each family takes as primary unknown
2. [**The incompressible limit**](#limit) — Where density stops carrying information
3. [**Problem 1: stiffness**](#stiffness) — Following sound waves you do not care about
4. [**Problem 2: dissipation scaled by the speed of sound**](#dissipation) — Numerical diffusion, 1/M times larger
5. [**Problem 3: round-off**](#roundoff) — Pressure differences lost in the last digits
6. [**The all-speed treatment**](#precond) — Preconditioning and dual time stepping
7. [**What this means in practice**](#practice) — Choosing the solver

Every general-purpose CFD code offers two families of solvers: **pressure-based** and **density-based**. Both solve the same Navier–Stokes equations, and both can handle compressible flow. Yet there is a rule every CFD user learns early: a density-based solver should not be used for incompressible or very low-speed flow, unless it has a special *low-Mach* or *all-speed* treatment. In this article I want to show *why*, with numbers, using the simplest possible model: the 1D Euler equations and a few hundred lines of C.

<a name="families"></a>
## 1. Two ways to solve the same equations

Both families start from conservation of mass, momentum and energy. They differ in which variable they treat as the main unknown and how they close the system.

| | Density-based | Pressure-based |
|---|---|---|
| Unknowns | $\rho$, $\rho\mathbf{u}$, $\rho E$ | $\mathbf{u}$, $p$ (and $T$) |
| Density | continuity eq. | equation of state |
| Pressure | equation of state | pressure equation |
| Coupling | all together | SIMPLE / PISO / PIMPLE |
| Born for | shocks, high speed | low speed |

*Pressure equation: built from continuity and momentum. Examples of density-based solvers: Fluent density-based, OpenFOAM `rhoCentralFoam`, SU2. Pressure-based: Fluent pressure-based, OpenFOAM `simpleFoam`, `pimpleFoam`, `rhoPimpleFoam`.*

In a density-based solver the pressure is a *consequence*: compute the new density and energy, then read the pressure off the equation of state. Information travels with the physical waves of the system, the flow speed $u$ and the acoustic waves $u \pm c$, which is exactly what makes these solvers excellent for shocks.

<a name="limit"></a>
## 2. The incompressible limit

Now let the flow get slow. The Mach number $M = |u|/c$ goes to zero, and two things happen.

First, density stops changing. For an isentropic flow the relative density change is about $M^2/2$: at $M = 0.1$ it is 0.5%, at $M = 0.01$ it is 0.005%. In the limit, the continuity equation $\partial\rho/\partial t + \nabla\cdot(\rho\mathbf{u}) = 0$ becomes the *constraint* $\nabla\cdot\mathbf{u} = 0$. It no longer tells you how density evolves, because density does not evolve.

Second, pressure splits into two very different parts:

$$
p(\mathbf{x}, t) = p_0(t) + M^2\, p_2(\mathbf{x}, t)
$$

a large thermodynamic pressure $p_0$, uniform in space, and a tiny dynamic part of order $\rho u^2$, which is the only part that pushes the flow. In a truly incompressible flow the speed of sound is infinite and $p_2$ is no longer linked to density at all: it becomes whatever it must be to keep $\nabla\cdot\mathbf{u} = 0$ (mathematically, a Lagrange multiplier).

This is why a density-based solver is in trouble: it wants to compute the pressure from the density, but in the limit the density carries **no information** about the pressure. A pressure-based solver does not have this problem: it builds an equation for the pressure directly from the divergence constraint, which is exactly the equation that survives in the limit.

Between "compressible" and "incompressible" there is no sharp wall, though. At low but finite Mach number a density-based solver still works in principle. The question is at what cost and with what accuracy, and that is what the next three sections measure.

### The test
[lowmach.c](/assets/projects/Low-Mach-Density-Based/lowmach.c) solves the 1D Euler equations on a periodic domain of length 1 with 200 cells. The flow is uniform ($u = M$, pressure constant, speed of sound $c = 1$) and carries a small density bump (+10%) once across the domain. Physically nothing happens except transport: the bump should come back unchanged.

- **Density-based:** explicit finite volumes on $(\rho, \rho u, \rho E)$ with the **Rusanov** (local Lax–Friedrichs) flux, time step from the fastest wave, CFL 0.9. Rusanov is the simplest upwind-type flux, and its dissipation, like that of the Kurganov–Tadmor central schemes used by OpenFOAM's `rhoCentralFoam`, is set by the largest local wave speed $|u| + c$.
- **Incompressible reference:** what any incompressible (pressure-based) solver reduces to for this flow: the velocity is uniform, the pressure is trivially constant, and the bump is simply transported at $u$ with first-order upwind, time step from $u$, CFL 0.9.

The two are deliberately comparable: same grid, same CFL number, both first order.

<a name="stiffness"></a>
## 3. Problem 1: stiffness

The 1D Euler equations have three waves, with speeds $u$, $u + c$ and $u - c$. An explicit density-based solver must resolve the fastest one, so its time step is $\Delta t \le \text{CFL}\,\Delta x / (|u| + c)$. But the flow itself, and everything you are interested in, moves at $u$. The number of steps to carry the flow across the domain is therefore proportional to $(|u| + c)/|u| = 1 + 1/M$.

![Time steps to carry the flow once across the domain, vs Mach number](/assets/projects/Low-Mach-Density-Based/steps.svg)

| Mach | DB steps | DB time (ms) | Inc. steps | Inc. time (ms) |
|---|---|---|---|---|
| 0.5 | 667 | 2.7 | 223 | 0.05 |
| 0.1 | 2 443 | 10.1 | 223 | 0.04 |
| 0.01 | 22 329 | 90.8 | 223 | 0.05 |
| 0.001 | 221 088 | 917 | 223 | 0.05 |

*DB: density-based, Inc.: incompressible reference. Steps: explicit time steps needed to carry the flow once across the domain. Time: wall-clock time in milliseconds on one core (i7-1165G7).*

At $M = 0.001$ (air at 0.34 m/s, a slow breeze in a room) the density-based solver needs **about 1000 times more time steps** to simulate the same physical time. It spends almost all its effort following sound waves that the physics of the problem does not care about. Implicit time integration relaxes the stability limit, but the same disparity shows up as an ill-conditioned linear system: the ratio between fastest and slowest wave speed, $1 + 1/M$, is the condition number that slows down every iterative method.

<a name="dissipation"></a>
## 4. Problem 2: dissipation scaled by the speed of sound

Stiffness only costs time. The second problem costs **accuracy**, and connects directly to the [previous article on numerical diffusion](/projects/numerical-diffusion/).

An upwind-type flux stabilizes the scheme by adding dissipation proportional to the wave speeds it upwinds. The Rusanov flux uses the largest one, $|u| + c$, for every variable. Following the same reasoning as for the upwind scheme, its modified equation for the transported density gets a numerical viscosity

$$
\nu_{num} \approx \frac{\Delta x}{2}\,(|u| + c)
$$

while an incompressible upwind scheme only adds $\Delta x\,|u|/2$. At low Mach the speed of sound dominates, and the ratio between the two is $1 + 1/M$. I measured it exactly as before, from how fast the variance of the density bump grows:

![Ratio of numerical viscosities vs Mach number](/assets/projects/Low-Mach-Density-Based/viscosity_ratio.svg)

| Mach | DB measured | DB predicted | Inc. upwind | Ratio |
|---|---|---|---|---|
| 0.5 | 3.34e-3 | 3.38e-3 | 1.25e-3 | 2.7× |
| 0.1 | 2.69e-3 | 2.73e-3 | 2.50e-4 | 11× |
| 0.01 | 2.49e-3 | 2.52e-3 | 2.50e-5 | 100× |
| 0.001 | 2.47e-3 | 2.50e-3 | 2.50e-6 | 986× |

*Numerical viscosity, dimensionless (domain length 1, speed of sound 1). DB measured: from the spreading of the bump with the density-based scheme. DB predicted: $(|u| + c)\,\Delta x/2$ with the small time-integration correction; they agree within 1.5%. Inc. upwind: $|u|\,\Delta x/2$. Ratio: DB measured / Inc. upwind.*

The key observation: the numerical viscosity of the density-based scheme **does not go down with the flow speed**. It stays at about $c\,\Delta x/2$ whatever the Mach number, while the physics it should transport gets slower and slower. The time needed to cross the domain grows like $1/M$, so the total smearing grows like $1/M$ too:

![Density bump after one crossing of the domain at Mach 0.5, 0.1 and 0.01](/assets/projects/Low-Mach-Density-Based/profiles.svg)

| Mach | DB peak left | Inc. peak left |
|---|---|---|
| 0.5 | 40% | 91% |
| 0.1 | 21% | 91% |
| 0.02 | 13% | 91% |
| 0.01 | 12.5% | 91% |

*Height of the density bump after it crossed the domain once, as a percentage of its initial height. 12.5% means completely flat: the mass of the bump spread uniformly over the whole domain.*

From $M = 0.02$ down, the density-based solution after one crossing is a **flat line**: the bump has been completely diffused away. The incompressible solver, with the same grid and the same order, keeps 91% of it at every Mach number.

A note on fairness. Rusanov is the most dissipative choice; a Roe or HLLC flux upwinds each wave separately and transports this particular density bump (a contact wave) with a dissipation proportional to $|u|$ only. But the problem does not disappear, it moves: in two and three dimensions, the dissipation on the acoustic waves still scales with $c$ and acts on velocity and pressure. Guillard and Viozat showed in 1999 that, as a result, upwind schemes produce pressure fluctuations of order $M$ instead of the physical $M^2$: at low Mach the computed pressure field is simply wrong, unless the dissipation is rescaled. That rescaling is part of the all-speed treatment of section 6.

<a name="roundoff"></a>
## 5. Problem 3: round-off

The third problem is more subtle and has nothing to do with the scheme. A density-based solver stores the total energy $\rho E = p/(\gamma - 1) + \rho u^2/2$ and recovers the pressure by subtraction:

$$
p = (\gamma - 1)\left(\rho E - \frac{(\rho u)^2}{2\rho}\right)
$$

The part of the pressure that drives a low-Mach flow is $p_2 \sim M^2 p_0$, a tiny difference on top of a large number. With 7 significant digits (single precision) or 16 (double), anything smaller than about $10^{-7}$ or $10^{-16}$ of $p_0$ is lost. The program builds the exact conserved variables for a pressure fluctuation $\delta p = M^2 p_0$, rounds them to float and double, and recovers $\delta p$ again:

![Relative error on the recovered pressure fluctuation vs Mach number](/assets/projects/Low-Mach-Density-Based/roundoff.svg)

| Mach | Error, single | Error, double |
|---|---|---|
| 1e-2 | 0.08% | 2e-11 % |
| 1e-3 | 2.5% | 3e-9 % |
| 1e-4 | 138% | 7e-7 % |
| 1e-6 | 2.4e6 % | 0.004% |
| 1e-8 | 2.4e10 % | 82% |

*Relative error on the recovered pressure fluctuation $\delta p$, in percent, in single (float) and double precision. Best case: exact values rounded once, no accumulation over iterations.*

In single precision the dynamic pressure is already wrong by 2.5% at $M = 10^{-3}$ and pure noise at $M = 10^{-4}$. Double precision pushes the limit far away, which is one reason why low-speed compressible simulations are normally run in double precision. The other, more elegant fix is to never store the absolute pressure: pressure-based solvers work with a **gauge pressure** relative to a reference. This is what Fluent's *operating pressure* is for, and why OpenFOAM's incompressible solvers use a kinematic pressure defined up to a constant.

<a name="precond"></a>
## 6. The all-speed treatment

So can a density-based solver handle incompressible flow at all? Yes, with **low-Mach preconditioning**, the "all-speed" treatment introduced by Turkel, Choi and Merkle, and Weiss and Smith in the 1980s and 90s. The idea is to multiply the time derivative by a matrix $\Gamma$:

$$
\Gamma\,\frac{\partial Q}{\partial \tau} + \frac{\partial F(Q)}{\partial x} = 0
$$

$\Gamma$ is chosen so that the acoustic waves of the modified system travel at a speed comparable to $u$ instead of $c$. For a steady solution the time derivative vanishes, so the result is the same steady state, but the system that gets there is no longer stiff. In 1D, with Turkel's preconditioner and parameter $\varepsilon = M^2$, the acoustic wave speeds become

$$
\lambda_\pm = \frac{1}{2}\left[(1 + \varepsilon)\,u \pm \sqrt{(1 - \varepsilon)^2 u^2 + 4\,\varepsilon\,c^2}\right]
$$

which for small $M$ are about $1.62\,u$ and $-0.62\,u$: of the same order as the flow speed.

![Ratio of fastest to slowest wave speed, with and without preconditioning](/assets/projects/Low-Mach-Density-Based/eigen.svg)

| Mach | Plain | Preconditioned |
|---|---|---|
| 0.5 | 3.0 | 3.8 |
| 0.1 | 11 | 2.65 |
| 0.01 | 101 | 2.62 |
| 0.001 | 1001 | 2.62 |

*Ratio between the largest and smallest absolute wave speed, without and with preconditioning: the stiffness of the system, and the condition number seen by an iterative solver.*

With preconditioning the stiffness stops depending on the Mach number: **2.62 at every low Mach**, instead of 1001 at $M = 0.001$. At $M = 0.5$ it is slightly worse than no preconditioning: the treatment is designed for low speed, and real codes switch it off gradually as $M$ approaches 1 (typically by letting $\varepsilon$ go to 1).

The full all-speed treatment has three parts:

1. **Preconditioning the time derivative**, as above: it fixes the stiffness.
2. **Preconditioning the dissipation**: the numerical flux uses the wave speeds of the preconditioned system instead of $|u| \pm c$. This fixes the accuracy problem of section 4: the dissipation now scales with $u$, and the pressure fluctuations recover the correct $M^2$ scaling.
3. **Dual time stepping** for unsteady flows: since $\Gamma$ changes the time derivative, it is applied only to a *pseudo-time* $\tau$. Each physical time step becomes an inner steady problem in $\tau$, solved with preconditioning, and the physical time derivative is kept unmodified. More expensive per step, but accurate in time.

It is not free: the preconditioner needs a cutoff to avoid dividing by zero at stagnation points (where $u \to 0$), and the parameters affect robustness. But with it, a density-based solver can work from almost incompressible flow up to hypersonic speeds. That is what *all-speed* means. Round-off (section 5) still has to be handled separately, with double precision and a gauge pressure.

<a name="practice"></a>
## 7. What this means in practice

- **Incompressible or low-speed flow: use a pressure-based solver.** As a rule of thumb, below $M \approx 0.3$ density changes are under about 5% and the incompressible assumption is common. In OpenFOAM that means `simpleFoam`, `pimpleFoam` and friends, or `rhoSimpleFoam` / `rhoPimpleFoam` when density changes with temperature. In Fluent, the pressure-based solver.
- **Do not use a density-based solver without low-Mach treatment at low speed.** OpenFOAM's `rhoCentralFoam` is a density-based solver with central-upwind fluxes and no preconditioning: its dissipation scales with $|u| + c$ exactly like the Rusanov flux above, so at low Mach it is slow and diffusive. It is built for high-speed flows with shocks. Commercial density-based solvers apply low-speed preconditioning (Fluent uses the Weiss–Smith formulation) precisely to make them usable at low speed.
- **High speed with strong shocks: density-based is the natural choice.** Pressure-based solvers can handle transonic and supersonic flow too (OpenFOAM's `transonic` option makes the pressure equation account for compressibility), but density-based solvers with upwind or central-upwind fluxes are usually more robust and sharper on strong shocks.
- **Mixed regimes**, like a low-speed flow that accelerates through a nozzle, or a combustion chamber with large density changes but low velocity, are exactly where the choice matters: either a pressure-based solver with full compressibility, or a density-based solver with an all-speed treatment.
- **Always run low-Mach compressible cases in double precision**, and with a sensible reference (operating) pressure.

The pattern is the same as in the article on numerical diffusion: a scheme is not just "more or less accurate", it carries the physics it was designed for. A density-based scheme is built around sound waves; when the sound waves stop mattering, the scheme does not know it unless you tell it.

### Reproduce it
All the files are available: [lowmach.c](/assets/projects/Low-Mach-Density-Based/lowmach.c), [plot_lm.py](/assets/projects/Low-Mach-Density-Based/plot_lm.py) and the data [crossing.csv](/assets/projects/Low-Mach-Density-Based/crossing.csv), [viscosity.csv](/assets/projects/Low-Mach-Density-Based/viscosity.csv), [roundoff.csv](/assets/projects/Low-Mach-Density-Based/roundoff.csv), [eigen.csv](/assets/projects/Low-Mach-Density-Based/eigen.csv), [profiles.csv](/assets/projects/Low-Mach-Density-Based/profiles.csv).

```bash
gcc -O3 -march=native lowmach.c -o lowmach -lm
mkdir -p out && ./lowmach out                     # all experiments, ~2 s
mkdir -p figures && python3 plot_lm.py out figures
```
