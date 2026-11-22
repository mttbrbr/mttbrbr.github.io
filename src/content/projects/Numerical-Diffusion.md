---
title: "What Numerical Diffusion Really Is"
date: 2026-11-22
tech: ["C", "CFD", "Numerical Methods"]
excerpt: "Numerical diffusion is not a vague error: it is a real viscosity that the scheme adds to your equations, and it can be predicted and measured to three digits. A few lines of C show how big it is and where it comes from."
---

## Technical Index
1. [**The equation you are really solving**](#modified) — The modified equation of the upwind scheme
2. [**Three schemes, one profile**](#schemes) — Diffusion, dispersion and limiters
3. [**Measuring the viscosity**](#measure) — Theory and experiment, to three digits
4. [**Accuracy and order**](#order) — What "second order" means near a discontinuity
5. [**False diffusion in 2D**](#false) — When the flow is not aligned with the grid
6. [**What this means for real CFD**](#realcfd) — Effective Reynolds number and scheme choice

Everybody who has run a CFD simulation has heard of *numerical diffusion*: "upwind is diffusive", "use a second-order scheme", "the result is smeared". It is often described as a generic inaccuracy, something like noise. It is not. Numerical diffusion is a **physical viscosity** that the discretization adds to the equations, with a precise value that depends on the mesh size, the velocity and the time step. In this article I derive it, then measure it with a small C program, and show that the two agree to three significant digits.

<a name="modified"></a>
## 1. The equation you are really solving

Take the simplest transport equation there is, pure advection of a quantity $u$ with constant velocity $a > 0$:

$$
\frac{\partial u}{\partial t} + a\,\frac{\partial u}{\partial x} = 0
$$

Its exact solution just moves the initial profile at speed $a$ without changing its shape: nothing spreads, nothing decays. There is no diffusion in this equation.

The first-order **upwind** scheme approximates the space derivative with the value on the upstream side:

$$
u_i^{n+1} = u_i^n - C\,\left(u_i^n - u_{i-1}^n\right), \qquad C = \frac{a\,\Delta t}{\Delta x}
$$

where $C$ is the CFL (Courant) number. Expanding $u_{i-1}$ and $u^{n+1}$ in Taylor series and keeping the first terms that the scheme gets wrong, one finds that the discrete solution satisfies, to leading order, a *different* equation:

$$
\frac{\partial u}{\partial t} + a\,\frac{\partial u}{\partial x} = \nu_{num}\,\frac{\partial^2 u}{\partial x^2}, \qquad \nu_{num} = \frac{a\,\Delta x}{2}\,(1 - C)
$$

This is the **modified equation**. The scheme does not solve pure advection: it solves an **advection–diffusion** equation, with a diffusion coefficient $\nu_{num}$ that does not exist in the physics. That is numerical diffusion, and the formula already says a lot:

- it is proportional to $\Delta x$: halving the cell size halves it, but never removes it;
- it is proportional to the velocity: faster flow, more numerical viscosity;
- it depends on the time step through $(1 - C)$: it vanishes at $C = 1$, where upwind moves the profile exactly one cell per step, and it is **largest for small time steps**. For an implicit steady-state solver ($C \to 0$ in practice) it is simply $a\,\Delta x / 2$.

The last point is counterintuitive: reducing the time step makes this scheme *more* diffusive, not less.

<a name="schemes"></a>
## 2. Three schemes, one profile

To see it, I wrote [numdiff.c](/assets/projects/Numerical-Diffusion/numdiff.c), which advects a profile on a periodic domain $[0, 1)$ with $a = 1$ for exactly one revolution ($t = 1$): the exact solution at the end is the initial profile. The profile contains a square wave (a discontinuity) and a Gaussian (a smooth bump, standard deviation $\sigma_0 = 0.04$).

I compare three schemes, all written as one **flux-limited** family, which makes the difference between them a single function $\phi$:

$$
F_{i+1/2} = a\,u_i + \frac{a}{2}\,(1 - C)\,\phi(r_i)\,\left(u_{i+1} - u_i\right), \qquad r_i = \frac{u_i - u_{i-1}}{u_{i+1} - u_i}
$$

| Scheme | Limiter $\phi(r)$ | Order | Modified equation adds |
|---|---|---|---|
| upwind | 0 | 1 | diffusion $\propto \partial^2 u/\partial x^2$ |
| Lax-Wendroff | 1 | 2 | dispersion $\propto \partial^3 u/\partial x^3$ |
| van Leer (TVD) | $(r + \lvert r\rvert)/(1 + \lvert r\rvert)$ | 2 where smooth | switches between the two |

The second term of the flux is an *anti-diffusive* correction: Lax-Wendroff adds it everywhere and cancels the numerical diffusion exactly, van Leer adds it only where the solution is smooth (where $r \approx 1$) and switches it off near jumps and extrema.

![Profiles after one revolution with the three schemes](/assets/projects/Numerical-Diffusion/profiles.svg)

| Scheme | Minimum of $u$ | Maximum of $u$ | Gaussian peak |
|---|---|---|---|
| exact | 0 | 1 | 1.00 |
| upwind | 0.0002 | 0.954 | 0.62 |
| Lax-Wendroff | −0.231 | 1.232 | 0.99 |
| van Leer | 0.0000 | 1.000 | 0.94 |

*200 cells, CFL = 0.5, one revolution. Values of $u$ are dimensionless (the initial profile goes from 0 to 1).*

- **Upwind** is smooth and never creates new maxima or minima, but it *diffuses*: the square wave becomes a bell, and the Gaussian loses about 37% of its peak while getting wider. This is exactly what a physical viscosity would do.
- **Lax-Wendroff** has no numerical diffusion, but its error is *dispersive*: waves of different wavelength travel at slightly different speeds, and around the discontinuity they separate into oscillations, with values down to −0.23 and up to 1.23. For a concentration, a temperature or a turbulent kinetic energy, these are unphysical values that can make a solver diverge.
- **van Leer** keeps the square wave sharp and bounded between 0 and 1. The price is visible on the Gaussian: near the peak the limiter switches back to upwind, so it loses 6% of its height ("clipping").

This is **Godunov's theorem** at work: a *linear* scheme that never creates new extrema can be at most first order. You can have monotone or accurate, not both, unless the scheme is nonlinear, i.e. it adapts to the solution like a limiter does. Every TVD/NVD scheme in a CFD code (`limitedLinear`, `vanLeer`, `MUSCL`, `SuperBee` in OpenFOAM) is a version of this idea.

<a name="measure"></a>
## 3. Measuring the viscosity

If numerical diffusion is really a viscosity, it should be measurable like one. A Gaussian under pure diffusion stays Gaussian, and its variance grows linearly in time:

$$
\sigma^2(t) = \sigma_0^2 + 2\,\nu\,t \quad\Rightarrow\quad \nu = \frac{\sigma^2(t) - \sigma_0^2}{2\,t}
$$

So the program advects the Gaussian with upwind for one revolution, computes the variance of the discrete solution before and after, and turns the growth into a viscosity. No fitting, no free parameter. Then it compares it with the modified equation, $\nu_{num} = \Delta x\,(1 - C)/2$ (here $a = 1$ and the domain length is 1, so everything is dimensionless).

![Measured numerical viscosity vs CFL number for four grids](/assets/projects/Numerical-Diffusion/viscosity.svg)

| Cells $N$ | CFL $C$ | $\nu$ predicted | $\nu$ measured | Measured / predicted |
|---|---|---|---|---|
| 100 | 0.1 | 4.500e-3 | 4.559e-3 | 1.013 |
| 100 | 0.5 | 2.500e-3 | 2.501e-3 | 1.000 |
| 100 | 0.9 | 5.000e-4 | 4.995e-4 | 0.999 |
| 100 | 1.0 | 0 | 0 | exact |
| 200 | 0.1 | 2.250e-3 | 2.251e-3 | 1.000 |
| 200 | 0.5 | 1.250e-3 | 1.250e-3 | 1.000 |
| 400 | 0.1 | 1.125e-3 | 1.125e-3 | 1.000 |
| 800 | 0.1 | 5.625e-4 | 5.625e-4 | 1.000 |
| 800 | 0.9 | 6.250e-5 | 6.251e-5 | 1.000 |

*Dimensionless viscosity ($a = 1$, domain length 1). Full table in [viscosity.csv](/assets/projects/Numerical-Diffusion/viscosity.csv).*

The measured value matches the formula to **three or four significant digits** on every grid. The only visible deviation, 1.3% at $N = 100$, $C = 0.1$, is the coarsest and most diffused case, where the next terms of the Taylor expansion start to count. At $C = 1$ the Gaussian comes back unchanged: zero viscosity, as predicted.

This is what "numerical diffusion" really means. It is not an unpredictable error spread over the solution: the upwind scheme **solves the advection–diffusion equation with a well-defined viscosity**, and the solution you get is the exact physical behaviour of a fluid with that viscosity. The plot also shows the time-step paradox: on every grid, the smaller the CFL number, the more viscous the result.

<a name="order"></a>
## 4. Accuracy and order

The same program measures the $L_1$ error after one revolution on grids from 50 to 1600 cells, separately for the smooth Gaussian and for the square wave.

![L1 error vs number of cells](/assets/projects/Numerical-Diffusion/convergence.svg)

| Profile | Scheme | L1 error, N = 100 | L1 error, N = 1600 | Observed order |
|---|---|---|---|---|
| Gaussian | upwind | 6.6e-2 | 8.7e-3 | 0.88 |
| Gaussian | Lax-Wendroff | 2.6e-2 | 1.2e-4 | 2.00 |
| Gaussian | van Leer | 9.5e-3 | 2.8e-5 | 2.19 |
| square | upwind | 1.1e-1 | 2.8e-2 | 0.50 |
| square | Lax-Wendroff | 7.4e-2 | 1.5e-2 | 0.60 |
| square | van Leer | 3.4e-2 | 4.4e-3 | 0.74 |

*L1 error: average absolute difference from the exact solution (dimensionless). Observed order: slope between N = 800 and N = 1600, i.e. the error is divided by $2^{p}$ when the cells are doubled.*

Two lessons:

- On a smooth profile the formal order is what you get: doubling the cells divides the upwind error by 2 and the second-order errors by 4. At 1600 cells upwind is still **300 times less accurate** than van Leer.
- On a discontinuity **every scheme drops below first order**. The error is concentrated in the few cells around the jump, and a smeared jump of width $w$ costs an error proportional to $w$. For upwind the smearing width grows like $\sqrt{\nu_{num}\,t} \propto \sqrt{\Delta x}$, hence order 0.5. "Second order" is a promise about smooth regions, not about shocks, interfaces or sharp fronts.

<a name="false"></a>
## 5. False diffusion in 2D

In more than one dimension a new effect appears. Consider a steady 2D flow at an angle $\theta$ to the grid, carrying a step: $u = 1$ on one side of a streamline, $u = 0$ on the other. The exact solution keeps the step perfectly sharp along the streamline. With upwind on a square grid of side $h$ and velocity components $(U\cos\theta, U\sin\theta)$, the modified equation has a diffusion term in each grid direction, and part of it acts **across** the streamlines. This is the classic *false diffusion* described by de Vahl Davis and Mallinson in the 1970s.

The program solves this problem on the unit square, with the step entering at the corner. With positive velocity components, upwind can be solved in a single sweep through the grid.

![Steady upwind solution of a step advected at 0, 22.5 and 45 degrees](/assets/projects/Numerical-Diffusion/false_diffusion_fields.png)

*200 × 200 cells. Light: $u = 1$, dark: $u = 0$. Dashed line: position of the exact (sharp) step.*

- At **0°**, with the flow aligned with the grid, the step stays perfectly sharp: zero cells smeared at the outlet. The numerical diffusion is all along the flow, where the solution does not change.
- At **22.5°** and **45°** the step spreads more and more as it travels. The same scheme, the same mesh, the same physics: only the angle between the flow and the grid changed.

At 45° the update of every cell is simply the average of the cell on its left and the one below it. That is exactly the rule of a **random walk** moving left or down with equal probability: the value in a cell is the probability that such a walk reaches the "$u = 1$" boundary first. Numerical diffusion here is literally a random walk on the mesh, and like any random walk its spread grows with the square root of the number of steps. The modified equation predicts it: at 45° the numerical viscosity is $U h / (2\sqrt{2})$ in every direction, so after a distance $s$ the step should have a Gaussian profile with variance $s\,h/\sqrt{2}$ and a 10–90% width of $2.563\,\sqrt{s\,h/\sqrt{2}}$.

![Width of the smeared step vs distance, measured and predicted](/assets/projects/Numerical-Diffusion/falsediff.svg)

| Cells per side | Width at $s = 0.1$ | Width at $s = 0.4$ | Width at $s = 1.0$ | Predicted at $s = 1.0$ |
|---|---|---|---|---|
| 100 | 0.068 | 0.137 | 0.216 | 0.215 |
| 200 | 0.048 | 0.097 | 0.153 | 0.152 |
| 400 | 0.034 | 0.068 | 0.108 | 0.108 |

*Width: distance across the flow between the points where $u = 0.1$ and $u = 0.9$. $s$: distance travelled from the corner. Both as a fraction of the domain size.*

Again, prediction and measurement agree within 1%. And the numbers are large: on a 100 × 100 grid a sharp interface becomes **more than a fifth of the domain wide** after crossing it. Refining the mesh helps slowly: 4 times more cells per side (16 times more cells) only halves the width, because it scales with $\sqrt{h}$.

<a name="realcfd"></a>
## 6. What this means for real CFD

**Compare the numerical viscosity with the physical one.** Take air ($\nu \approx 1.5\cdot10^{-5}$ m²/s) at $U = 10$ m/s on a mesh with 1 mm cells, solved with first-order upwind in a steady (implicit) solver:

$$
\nu_{num} \approx \frac{U\,\Delta x}{2} = \frac{10 \cdot 0.001}{2} = 5\cdot10^{-3}\ \text{m}^2/\text{s} \approx 330\,\nu_{air}
$$

On a car 4 m long, the physical Reynolds number is $U L / \nu \approx 2.7\cdot10^{6}$. The flow that pure upwind would actually compute, wherever the velocity varies along the flow, behaves like one at about $U L/\nu_{num} \approx 8\,000$. These are order-of-magnitude numbers (in a RANS simulation the turbulent viscosity is also hundreds of times the molecular one, and numerical diffusion only acts where the solution has gradients), but they show the point: with a first-order scheme, the mesh can easily be more viscous than the physics you want to simulate. Wakes decay too fast, vortices fade, separation moves, recirculation zones shrink.

What follows in practice:

- **Use upwind to start, not to finish.** First-order upwind is extremely robust precisely because it is diffusive: the numerical viscosity damps everything. A common strategy is to start a difficult case with upwind and then switch to a second-order bounded scheme (`linearUpwind`, `limitedLinear`, `LUST` in OpenFOAM) for the results.
- **Bounded quantities need limiters.** For quantities that must stay in a range (volume fractions, $k$, $\varepsilon$, $\omega$, species), the overshoots of an unlimited second-order scheme like the −0.23 above can make the solver diverge. That is why the turbulence equations are often discretized with bounded or upwind-biased schemes even when the momentum equation is second order.
- **Align the mesh with the flow when you can.** False diffusion is zero along the grid lines and maximal at 45°. Prism layers and hexahedral cells aligned with a boundary layer or a wake add little cross-stream diffusion; a tetrahedral mesh or a structured mesh crossed at an angle adds much more.
- **Refinement is slow medicine.** With a first-order scheme the error falls like $\Delta x$, and the width of a smeared front only like $\sqrt{\Delta x}$. A second-order scheme on a coarser mesh is usually far cheaper than first-order on a finer one.
- **Smaller time steps do not fix it.** For explicit upwind the numerical viscosity even grows as the CFL number decreases. Accuracy in space and in time are separate problems.
- **Numerical diffusion is not always the enemy.** Some approaches use it on purpose: in *implicit LES* the dissipation of the scheme plays the role of the subgrid-scale model. The difference is that it is chosen and controlled, not suffered.

The common thread is the same as in the [previous](/projects/data-locality-cfd/) [articles](/projects/data-between-cores/): a simulation is a model of a model. The equations are discretized, and the discretization changes the equations in a precise way. Knowing *how* is what turns a colourful picture into a result you can trust.

### Reproduce it
All the files are available: [numdiff.c](/assets/projects/Numerical-Diffusion/numdiff.c), [plot_nd.py](/assets/projects/Numerical-Diffusion/plot_nd.py) and the data [profiles.csv](/assets/projects/Numerical-Diffusion/profiles.csv), [convergence.csv](/assets/projects/Numerical-Diffusion/convergence.csv), [viscosity.csv](/assets/projects/Numerical-Diffusion/viscosity.csv), [falsediff.csv](/assets/projects/Numerical-Diffusion/falsediff.csv).

```bash
gcc -O3 -march=native numdiff.c -o numdiff -lm
mkdir -p out && ./numdiff out        # all three experiments, ~0.3 s
mkdir -p figures && python3 plot_nd.py out figures   # SVG charts + 2D field images (needs ImageMagick)
```
