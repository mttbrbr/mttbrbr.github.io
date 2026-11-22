---
title: "Hyperviscosity in DNS: More Reynolds Number for Free?"
date: 2026-11-22
tech: ["C", "DNS", "Turbulence", "FFTW"]
excerpt: "Replacing the viscous term with a higher power of the Laplacian pushes all the dissipation to the smallest resolved scales: on the same 64³ grid, the equivalent Reynolds number goes from 47 to 123. A pseudo-spectral 3D turbulence code shows what is gained, and the price: a bottleneck that piles energy up near the cutoff, and equations that are no longer Navier-Stokes."
---

## Technical Index
1. [**Why DNS is expensive**](#cost) — The Kolmogorov scale sets the grid
2. [**Hyperviscosity**](#hyper) — A different dissipation term
3. [**The experiment**](#experiment) — Forced isotropic turbulence, pseudo-spectral
4. [**What is gained**](#gain) — A longer inertial range on the same grid
5. [**The price: the bottleneck**](#bottleneck) — Energy piling up at the cutoff
6. [**What this means in practice**](#practice) — When it is legitimate, and when it is not

<a name="cost"></a>
## 1. Why DNS is expensive

A Direct Numerical Simulation resolves every scale of a turbulent flow, down to the smallest eddies where viscosity turns kinetic energy into heat. In three-dimensional turbulence energy is injected at a large scale $L$, cascades to smaller and smaller eddies, and is dissipated around the **Kolmogorov scale**

$$
\eta = \left(\frac{\nu^3}{\varepsilon}\right)^{1/4}
$$

where $\nu$ is the viscosity and $\varepsilon$ the rate at which energy is dissipated. Between $L$ and $\eta$ lies the **inertial range**, where energy is only transferred, not dissipated, and the spectrum follows Kolmogorov's law $E(k) = C_K\,\varepsilon^{2/3} k^{-5/3}$, with $C_K \approx 1.6$.

A DNS needs cells of the order of $\eta$; in a spectral code the usual criterion is $k_{max}\,\eta \gtrsim 1$, where $k_{max}$ is the largest resolved wavenumber. The ratio $L/\eta$ grows like $Re^{3/4}$, so the number of points grows like $Re^{9/4}$, and since the time step shrinks with the cell size, the total cost grows roughly like $Re^3$. In terms of the Taylor-scale Reynolds number $Re_\lambda$ often used in turbulence research, the grid size grows like $Re_\lambda^{3/2}$ per direction. Doubling the range of scales is expensive.

<a name="hyper"></a>
## 2. Hyperviscosity

Look at what the viscous term does to a single Fourier mode of wavenumber $k$: it damps it at a rate $\nu k^2$. The damping grows with $k$, but slowly: it already acts at the large and intermediate scales, eating into the inertial range.

**Hyperviscosity** replaces the Laplacian with a power of it:

$$
\nu \nabla^2 \mathbf{u} \quad\longrightarrow\quad -\nu_p\,(-\nabla^2)^p\,\mathbf{u}, \qquad \text{damping rate } \nu_p k^{2p}
$$

With $p = 1$ it is ordinary viscosity. With $p = 2, 4, 8$ the damping rate grows like $k^4$, $k^8$, $k^{16}$: almost zero over most of the spectrum, then a wall right at the smallest resolved scales. The hope is a longer inertial range on the same grid, as if the Reynolds number were much higher.

The equations are then **no longer the Navier–Stokes equations**: there is no physical fluid with this dissipation. A simulation with hyperviscosity is, strictly speaking, a model of the small scales, not a DNS. The question is how much this matters.

<a name="experiment"></a>
## 3. The experiment

[hyper3d.c](/assets/projects/Hyperviscosity/hyper3d.c) is a compact pseudo-spectral solver for **forced, homogeneous, isotropic turbulence** in a periodic cube, the standard setting of turbulence DNS:

- velocity represented by Fourier modes, nonlinear term computed in physical space with FFTW (rotational form), 2/3-rule dealiasing, projection onto divergence-free fields;
- time integration with fourth-order Runge–Kutta and an integrating factor, so the dissipation term is integrated exactly and even $k^{16}$ does not limit the time step;
- forcing at the largest scales ($|k| \le 2.5$) with a constant energy injection rate $\varepsilon = 0.1$;
- spectra averaged over $15 \le t \le 30$, about five large-eddy turnover times after the flow has reached a statistically steady state.

Five runs, same initial condition and same forcing:

- **p = 1 on 64³**, with the viscosity chosen so that $k_{max}\eta \approx 1$: a properly resolved, low-Reynolds-number DNS;
- **p = 2, 4, 8 on 64³**, with $\nu_p$ chosen so that every run damps the largest resolved wavenumber at the same rate ($\nu_p k_{max}^{2p} = \nu_1 k_{max}^2$);
- **p = 1 on 96³** with a lower viscosity, again at $k_{max}\eta \approx 1$: a real DNS at a higher Reynolds number, as a reference.

![Kinetic energy vs time](/assets/projects/Hyperviscosity/energy.svg)

All runs reach a steady state in which the dissipated energy matches the injected one within 2–5%.

<a name="gain"></a>
## 4. What is gained

![Energy spectra](/assets/projects/Hyperviscosity/spectra.svg)

| Run | Re_λ (equivalent) | Mean energy | k_max η (equivalent) |
|---|---|---|---|
| p = 1, 64³ | 47 | 0.53 | 1.06 |
| p = 2, 64³ | 89 | 0.59 | 0.47 |
| p = 4, 64³ | 110 | 0.61 | 0.36 |
| p = 8, 64³ | 123 | 0.62 | 0.31 |
| p = 1, 96³ | 61 | 0.51 | 1.03 |

*For the hyperviscous runs, "equivalent" quantities use the ordinary viscosity that would dissipate the same energy with the same spectrum, $\nu_{eff} = \varepsilon / (2\int k^2 E\,dk)$.*

The ordinary DNS on 64³ has a Taylor Reynolds number of 47 and **no inertial range at all**: its spectrum bends down from $k \approx 4$, well before the cutoff. With hyperviscosity, on exactly the same grid, the spectrum stays close to $k^{-5/3}$ almost up to $k_{max}$, and the equivalent Reynolds number rises to 89, 110 and 123.

The most telling picture is *where* the energy is dissipated:

![Dissipation spectrum, normalized](/assets/projects/Hyperviscosity/dissipation.svg)

With ordinary viscosity the dissipation peaks at $k \approx 5$ and is spread over the whole resolved range, overlapping the scales that should only transfer energy. On a grid this small, viscosity is everywhere. With $p = 8$ the dissipation is essentially zero below $k = 15$ and concentrated in the last few shells before the cutoff.

How expensive would the same Reynolds number be with real viscosity? At fixed $k_{max}\eta$ the grid scales like $Re_\lambda^{3/2}$, and this run confirms it: the 96³ reference reached $Re_\lambda = 61$, against the 62 predicted by scaling up the 64³ run. Reaching the equivalent $Re_\lambda = 123$ of the $p = 8$ run would need a grid of about $270^3$, **75 times more points**, and with the smaller time step roughly **300 times the cost**. That is the attraction.

<a name="bottleneck"></a>
## 5. The price: the bottleneck

The spectrum multiplied by $\varepsilon^{-2/3} k^{5/3}$ should be flat at the Kolmogorov constant in an ideal inertial range. It is not:

![Compensated spectra](/assets/projects/Hyperviscosity/compensated.svg)

| Run | Peak of the compensated spectrum | at k ≈ |
|---|---|---|
| p = 1, 64³ | no plateau, falls from k ≈ 4 | – |
| p = 2, 64³ | 2.9 | 8 |
| p = 4, 64³ | 3.3 | 12 |
| p = 8, 64³ | **4.7** | 16–18 |
| p = 1, 96³ | 2.4 (short plateau) | 4–5 |

*Kolmogorov constant ≈ 1.6. Peaks below k ≈ 3 are in the forcing range.*

With hyperviscosity the spectrum rises **above** the Kolmogorov level before the cutoff, and the higher $p$, the larger and sharper the bump: with $p = 8$ the energy near $k_{max}$ is about three times what an ideal cascade would carry. This is the **bottleneck effect**. The cascade transfers energy to smaller scales through interactions between eddies of similar size. When dissipation switches on abruptly, the modes just below the cutoff lose the smaller partners they would pass energy on to, the transfer slows down there, and energy accumulates.

A milder bottleneck exists in ordinary DNS and in experiments too, just before the dissipation range. Hyperviscosity amplifies it strongly. In the limit $p \to \infty$ the dissipation becomes a sharp truncation of the Fourier series, and the truncated system tends to *thermalize*: the smallest scales fill up towards energy equipartition instead of following a cascade (Frisch et al., 2008, studied exactly this link between hyperviscosity, Galerkin truncation and bottlenecks).

So the extended "inertial range" is only partly genuine. Its first part is close to Kolmogorov; its last part is a pile-up that real turbulence at that Reynolds number would not show.

<a name="practice"></a>
## 6. What this means in practice

- **A simulation with hyperviscosity is not a DNS.** It solves modified equations. It can be a very useful model, but results should be reported as such, not as a Navier–Stokes solution at the equivalent Reynolds number.
- **It is useful when the question is about the large and intermediate scales.** Studies of the inertial range, of large-scale dynamics, of mixing and transport at scales well above the cutoff have used hyperviscosity to get more scale separation for the same cost (Borue and Orszag in the 1990s, many later works in astrophysical and geophysical turbulence).
- **It is wrong for small-scale statistics.** Dissipation-range spectra, velocity-gradient statistics, intermittency of the dissipation, anything that depends on how energy is finally removed, is controlled by the artificial term, not by physics.
- **Moderate powers are a compromise.** $p = 2$ gives less extension than $p = 8$ but a much smaller bottleneck. Higher is not better.
- **Check the bottleneck.** A compensated spectrum is the simplest diagnostic: a plateau is the goal, a bump near the cutoff is the artefact. Discard the affected wavenumbers when fitting scaling laws.
- **The same idea appears elsewhere under other names.** Weather and climate models routinely use $\nabla^4$ diffusion as a numerical filter, spectral vanishing viscosity adds dissipation only to the highest modes, and many finite-volume schemes for compressible flow include a fourth-difference artificial dissipation. In all of them the trade-off is the same: protect the resolved scales, accept artefacts at the smallest ones.

The thread of the series again: the cheapest way to get a better-looking result is often to change the equations a little. Sometimes that is exactly the right model; it is never free, and it has to be said.

### Reproduce it
All the files are available: [hyper3d.c](/assets/projects/Hyperviscosity/hyper3d.c) (needs FFTW 3), [plot3d.py](/assets/projects/Hyperviscosity/plot3d.py) with [svgchart.py](/assets/projects/Hyperviscosity/svgchart.py), the averaged spectra and energy histories (`spectrum3d_*.csv`, `history3d_*.csv`) and a [summary](/assets/projects/Hyperviscosity/summary.txt) of the five runs.

```bash
gcc -O3 -march=native hyper3d.c -o hyper3d -lfftw3 -lm
mkdir -p out
./hyper3d out 64 1 8e-3 30 15      # Navier-Stokes, 64^3          ~5 min
./hyper3d out 64 8 8e-3 30 15      # hyperviscosity p = 8, 64^3  ~5 min (also p = 2, 4)
./hyper3d out 96 1 4.6e-3 30 15    # Navier-Stokes, 96^3          ~25 min
mkdir -p figures && python3 plot3d.py out figures
```
