---
layout: project
title: "Deep Parametric PINN: Von Kármán Wake"
date: 2026-05-05
tech: ["PyTorch", "Python", "DeepXDE"]
excerpt: "Parametric neural solver for Von Kármán vortex shedding using SIREN activations and RAR refinement ($Re \\in [20, 150]$)."
---

![Flow Field Visualization](/assets/images/pinn/hero_flow.png)

## Technical Index
1. [**4D Parametric Domain**](#parametric) — The Universal Solver approach.
2. [**SIREN vs. Tanh**](#siren) — Overcoming gradient saturation.
3. [**RAR Algorithm**](#rar) — Automated mesh-less refinement.
4. [**Optimization**](#optimization) — Hybrid Adam/L-BFGS strategy.
5. [**Hardware**](#hardware) — AMD ROCm acceleration.

<a name="parametric"></a>
## 1. 4D Parametric Domain: The Universal Solver
A standard CFD simulation solves the Navier-Stokes equations for a specific set of boundary conditions and physical constants. In this project, the **Reynolds number ($Re$)** is not a constant, but a fourth input dimension:

$$\mathcal{N}: [x, y, t, Re] \to [u, v, p]$$

### Motivation
By injecting $Re$ into the network, we transform a single-case solver into a **parametric surrogate model**.
* **Computational Efficiency:** Once trained, the network can predict the flow field for any $Re \in [20, 150]$ in milliseconds.
* **Design Optimization:** This allows for real-time sensitivity analysis and exploration of the transition from steady-state flow to unsteady vortex shedding without re-running the training pipeline.

<a name="siren"></a>
## 2. SIREN Architecture vs. Tanh: Breaking the Gradient Barrier
The choice of the activation function is critical for Physics-Informed machine learning. While many PINN implementations use `tanh`, this solver utilizes **Sinusoidal Representation Networks (SIREN)** (`sin` activation).

### Why `sin` is superior to `tanh` for PDEs:
1.  **Gradient Preservation:** The derivatives of a sine function are shifted sines (cosines). This means that the spatial and temporal gradients (and higher-order Hessians required by Navier-Stokes) maintain the same distribution and "energy" as the original signal. `tanh`, conversely, has derivatives that vanish as the input increases, leading to the **vanishing gradient problem** in deep PDE solvers.
2.  **Spectral Bias & High Frequencies:** Standard networks (ReLU/Tanh) suffer from "spectral bias," learning low-frequency components first and struggling with sharp gradients. The periodic nature of **SIREN** allows the network to represent the fine, high-frequency details of the **Von Kármán vortices** and sharp boundary layers much more accurately.
3.  **Hessian Stability:** Since Navier-Stokes involves second-order derivatives ($\nabla^2 \mathbf{u}$), having an activation function like `sin` that is non-zero in its second derivative is essential for stable backpropagation of the physical loss.


<a name="rar"></a>
## 3. Adaptive Training: The RAR Algorithm
In fluid dynamics, the most critical physics occur in small, high-gradient regions (the wake and the boundary layer). A uniform distribution of collocation points is computationally wasteful.

### Residual-based Adaptive Refinement (RAR)
The implementation uses a 4-cycle **RAR algorithm**:
1.  **Residual Evaluation:** Every 25,000 iterations, the model evaluates the Navier-Stokes residuals on 50,000 random spatio-temporal points.
2.  **Point Injection:** The top 1,000 points with the highest absolute error—representing zones where the physics is not yet satisfied—are added to the training set.
3.  **Focused Learning:** This effectively creates a "smart mesh" that dynamically densifies in the wake, ensuring high-fidelity results where the flow is most complex.


<a name="optimization"></a>
## 4. Hybrid Optimization & Loss Weighting

### The Two-Stage Pipeline
The training leverages two different optimization philosophies:
* **Phase 1 (Adam):** 100,000 iterations to explore the loss landscape and establish the global flow topology. Adam's stochastic nature is perfect for overcoming local minima during the initial RAR cycles.
* **Phase 2 (L-BFGS):** A second-order optimizer that uses curvature information (Hessian approximation). This is the "precision strike" that drives the physical residuals down to machine epsilon ($10^{-5}$ range), ensuring the solution is strictly physical.

### Loss Hierarchy
The loss function is heavily weighted to enforce boundary integrity:
* **Cylinder No-Slip (Weight 100):** Prioritizes zero-velocity at the cylinder wall, the primary source of vorticity.
* **Continuity (Weight 20):** Enforces incompressibility ($\nabla \cdot \mathbf{u} = 0$) as a hard constraint.


<a name="hardware"></a>
## 5. Hardware Acceleration & Performance
The solver is optimized for **AMD GPU (ROCm)** architectures, utilizing specialized environment overrides to ensure compatibility and performance for large-scale tensor operations:

```python
os.environ["DDE_BACKEND"] = "pytorch"
os.environ["HSA_OVERRIDE_GFX_VERSION"] = "10.3.0" # Target: AMD RDNA2/3
