---
layout: project
title: "OpenFOAM socket Optimization"
date: 2026-05-05
tech: ["Linux", "C++", "Hardware"]
excerpt: "An in-depth study on optimizing OpenFOAM's computational performance in a Linux environment."
---


## Technical Index
1. [**Compilation and hardware**](#parametric) — When bare metal integration matters
2. [**Single Precision vs Double Precision**](#parametric) — How to choose a compromise
3. [**Processor asymmetry**](#rar) — A comparison of the latest technologies


<a name="parametric"></a>
## 1. Compilation and Hardware

### Compilation
In Linux environments, we're accustomed to using package managers like apt or pacman to install software. These tools are very useful because they allow us to quickly and easily install applications and create dependencies. However, when the need is to achieve maximum performance, using these tools is no longer convenient. In the world of scientific computing, in fact, compiled languages ​​are mostly used—that is, languages ​​that require a phase called compilation before the code can be executed, which transforms our code into machine language.
Compiling is a very delicate process because transforming code, which in OpenFOAM is C++ code, into machine language is not a one-time operation; it's not a simple translation from one language to another, but rather an interpretation. The most widely used compilers are GCC and Clang, and when compiling code, they must balance three main factors: execution speed, binary size, and ease of debugging.
The pursuit of computational performance has led me to study the main types of optimization possible. In particular, it is clear that the compiled distributed package must meet the need for compatibility with as many hardware configurations as possible. While this is the best approach when distributing a package online, it clashes with the requirement that the calculation must be performed securely and in the shortest possible time. A "generalist" compilation prevents the software from being tailored to our hardware. Recent processors have available SIMD (single instruction, multiple data) units that allow the same operation to be performed on multiple data points simultaneously. A specialized compilation can therefore unlock AVX, AVX2, and AX512. These specific technologies allow the processing of 4, 8, or even 16 numbers in a single clock cycle. On the memory side, a hardware-based compilation allows you to avoid making assumptions about the size of the L1, L2, and L3 cache levels and to take advantage of Loop Tiling, reorganizing the code so that the processed data fits perfectly into the cache, reducing the CPU's data wait times. Another important detail concerns the processor structure itself: by compiling for your hardware, you can differentiate the use of cores geared toward computational power and cores geared toward efficiency, reducing the penalty of a highly asymmetric configuration.

### Test configuraton
To make the results as transparent as possible, this section illustrates the configuration I used for testing. The processor is a Ryzen 7900X3D with 12 cores spread across two different CCDs, 6 of which have an additional 64MB L3 cache. The RAM is DDR5 at 6000MT/s. The operating system used for testing is Ubuntu 22.04 LTS.

### Flags 


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
