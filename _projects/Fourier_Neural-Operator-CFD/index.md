---
layout: project
title: "Fourier Neural Operators for CFD"
date: 2026-05-18
tech: ["Python", "PyTorch", "Deep Learning", "CFD"]
excerpt: "Exploring the theoretical foundations and application of Fourier Neural Operators (FNO) and Geo-FNO for modeling complex fluid dynamics."
---

## Technical Index
1. [**The Theory of Neural Operators**](#fno-theory) — Mapping between infinite-dimensional function spaces
2. [**Spectral Convolution in FNO**](#spectral-conv) — The power of the Fourier domain
3. [**Handling Complex Geometries**](#geofno) — Latent space mapping with Geo-FNO
4. [**Masked Formulation for Obstacles**](#masked-fno) — Adapting FNO for Von Kármán vortex streets

<a name="fno-theory"></a>
## 1. The Theory of Neural Operators
### Operator Learning vs. Function Approximation
Standard deep learning architectures (like MLPs or CNNs) are designed to learn mappings between finite-dimensional Euclidean spaces ($\mathbb{R}^n \to \mathbb{R}^m$). In Computational Fluid Dynamics (CFD), solving the Navier-Stokes equations intrinsically means finding an operator $\mathcal{G}$ that maps an initial condition (or boundary condition) function $a(x)$ to a solution function $u(x)$.

**Neural Operators** extend the universal approximation theorem to infinite-dimensional function spaces. Instead of learning a mapping for a specific grid resolution, Neural Operators learn the continuous operator itself. This leads to the defining property of Neural Operators: **Mesh Independence**. 

Because the network learns the continuous physics mapping rather than a discrete grid representation, a model trained on a coarse $32 \times 32$ grid can seamlessly infer solutions on a $256 \times 256$ grid without any retraining or interpolation steps. This is a profound shift from traditional finite-volume or finite-element solvers.

<a name="spectral-conv"></a>
## 2. Spectral Convolution in FNO
### The Local vs. Global Receptive Field
Convolutional Neural Networks (CNNs) extract features using small, localized spatial kernels. To model large-scale fluid structures (like massive vortex shedding), a CNN requires dozens of layers to expand its receptive field, suffering from severe numerical dissipation and compounding errors along the way.

**Fourier Neural Operators (FNO)** bypass this by performing the convolution operation in the frequency domain. According to the Convolution Theorem, a convolution in physical space is equivalent to a pointwise multiplication in Fourier space. 

### The FNO Architecture
The core algorithm of the FNO layer consists of three steps:
1. **Fourier Transform ($\mathcal{F}$):** The input field is transformed into the frequency domain using the Fast Fourier Transform (FFT).
2. **Spectral Filtering & Multiplication ($R_\phi$):** In the Fourier domain, we truncate the higher-frequency modes, keeping only the first $k$ modes. The remaining lower modes (which represent the macro-structures of the fluid) are multiplied by learnable complex weight tensors $R_\phi$.
3. **Inverse Fourier Transform ($\mathcal{F}^{-1}$):** The manipulated frequency representation is transformed back into the physical domain.

Mathematically, the iterative update within an FNO layer is represented as:
$$ v^{(t+1)}(x) = \sigma \left( W v^{(t)}(x) + \mathcal{F}^{-1} \big( R_\phi \cdot \mathcal{F}(v^{(t)}) \big)(x) \right) $$
Where $W$ is a local linear transformation (acting as a residual bypass for high-frequency spatial details) and $\sigma$ is a non-linear activation function. By aggressively truncating high-frequency noise and focusing on global low-frequency dynamics, the FNO achieves unprecedented computational efficiency.

<a name="geofno"></a>
## 3. Handling Complex Geometries
### The Limitation of the Standard FNO
The standard FNO relies heavily on the discrete FFT, which strictly assumes a uniform, Cartesian grid with periodic boundaries. Real-world engineering problems (like flow around an airfoil or cylinder) are inherently non-uniform, requiring unstructured meshes that adapt to complex shapes to capture boundary layers accurately. 

### The Geo-FNO Architecture
To bridge the gap between spectral convolutions and arbitrary geometries, this project implements **Geo-FNO**. Geo-FNO introduces an explicit deformation mapping between the physical space and a computational latent space.

1. **Mapping Phase ($\Psi$):** The physical coordinates (from an unstructured or non-uniform mesh) are mapped onto a uniform computational grid in latent space. This mapping can be fixed mathematically or learned via an auxiliary neural network.
2. **Latent FNO:** Once the data is projected onto the uniform Cartesian grid, standard FFT-based FNO layers are applied optimally.
3. **Inverse Mapping ($\Psi^{-1}$):** The calculated solution is projected back to the original, complex physical geometry.

<a name="masked-fno"></a>
## 4. Masked Formulation for Obstacles
When simulating phenomena like the **Von Kármán vortex street** behind a blunt body (such as a cylinder), the physical domain contains "holes" (the obstacle itself) where the fluid equations do not apply. 

Standard FFT operations would artificially smear flow information directly through the solid obstacle, violating the physics. To enforce the impermeability of the obstacle, the architecture is modified into a **Masked FNO**. 

A binary mask $\mathcal{M}(x)$ is introduced, where $1$ represents the fluid domain and $0$ represents the solid geometry. This mask is concatenated with the input features (such as velocity components and pressure) and passed through the network. Furthermore, the loss function is mathematically weighted by this mask:
$$ \mathcal{L} = \frac{1}{N} \sum_{i} \mathcal{M}(x_i) \left\| \hat{u}(x_i) - u_{true}(x_i) \right\|^2 $$
This ensures that the optimizer completely ignores the unphysical region inside the cylinder during backpropagation, forcing the Neural Operator to strictly learn the shedding mechanics and pressure variations occurring in the wake.

### Example: Masked Loss Implementation (PyTorch)

```python
import torch
import torch.nn.functional as F

def masked_l2_loss(pred, target, mask):
    """
    Computes the L2 loss only on the fluid domain.
    mask: 1 for fluid, 0 for solid obstacle.
    """
    # Calculate squared error
    squared_error = (pred - target) ** 2
    
    # Apply the mask
    masked_error = squared_error * mask
    
    # Calculate mean only over the valid fluid points
    # (avoiding division by zero if mask is empty)
    fluid_points = torch.sum(mask) + 1e-8
    loss = torch.sum(masked_error) / fluid_points
    
    return loss
```