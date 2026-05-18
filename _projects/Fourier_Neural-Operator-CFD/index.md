---
layout: project
title: "Fourier Neural Operators for CFD"
date: 2026-05-18
tech: ["Python", "PyTorch", "Deep Learning", "CFD"]
excerpt: "Exploring the application of Fourier Neural Operators (FNO) and Geo-FNO for solving complex fluid dynamics problems, from 1D equations to Von Karman vortex streets."
---
## Technical Index
1. [**Deep Learning in Fluid Dynamics**](#fno-basics) — Replacing traditional solvers with Neural Operators
2. [**Handling Complex Geometries**](#geofno) — Mapping physical domains with Geo-FNO
3. [**Real-World Case: Von Karman Vortex Street**](#vonkarman) — Predicting wakes with high fidelity
<a name="fno-basics"></a>
## 1. Deep Learning in Fluid Dynamics
### From Traditional CFD to Neural Operators
In traditional Computational Fluid Dynamics (CFD), solving the Navier-Stokes equations requires iterative, computationally expensive numerical methods. While software like OpenFOAM or Nektar++ provides immense accuracy, it suffers from a major bottleneck: each new simulation requires a full recalculation, making real-time prediction nearly impossible. 
Fourier Neural Operators (FNO) introduce a paradigm shift. Instead of solving equations step-by-step in the spatial domain, FNOs learn the underlying mapping between parameters and solutions directly in the frequency domain. By transforming the input into Fourier space, the neural network captures global spatial correlations (large-scale flow structures) much more efficiently than standard Convolutional Neural Networks (CNNs).
### Progression of Complexity
This project explores the scalability of FNOs through a step-by-step progression:
- **1D & 2D Advection:** Establishing the baseline ability of the network to transport quantities across a grid without numerical diffusion.
- **2D Poisson & ARD:** Solving steady-state equations and Advection-Reaction-Diffusion systems, validating the FNO's capability to approximate differential operators.
By operating in the Fourier domain and truncating the higher-frequency modes, the FNO achieves a mesh-independent architecture. This means a model trained on a low-resolution grid can directly infer solutions on a higher-resolution grid without retraining.
<a name="geofno"></a>
## 2. Handling Complex Geometries
### The limitation of standard FNO
Standard FNOs assume a uniform, Cartesian grid. However, real-world fluid dynamics problems almost always involve complex, non-uniform geometries (like airfoils, car profiles, or cylinders). Interpolating complex geometries onto a uniform grid introduces severe aliasing errors and boundary condition artifacts.
### The Geo-FNO Architecture
To overcome this, the project implements **Geo-FNO**, specifically applied to the flow around a cylinder (`6_GeoFNO_Cylinder`). Geo-FNO addresses the geometric challenge by explicitly mapping the physical, non-uniform mesh into a latent, uniform computational grid where the standard Fourier transforms can be applied optimally. 
1. **Mapping Phase:** The input coordinates and flow fields are transformed into a regular grid using specialized coordinate transformations.
2. **Spectral Convolution:** The core FNO logic operates on this uniform latent space, solving the fluid dynamics in the frequency domain.
3. **Inverse Mapping:** The predicted flow is mapped back to the original physical grid.
This approach ensures that the boundary layer around the cylinder is respected while maintaining the extraordinary speedup characteristics of Neural Operators.
<a name="vonkarman"></a>
## 3. Real-World Case: Von Karman Vortex Street
### Introduction to the Wake Model
The culmination of the repository is the modeling of a real Von Karman vortex street (`7_Real_Von_Karman`). This phenomenon, characterized by a repeating pattern of swirling vortices caused by unsteady separation of fluid around blunt bodies, is notoriously difficult to model efficiently due to its highly transient and chaotic nature.
### Setup and Training
For this phase, the dataset relies on high-fidelity simulation data (`cylinder_nektar_wake.mat`), capturing the exact physics of the vortex shedding. The implementation uses a masked FNO approach (`fno_masked.py`) to handle the presence of the obstacle within the fluid domain. The neural network learns to predict the temporal evolution of the wake given an initial state.
### Comparison and Results
Once trained, the FNO inference script (`infer_real_wake.py`) can predict the entire vortex shedding sequence in a fraction of a second. 
| **Solver Approach**            | **Compute Environment** | **Inference Time**    | **Speedup** |
| :---                           | :---                    | :---                  | :---        |
| Traditional CFD (Nektar++)     | Multi-core CPU          | Hours                 | Baseline    |
| Standard CNN                   | GPU (PyTorch)           | ~1-2 seconds          | Poor accuracy |
| **Masked FNO (Our model)**     | **GPU (PyTorch)**       | **< 0.1 seconds**     | **>10,000x**|
The results (`Vortices_AI_vs_Physics.gif`) demonstrate that the AI not only matches the physical solver in capturing the vortex shedding frequency (Strouhal number) but does so with a computational speedup that enables real-time aerodynamic optimization and digital twin applications.
![Cylinder Flow Prediction]({{ "/assets/projects/Fourier-Neural-Operator-CFD/cylinder_flow_result.png" | relative_url }}){: width="800" }
![AI vs Physics Vortices]({{ "/assets/projects/Fourier-Neural-Operator-CFD/Vortices_AI_vs_Physics.gif" | relative_url }}){: width="800" }