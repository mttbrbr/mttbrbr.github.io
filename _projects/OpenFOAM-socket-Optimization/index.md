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
To make the results as transparent as possible, this section illustrates the configuration I used for testing. The processor is a Ryzen 7900X3D with 12 cores spread across two different CCDs, 6 of which have an additional 64MB L3 cache. The RAM is DDR5 at 6000MT/s. The operating system used for testing is Ubuntu 22.04 LTS e OpenFOAM ESI v2312.
The case used has 17M cells and is based on the incompressible and stationary motorbike tutorial. It can be found as a media mesh at the following link: https://develop.openfoam.com/committees/hpc/-/tree/develop/incompressible/simpleFoam/HPC_motorbike. The simulation ran for 500 iterations and only one writeStep at 250 in the middle of the simulation.

### Flags
#### Flag -O
We now introduce the general optimization flag "-O", which, followed by a number, introduces general optimizations in the compilation.
- -O0: The compiler does not introduce any optimizations and is normally used in the code debugging phase;
- -O2: This is the standard optimization flag for normal distributed binaries; it is a compromise between performance and the size of the final binary;
- -O3: This is the maximum level of safe optimization. Loop unrolling and predictive vectorization are enabled, and it is the standard for OpenFOAM.
- -Ofast: Removes some IEEE floating-point standard compliance. This type of optimization is not recommended because it can introduce rounding errors and numerical instability.

#### Architectural Flags
These flags allow the compiler to take full advantage of our hardware:
- -march=native: Allows the compiler to read hardware information and use what it finds. In my experience, the compiler version should be taken into account; sometimes, with older operating systems, the included compiler does not contain the information for the latest technologies.
- -mtune=native: This flag allows you to optimize performance on the local CPU without introducing tricks that could make the binary incompatible with other machines.
- -msse4.2, -mavx, -mavx2, -mavx512f: Enable the corresponding technology's instruction sets.

#### Specific Numerical Flags
These flags influence how the processor solves math.
- -ffast-math: Allows the compiler to reorder operations and transform $a/b$ into $a * (1/b)$. This allows for maximum performance but introduces risks due to the required decimal precision.
- -fopenmp: Enables OpenMP, which allows processes to be managed as subprocesses. This may be useful for managing hybrid architectures.

### Setup Clang

<a name="siren"></a>
## 2. Single Precision vs Double Precision
### Introduction
In scientific computing, the choice between Single Precision and Double Precision defines the tradeoff between simulation speed and the accuracy of numerical results. It all depends on the numerical precision with which we store the numbers. OpenFOAM typically uses 64-bit precision, which guarantees 15-17 significant decimal places, which ensures numerical stability and allows for better handling of very high gradients and highly distorted or very small grids. The tradeoff is the need to move a larger amount of data.

### Memory-bound processes and the benefits of single precision
In computational fluid dynamics, the main limitation is bandwidth. Moving data is a very expensive process; the CPU performs calculations every 0.5 nanoseconds, while RAM takes up to 50-100 nanoseconds to provide the necessary data. Furthermore, CPUs benefit from the methods described above to perform multiple operations simultaneously.
Using single precision reduces the size of numbers, which go from 64 bits to 32 bits, significantly increasing the amount of data moved between RAM and the CPU. Problems arise when the solver doesn't have enough significant digits to continue iterating and converging.

### The best of both worlds
OpenFOAM allows the use of a "-spdp" compiler flag, which is a compromise between the two. Data stored in RAM is in 32-bit format, but the solver solves linear systems using 64-bit precision. This allows for good numerical stability in solvers while optimizing data transfer between RAM and the CPU.

### Setup

```
export WM_PRECISION_OPTION=SP

export WM_LABEL_SIZE=32
```
### Comparison and Results

| **Compiler Configuration**     | **Time (s)**  | Average Iteration Time| **Speedup** | **Notes**                |
| Standard (APT)                 | 8825.87       |                       | Baseline    | Compatibilità universale |
| gcc -znver4                    | 8656.06       |                       |             |                          |
| gcc -Ofast                     |               |                       |             |                          |
| gcc -znver4 -ffast-math        |               |                       |             |                          |
| Clang -O3                      |               |                       |             |                          |
| Clang -                        |               |                       |             |                          |
| gcc (spdp)                     | 6398.86       |                       |             |                          |


![Residuals baseline]({{ "/assets/projects/OpenFOAM-socket-Optimization/residuals_combined.png" | relative_url }})
![Iterations baseline]({{ "/assets/projects/OpenFOAM-socket-Optimization/iterations_grid.png" | relative_url }})