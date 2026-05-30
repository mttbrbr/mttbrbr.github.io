---
layout: project
title: "Monolithic Fluid Dynamics: Architecting High-Performance Coupled Solvers via OpenFOAM and PETSc Integration"
date: 2026-05-30
tech: ["C++", "OpenFOAM", "PETSc", "Linear Algebra", "MPI", "HPC"]
excerpt: "An advanced software-engineering analysis of monolithic CFD architectures. This project bridges the face-based geometric discretization of OpenFOAM with the algebraic distributed robustness of PETSc to bypass the memory bandwidth wall."
---

![Coupled Solver Concept](/assets/images/coupled-foam/hero_wide.png)

## Table of Contents
1. [Architectural Bottlenecks and the Memory Bandwidth Wall](#architectural-bottlenecks-and-the-memory-bandwidth-wall)
2. [Mathematical Formulation and Collocated Grid Discretization](#mathematical-formulation-and-collocated-grid-discretization)
3. [The Architectural Dilemma of Matrix-Free Shells and Explicit Blocks](#the-architectural-dilemma-of-matrix-free-shells-and-explicit-blocks)
4. [The Face-to-Cell Scattering Penalty and Cache Line Invalidation](#the-face-to-cell-scattering-penalty-and-cache-line-invalidation)
5. [Topology Mapping and Strict Non-Zero Preallocation](#topology-mapping-and-strict-non-zero-preallocation)
6. [Distributed Jacobian Assembly and Implicit Under-Relaxation Mechanics](#distributed-jacobian-assembly-and-implicit-under-relaxation-mechanics)
7. [Asynchronous Parallel Synchronicity and MPI Deadlock Elimination](#asynchronous-parallel-synchronicity-and-mpi-deadlock-elimination)
8. [Flux Reconstruction and Segregated Turbulence Closures](#flux-reconstruction-and-segregated-turbulence-closures)
9. [Advanced Block Preconditioning and Schur Complement Approximations](#advanced-block-preconditioning-and-schur-complement-approximations)
10. [Theoretical Hardware Profiling Strategy and L3 V-Cache Mapping](#theoretical-hardware-profiling-strategy-and-l3-v-cache-mapping)

---

## Architectural Bottlenecks and the Memory Bandwidth Wall

In high-performance Computational Fluid Dynamics, performance scalability is strictly governed by the Roofline model boundary, where execution speed is limited by memory bandwidth rather than raw floating-point computing capacity. Standard segregated algorithms, such as the widely implemented SIMPLE or PISO loops, resolve the incompressible Navier-Stokes equations by completely decoupling the momentum equations from the continuity constraint. This mathematical splitting treats pressure and velocity as isolated linear systems solved in a sequential, staggered fashion. 

The primary consequence of this decoupling is a severe numerical lag in the kinematic coupling. Because the velocity field is updated using an outdated pressure field from the previous iteration, engineers must apply aggressive empirical under-relaxation factors to prevent catastrophic numerical divergence, especially when dealing with non-orthogonal or highly skewed meshes. 

From an architectural standpoint, these nested loops exhibit a very low arithmetic intensity, which represents the ratio of floating-point operations to data bytes transferred. The CPU core complexes waste billions of clock cycles waiting for fragmented matrix arrays, velocity vectors, and pressure fields to stream across the memory bus, which starves the floating-point units. 

The monolithic approach natively addresses this hardware inefficiency. By assembling velocity and pressure fields simultaneously inside a single, unified global algebraic system, the exact physical coupling is embedded directly within the global Jacobian matrix. This eliminates the staggered outer iterations and decreases the non-linear iterations required to reach strict convergence, shifting the computation toward higher arithmetic intensity.

## Mathematical Formulation and Collocated Grid Discretization

The continuous incompressible Navier-Stokes equations for a Newtonian fluid govern the conservation of momentum and mass through a system of partial differential equations:

$$\nabla \cdot (\mathbf{u} \otimes \mathbf{u}) - \nabla \cdot (\nu \nabla \mathbf{u}) + \nabla p = 0$$

$$\nabla \cdot \mathbf{u} = 0$$

When discretized via the Finite Volume Method on a collocated grid layout where both velocity and pressure are stored at the cell centers, assembling a naive saddle-point system yields a mathematically unstable linear system. The presence of an exact zero block on the main diagonal creates a highly indefinite matrix and triggers spurious spatial oscillations known as the checkerboard effect. This happens because the central-differencing scheme decouples the pressure values at alternate cells, leading to a singular or poorly bounded system.

To resolve this issue within our monolithic infrastructure, we leverage OpenFOAM's spatial discretization modules to extract stabilized operators. Instead of assembling the matrix from scratch, our application constructs underlying segregated matrix equations. This allows us to pull coefficients from the pressure equation, which incorporates the Laplacian operator and explicitly injects Rhie-Chow stabilization directly into the main diagonal block of the monolithic system. 

```cpp
fvVectorMatrix UEqn = fvm::div(phi, U) + turbulence->divDevReff(U);
volScalarField rAU = 1.0 / UEqn.A();
surfaceScalarField rAUf = fvc::interpolate(rAU);
fvScalarMatrix pEqn = fvm::laplacian(rAUf, p);

```

This mathematical transformation changes the continuity equation into a well-posed elliptic formulation:

$$\begin{bmatrix} \mathbf{A}_{uu} & \mathbf{A}_{up} \\ \mathbf{A}_{pu} & \mathbf{A}_{pp} \end{bmatrix} \begin{bmatrix} \mathbf{u} \\ p \end{bmatrix} = \begin{bmatrix} \mathbf{b}_u \\ \mathbf{b}_p \end{bmatrix}$$

In this block framework, $$\mathbf{A}_{uu}$$ represents the unrolled diagonal and off-diagonal coefficients of the velocity transport matrix, while $$\mathbf{A}_{pp}$$ represents the stabilized pressure Laplacian coefficients. The off-diagonal blocks $$\mathbf{A}_{up}$$ and $$\mathbf{A}_{pu}$$ represent the explicit geometric face-interpolated cross-couplings derived from the surface normal vectors.

## The Architectural Dilemma of Matrix-Free Shells and Explicit Blocks

Interfacing OpenFOAM with PETSc forces a critical software engineering trade off between memory footprint and preconditioning capability. For standard Krylov subspace accelerators like GMRES, PETSc does not strictly require explicit knowledge of matrix elements; it only requires the evaluation of the matrix-vector product action. By wrapping OpenFOAM's native Lower-Diagonal-Upper storage inside PETSc's matrix-shell type and remapping the operational callback, we can redirect PETSc's vector pointers directly to OpenFOAM's optimized internal matrix multiplication routine.

This matrix-free path offers the advantage of zero memory duplication since OpenFOAM retains exclusive ownership of its native data layout without any translation overhead. However, this path completely breaks down when applying advanced preconditioning. Modern, scalable Algebraic Multigrid methods and physics-based block field splits are inherently algebraic. They must explicitly inspect the matrix coefficients to construct strength-of-connection graphs, aggregate coarse nodes, and compute coarse-grid triple products. Passing a matrix shell to these systems causes an immediate runtime crash because the preconditioner cannot access the underlying topology.

To unlock multi-level algebraic preconditioners, we must accept the explicit duplication of the matrix coefficients via the explicit block path. The decoupled, face-based Lower-Diagonal-Upper matrix arrays must be completely unrolled and re-mapped into an interleaved cell-based Block Compressed Sparse Row format managed by PETSc. This layout uses a fixed block size of four to contain the interleaved degrees of freedom for the three velocity components and pressure.

This choice introduces a memory duplication overhead of approximately two times the total matrix storage. Rewriting OpenFOAM’s core to natively utilize block storage would require refactoring the entire discretization pipeline, dynamic mesh handling, and boundary conditions, which is an impractical task. Therefore, the monolithic architecture intentionally trades memory capacity for convergence speed.

## The Face-to-Cell Scattering Penalty and Cache Line Invalidation

Choosing the explicit block route introduces a significant operational bottleneck known as the Face-to-Cell Scattering Penalty. OpenFOAM's discretization loops are fundamentally face-oriented, meaning the framework iterates over the mesh faces, computes the flux and transport coefficients once per face, and streams these values directly into sequential array blocks.

Conversely, PETSc's algebraic structures are cell-oriented and row-oriented. Transforming OpenFOAM's coefficients into block entries requires an intense indexing permutation. For every internal face, the solver must look up the owner and neighbour cell IDs, calculate their global interprocess offsets, map the velocity components and pressure component into a contiguous sixteen-element scalar sub-matrix, and invoke the blocked value insertion routine.

This scattering operation introduces non-contiguous memory writes and pointer-chasing across the heap. When coefficients are written to scattered memory locations corresponding to distant cell IDs, it triggers frequent cache line invalidations and Translation Lookaside Buffer misses. To minimize this penalty and make the monolithic solver viable, this conversion overhead must be highly optimized using strict memory preallocation and vectorized block memory transfers during the assembly loop.

## Topology Mapping and Strict Non-Zero Preallocation

Dynamic memory allocation during distributed matrix assembly severely degrades performance. If PETSc is forced to dynamically allocate memory for unexpected non-zero entries during assembly, it triggers heap fragmentation, expensive internal array copies, and implicit synchronization stalls across MPI ranks.

Our implementation addresses this via a dedicated topology mapper class. This utility scans the mesh connectivity, accounting for internal faces, processor-to-processor halo boundaries, and cyclic patches to compute the exact number of diagonal and off-diagonal non-zero blocks for every single cell before any matrix values are inserted.

**File:** `PetscTopologyMapper.H`

```cpp
#ifndef PetscTopologyMapper_H
#define PetscTopologyMapper_H

#include "fvMesh.H"
#include "processorPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "labelList.H"

namespace Foam
{
struct PreallocData
{
    labelList d_nnz; 
    labelList o_nnz; 
};

class PetscTopologyMapper
{
public:
    static PreallocData computeBlockPreallocation(const fvMesh& mesh)
    {
        const label nCells = mesh.nCells();
        PreallocData allocData;
        allocData.d_nnz.setSize(nCells, 1); 
        allocData.o_nnz.setSize(nCells, 0);

        const labelUList& own = mesh.owner();
        const labelUList& nei = mesh.neighbour();

        forAll(nei, faceI)
        {
            allocData.d_nnz[own[faceI]]++;
            allocData.d_nnz[nei[faceI]]++;
        }

        const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
        forAll(boundaryMesh, patchI)
        {
            const polyPatch& patch = boundaryMesh[patchI];
            
            if (isA<processorPolyPatch>(patch))
            {
                const labelUList& faceCells = patch.faceCells();
                forAll(faceCells, i) 
                {
                    allocData.o_nnz[faceCells[i]]++; 
                }
            }
            else if (isA<cyclicPolyPatch>(patch))
            {
                const labelUList& faceCells = patch.faceCells();
                forAll(faceCells, i) 
                {
                    allocData.d_nnz[faceCells[i]]++; 
                }
            }
        }
        return allocData;
    }
};
}
#endif

```

## Distributed Jacobian Assembly and Implicit Under-Relaxation Mechanics

The assembly loop unrolls OpenFOAM's geometric fields and maps them into PETSc's linear structures. To maintain stability, under-relaxation must be handled implicitly within the monolithic matrix itself rather than as a post-processing step. For a cell with a diagonal momentum coefficient $\mathcal{A}_{uu}$ and source term $\mathbf{b}_u$, implicit under-relaxation alters the matrix diagonal by dividing it by the relaxation factor $\alpha_U$, while the right-hand side vector incorporates historical time-step values scaled by the remaining fraction.

The main solver file loops over internal faces, applies interprocess boundary values, and embeds relaxation directly onto the cell diagonals.

**File:** `coupledSimpleFoam.C`

```cpp
// ... Initialization and Setup inside main() ...
const PetscInt bs = 4; 
PetscInt localEqs = mesh.nCells() * bs;
PetscInt globalEqs = indexer.totalCells() * bs;

Mat A;
MatCreate(PETSC_COMM_WORLD, &A);
MatSetSizes(A, localEqs, localEqs, globalEqs, globalEqs);
MatSetType(A, MATBAIJ); 
// ... Preallocation calls using d_nnz and o_nnz arrays ...
MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

forAll(own, faceI)
{
    PetscInt global_o = indexer.global(own[faceI]);
    PetscInt global_n = indexer.global(nei[faceI]);
    vector nA = Sf[faceI];
    scalar w  = weights[faceI];

    vector Aup_upper = (1.0 - w) * nA; vector Apu_upper = (1.0 - w) * nA; 
    vector Aup_lower = -w * nA;        vector Apu_lower = -w * nA;

    Aup_diag[own[faceI]] += w * nA; Aup_diag[nei[faceI]] -= (1.0 - w) * nA;
    Apu_diag[own[faceI]] += w * nA; Apu_diag[nei[faceI]] -= (1.0 - w) * nA;

    PetscScalar b_upper[16] = {0.0};
    b_upper[0] = U_upper[faceI]; b_upper[5] = U_upper[faceI]; b_upper[10] = U_upper[faceI]; 
    b_upper[3] = Aup_upper.x();  b_upper[7] = Aup_upper.y();  b_upper[11] = Aup_upper.z(); 
    b_upper[12]= Apu_upper.x();  b_upper[13]= Apu_upper.y();  b_upper[14]= Apu_upper.z(); 
    b_upper[15]= p_upper[faceI];                                                           
    MatSetValuesBlocked(A, 1, &global_o, 1, &global_n, b_upper, ADD_VALUES);

    PetscScalar b_lower[16] = {0.0};
    b_lower[0] = U_lower[faceI]; b_lower[5] = U_lower[faceI]; b_lower[10] = U_lower[faceI];
    b_lower[3] = Aup_lower.x();  b_lower[7] = Aup_lower.y();  b_lower[11] = Aup_lower.z(); 
    b_lower[12]= Apu_lower.x();  b_lower[13]= Apu_lower.y();  b_lower[14]= Apu_lower.z();  
    b_lower[15]= p_lower[faceI];
    MatSetValuesBlocked(A, 1, &global_n, 1, &global_o, b_lower, ADD_VALUES);
}

forAll(boundaryMesh, patchI)
{
    const polyPatch& patch = boundaryMesh[patchI];
    if (isA<processorPolyPatch>(patch))
    {
        const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(patch);
        label bndFaceOffset = procPatch.start() - mesh.nInternalFaces();

        const vectorField& U_patchCoeffs = UEqn.internalCoeffs()[patchI];
        const scalarField& p_patchCoeffs = pEqn.internalCoeffs()[patchI];
        const vectorField& patchSf = mesh.Sf().boundaryField()[patchI];

        forAll(procPatch.faceCells(), i)
        {
            label localCell = procPatch.faceCells()[i];
            PetscInt global_o = indexer.global(localCell);
            PetscInt global_rem = indexer.remoteGlobal(bndFaceOffset + i);
            vector nA = patchSf[i];
            
            PetscScalar b_mpi[16] = {0.0};
            b_mpi[0] = U_patchCoeffs[i].x(); b_mpi[5] = U_patchCoeffs[i].y(); b_mpi[10]= U_patchCoeffs[i].z();
            b_mpi[3] = nA.x(); b_mpi[7] = nA.y(); b_mpi[11]= nA.z();
            b_mpi[12]= nA.x(); b_mpi[13]= nA.y(); b_mpi[14]= nA.z();
            b_mpi[15]= p_patchCoeffs[i];
            
            MatSetValuesBlocked(A, 1, &global_o, 1, &global_rem, b_mpi, ADD_VALUES);
        }
    }
    else 
    {
        const labelUList& faceCells = patch.faceCells();
        const vectorField& Up_int = UEqn.internalCoeffs()[patchI];
        const vectorField& Up_bnd = UEqn.boundaryCoeffs()[patchI];
        const scalarField& pp_int = pEqn.internalCoeffs()[patchI];
        const scalarField& pp_bnd = pEqn.boundaryCoeffs()[patchI];

        forAll(faceCells, i)
        {
            label cellI = faceCells[i];
            U_bndDiag[cellI] += Up_int[i]; U_bndRHS[cellI]  += Up_bnd[i];
            p_bndDiag[cellI] += pp_int[i]; p_bndRHS[cellI]  += pp_bnd[i];
        }
    }
}

for (label cellI = 0; cellI < mesh.nCells(); ++cellI)
{
    PetscInt global_c = indexer.global(cellI);
    PetscScalar b_diag[16] = {0.0};
    
    scalar Auu_x = U_diag[cellI] + U_bndDiag[cellI].x();
    scalar Auu_y = U_diag[cellI] + U_bndDiag[cellI].y();
    scalar Auu_z = U_diag[cellI] + U_bndDiag[cellI].z();
    scalar App   = p_diag[cellI] + p_bndDiag[cellI];

    b_diag[0]  = Auu_x / alpha_U; b_diag[5]  = Auu_y / alpha_U; b_diag[10] = Auu_z / alpha_U;
    b_diag[3]  = Aup_diag[cellI].x(); b_diag[7]  = Aup_diag[cellI].y(); b_diag[11] = Aup_diag[cellI].z();
    b_diag[12] = Apu_diag[cellI].x(); b_diag[13] = Apu_diag[cellI].y(); b_diag[14] = Apu_diag[cellI].z();
    b_diag[15] = App / alpha_p;

    MatSetValuesBlocked(A, 1, &global_c, 1, &global_c, b_diag, ADD_VALUES);

    PetscInt rows[4] = {global_c*4, global_c*4+1, global_c*4+2, global_c*4+3};
    vector U_old = U[cellI]; scalar p_old = p[cellI];

    PetscScalar vals[4] = {
        U_source[cellI].x() + U_bndRHS[cellI].x() + ((1.0 - alpha_U) / alpha_U) * Auu_x * U_old.x(),
        U_source[cellI].y() + U_bndRHS[cellI].y() + ((1.0 - alpha_U) / alpha_U) * Auu_y * U_old.y(),
        U_source[cellI].z() + U_bndRHS[cellI].z() + ((1.0 - alpha_U) / alpha_U) * Auu_z * U_old.z(),
        p_source[cellI]     + p_bndRHS[cellI]     + ((1.0 - alpha_p) / alpha_p) * App * p_old
    };
    VecSetValues(b, 4, rows, vals, ADD_VALUES);
}
// ... Assembly termination and solver execution paths ...

```

## Asynchronous Parallel Synchronicity and MPI Deadlock Elimination

When running large-scale cases across thousands of MPI ranks, exchanging global cell indices across processor boundaries is highly prone to communication deadlocks. A standard blocking approach where a rank waits for a neighbor to receive data while that neighbor is simultaneously blocked trying to send to a third rank will freeze execution.

To guarantee parallel scalability, our integration relies on OpenFOAM's asynchronous non-blocking streaming buffers. This framework separates data serialization from physical transmission by aggregating the global indexing data into contiguous memory buffers and scheduling asynchronous, non-blocking point-to-point exchanges underneath the hood.

**File:** `PetscIndexer.H`

```cpp
#ifndef PetscIndexer_H
#define PetscIndexer_H

#include "fvMesh.H"
#include "globalIndex.H"
#include "processorPolyPatch.H"
#include "PstreamBuffers.H"

namespace Foam
{
class PetscIndexer
{
    globalIndex globalCells_;
    labelList localToGlobal_;
    labelList bndFaceToRemoteGlobal_;

public:
    PetscIndexer(const fvMesh& mesh)
    :
        globalCells_(mesh.nCells()),
        localToGlobal_(mesh.nCells()),
        bndFaceToRemoteGlobal_(mesh.nBoundaryFaces(), -1)
    {
        forAll(localToGlobal_, cellI)
        {
            localToGlobal_[cellI] = globalCells_.toGlobal(cellI);
        }

        if (Pstream::parRun())
        {
            const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

            forAll(boundaryMesh, patchI)
            {
                if (isA<processorPolyPatch>(boundaryMesh[patchI]))
                {
                    const processorPolyPatch& procPatch = 
                        refCast<const processorPolyPatch>(boundaryMesh[patchI]);
                    const labelUList& faceCells = procPatch.faceCells();
                    
                    labelList myGlobalIDsToSend(faceCells.size());
                    forAll(faceCells, i)
                    {
                        myGlobalIDsToSend[i] = localToGlobal_[faceCells[i]];
                    }
                    UOPstream toNeighbour(procPatch.neighbProcNo(), pBufs);
                    toNeighbour << myGlobalIDsToSend;
                }
            }

            pBufs.finishedSends();

            forAll(boundaryMesh, patchI)
            {
                if (isA<processorPolyPatch>(boundaryMesh[patchI]))
                {
                    const processorPolyPatch& procPatch = 
                        refCast<const processorPolyPatch>(boundaryMesh[patchI]);

                    UIPstream fromNeighbour(procPatch.neighbProcNo(), pBufs);
                    labelList remoteGlobalIDsReceived;
                    fromNeighbour >> remoteGlobalIDsReceived;

                    label startBndFace = procPatch.start() - mesh.nInternalFaces();
                    forAll(remoteGlobalIDsReceived, i)
                    {
                        bndFaceToRemoteGlobal_[startBndFace + i] = remoteGlobalIDsReceived[i];
                    }
                }
            }
        }
    }

    inline label global(const label localCellI) const { return localToGlobal_[localCellI]; }
    inline label remoteGlobal(const label bndFaceI) const { return bndFaceToRemoteGlobal_[bndFaceI]; }
    inline label totalCells() const { return globalCells_.totalSize(); }
};
}
#endif

```

## Flux Reconstruction and Segregated Turbulence Closures

A key advantage of using OpenFOAM as the physics engine is leveraging its library of advanced RANS and LES turbulence closures without modifying the underlying linear algebra infrastructure. During the initial setup phase inside the fields initialization header, the solver registers the fields and instantiates the chosen transport and turbulence closures.

Within the main execution loop, our architecture applies a physics co-design strategy that balances monolithic execution with segregated property updates through three distinct phases.

During the initial execution phase, the monolithic core tackles velocity and pressure coupling. The effective viscosity field, which combines laminar and turbulent eddy viscosity, is updated by the turbulence closure and fed into the effective stress tensor divergence. The spatial discretization derivatives are computed and passed to PETSc, which updates velocity and pressure simultaneously.

Following the resolution of the coupled system, the solver initiates a flux reconstruction phase. After unpacking the coupled solution vector, the face flux field must be explicitly reconstructed to enforce strict mass conservation across cell faces. This is achieved by applying a localized Rhie-Chow correction back onto the geometric field using the updated pressure gradient terms.

Finally, the framework handles the segregated turbulence equations. High-fidelity turbulence model transport equations, such as those for turbulent kinetic energy or specific dissipation rates, are highly non-linear but weakly coupled to the instantaneous kinematic fields. Solving them inside the monolithic block matrix would increase the system block size to six, exponentially increasing memory utilization and preconditioning costs. Therefore, the turbulence model equations are resolved in a segregated fashion via the turbulence correct method after the kinematic system has converged.

## Advanced Block Preconditioning and Schur Complement Approximations

Because the stabilized coupled matrix is highly indefinite, treating it with a standard monolithic Krylov solver leads to stagnation or divergence. We resolve this by applying a physics-based FieldSplit preconditioner. This operator algebraically splits the monolithic matrix into velocity and pressure spaces exclusively during the preconditioning phase.

The global system is transformed into an upper triangular block form by isolating the pressure space via the exact Schur complement matrix:

$$\mathbf{S} = \mathbf{A}_{pp} - \mathbf{A}_{pu} \mathbf{A}_{uu}^{-1} \mathbf{A}_{up}$$

Evaluating the exact inverse of the momentum block is computationally prohibitive due to its cubic complexity. Instead, we rely on sparse physical scaling operators, approximating the inverse using the inverse of the matrix diagonal block. This yields the scalable Least-Squares Commutator approximation of the Schur complement:

$$\hat{\mathbf{S}} \approx \mathbf{A}_{pp} - \mathbf{A}_{pu} \left[ \text{diag}(\mathbf{A}_{uu}) \right]^{-1} \mathbf{A}_{up}$$

Within this splitting framework, the velocity sub-block is isolated and accelerated via a flexible GMRES loop preconditioned with a single, aggressive V-cycle of Algebraic Multigrid to rapidly smooth high-frequency velocity errors. Conversely, the pressure sub-block operates under a Conjugate Gradient loop preconditioned via a dedicated multigrid cycle, resolving the elliptic pressure equation across distributed processors.

## Theoretical Hardware Profiling Strategy and L3 V-Cache Mapping

The ultimate theoretical optimization metric for this architecture involves its performance on modern split-cache processors that feature a stacked vertical cache layer. These platforms utilize a heterogeneous multi-die architecture where a massive static random-access memory block is stacked directly on top of a single core complex die, providing an expanded L3 cache pool accessible by a local group of high-performance cores.

Standard segregated solvers continuously cycle through entire cell data spaces, dropping the L3 cache hit rate below optimal levels and bottlenecking on the main memory channels due to continuous data streaming.

Our monolithic implementation counters this behavior through a strict alignment with the hardware cache layout. By utilizing the field splitting architecture, the dimensionality of the approximated Schur complement operator is limited exclusively to the pressure degree of freedom, reducing data requirements to one degree of freedom per cell instead of four.

For a targeted mesh partition size optimized for the cache limits, the structural footprint of this pressure sub-matrix remains small enough to reside entirely within the ultra-low-latency stacked L3 cache lines. This approach shifts the most demanding phase of the CFD calculation from a memory-bound streaming pattern into a compute-bound cache loop, effectively bypassing the memory bandwidth wall during multigrid cycles.

