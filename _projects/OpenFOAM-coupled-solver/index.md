---
layout: project
title: "Monolithic Fluid Dynamics: Architecting High Performance Coupled Solvers via OpenFOAM and PETSc Integration"
date: 2026-05-30
tech: ["C++", "OpenFOAM", "PETSc", "Linear Algebra", "MPI"]
excerpt: "A comprehensive engineering analysis of monolithic CFD architectures leveraging the mathematical robustness of PETSc and the physical discretization of OpenFOAM to overcome memory bandwidth bottlenecks."
---

![Coupled Solver Concept](/assets/images/coupled-foam/hero_wide.png)

## Table of Contents
1. [Architectural Bottlenecks in Segregated Solvers](#architectural-bottlenecks-in-segregated-solvers)
2. [Mathematical Formulation of the Coupled System](#mathematical-formulation-of-the-coupled-system)
3. [Topology Mapping and Strict Preallocation](#topology-mapping-and-strict-preallocation)
4. [PETSc Integration and Distributed Matrix Assembly](#petsc-integration-and-distributed-matrix-assembly)
5. [Resolving the Parallel MPI Boundary Deadlock](#resolving-the-parallel-mpi-boundary-deadlock)
6. [Flux and Turbulence Integration](#flux-and-turbulence-integration)
7. [Advanced Preconditioning via the Schur Complement](#advanced-preconditioning-via-the-schur-complement)
8. [Performance Benchmarking](#performance-benchmarking)

---

## Architectural Bottlenecks in Segregated Solvers

The pursuit of high fidelity Computational Fluid Dynamics is strictly bounded by the physical limitations of modern hardware architectures, specifically the memory bandwidth wall. Standard segregated algorithms like SIMPLE fundamentally decouple the Navier Stokes equations, dividing the momentum and continuity constraints into isolated linear systems. This mathematical splitting introduces a severe numerical lag, requiring heavy under relaxation factors to prevent divergence in highly skewed meshes or complex flow regimes. 

The algorithmic necessity to continuously update and transfer intermediate velocity and pressure fields across the memory bus starves the floating point units of modern processors. The bottleneck is no longer raw compute power, but the latency involved in fetching fragmented matrix entries from main memory.

The monolithic approach directly addresses this inefficiency by resolving the velocity and pressure fields simultaneously within a single block matrix system. By embedding the physical interdependencies directly within the Jacobian matrix, the solver bypasses the iterative outer loops characteristic of segregated methods. The Navier Stokes equations are treated as a unified saddle point problem, capturing the exact physical coupling in one cohesive step.

## Mathematical Formulation of the Coupled System

To ground the monolithic architecture in rigorous fluid dynamics, we must observe the steady state incompressible Navier Stokes equations. The physical system is governed by the conservation of momentum and mass.

$$\nabla \cdot (\mathbf{u} \otimes \mathbf{u}) - \nabla \cdot (\nu \nabla \mathbf{u}) + \nabla p = 0$$

$$\nabla \cdot \mathbf{u} = 0$$

When discretized using the Finite Volume Method, these continuous operators are transformed into a massive system of linear algebraic equations. In the monolithic framework, these discrete operators are assembled into a unified saddle point block matrix topology.

$$\begin{bmatrix} \mathbf{A}_{uu} & \nabla \\ \nabla \cdot & \mathbf{0} \end{bmatrix} \begin{bmatrix} \mathbf{u} \\ p \end{bmatrix} = \begin{bmatrix} \mathbf{b}_u \\ \mathbf{0} \end{bmatrix}$$

The $\mathbf{A}_{uu}$ block contains the discretized convection and diffusion operators representing the momentum auto coupling. The off diagonal blocks represent the pressure gradient $\nabla$ and the velocity divergence $\nabla \cdot$. The zero block on the main diagonal strictly enforces the incompressible continuity constraint.

## Topology Mapping and Strict Preallocation

Memory preallocation is the most critical phase in distributed linear algebra. Dynamic memory reallocation during matrix assembly will catastrophically degrade performance. The topology mapping relies on computing the exact number of diagonal (`d_nnz`) and off diagonal (`o_nnz`) non zero block entries for each cell before matrix allocation. 

To ensure maximum optimization, we enforce a block structured allocation using the MATBAIJ format with a block size of 4, containing the interleaved degrees of freedom for velocity and pressure.

**File:** `PetscTopologyMapper.H`
```cpp
static PreallocData computeBlockPreallocation(const fvMesh& mesh)
{
    const label nCells = mesh.nCells();
    PreallocData allocData;
    allocData.d_nnz = labelList(nCells, 1);
    allocData.o_nnz = labelList(nCells, 0);

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
            forAll(faceCells, i) {
                allocData.o_nnz[faceCells[i]]++;
            }
        }
        else if (isA<cyclicPolyPatch>(patch))
        {
            const labelUList& faceCells = patch.faceCells();
            forAll(faceCells, i) {
                allocData.d_nnz[faceCells[i]]++;
            }
        }
    }
    return allocData;
}

```

## PETSc Integration and Distributed Matrix Assembly

OpenFOAM serves strictly as a physics host, providing the raw operator coefficients. Developing a monolithic solver requires bypassing the high level wrappers to directly access the Lower Diagonal Upper storage format.

The following C++ implementation demonstrates the distributed assembly loop, extracting the raw arrays from the OpenFOAM matrices and translating them into the coupled PETSc system, applying under relaxation implicitly directly into the RHS source terms.

**File:** `coupledSimpleFoam.C`

```cpp
Mat A;
MatCreate(PETSC_COMM_WORLD, &A);
MatSetSizes(A, localEqs, localEqs, globalEqs, globalEqs);
MatSetType(A, MATBAIJ); 
MatMPIBAIJSetPreallocation(A, bs, 0, d_nnz.data(), 0, o_nnz.data());
MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);

forAll(own, faceI)
{
    PetscInt global_o = indexer.global(own[faceI]);
    PetscInt global_n = indexer.global(nei[faceI]);
    vector nA = Sf[faceI];
    scalar w  = weights[faceI];

    vector Aup_upper = (1.0 - w) * nA; 
    vector Apu_upper = (1.0 - w) * nA; 

    PetscScalar b_upper[16] = {0.0};
    b_upper[0] = U_upper[faceI]; b_upper[5] = U_upper[faceI]; b_upper[10] = U_upper[faceI];
    b_upper[3] = Aup_upper.x();  b_upper[7] = Aup_upper.y();  b_upper[11] = Aup_upper.z();
    b_upper[12]= Apu_upper.x();  b_upper[13]= Apu_upper.y();  b_upper[14]= Apu_upper.z();
    b_upper[15]= p_upper[faceI];
    
    MatSetValuesBlocked(A, 1, &global_o, 1, &global_n, b_upper, ADD_VALUES);
}

// Source Terms and Under Relaxation Assembly
for (label cellI = 0; cellI < mesh.nCells(); ++cellI)
{
    PetscInt rows[4] = {global_c*4, global_c*4+1, global_c*4+2, global_c*4+3};
    vector U_old = U[cellI];
    scalar p_old = p[cellI];

    PetscScalar vals[4] = {
        U_source[cellI].x() + U_bndRHS[cellI].x() + ((1.0 - alpha_U) / alpha_U) * Auu_x * U_old.x(),
        U_source[cellI].y() + U_bndRHS[cellI].y() + ((1.0 - alpha_U) / alpha_U) * Auu_y * U_old.y(),
        U_source[cellI].z() + U_bndRHS[cellI].z() + ((1.0 - alpha_U) / alpha_U) * Auu_z * U_old.z(),
        p_source[cellI]     + p_bndRHS[cellI]     + ((1.0 - alpha_p) / alpha_p) * App * p_old
    };
    
    VecSetValues(b, 4, rows, vals, ADD_VALUES);
}

```

## Resolving the Parallel MPI Boundary Deadlock

Interfacing distributed memory structures between OpenFOAM and PETSc requires a reliable synchronization engine for halo cells. A naive sequential reading of processor boundaries introduces standard communication deadlocks, where neighboring MPI ranks stall waiting for mutual blocking operations.

To overcome this, the solver utilizes the OpenFOAM non blocking PstreamBuffers. The local global IDs are packed asynchronously into outbox streams, sent via non blocking primitives, and safely unpacked from incoming buffers to accurately assign remote global column indices without execution stalls.

**File:** `PetscIndexer.H`

```cpp
PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

forAll(boundaryMesh, patchI)
{
    if (isA<processorPolyPatch>(boundaryMesh[patchI]))
    {
        const processorPolyPatch& procPatch = refCast<const processorPolyPatch>(boundaryMesh[patchI]);
        labelList myGlobalIDsToSend(procPatch.faceCells().size());
        
        forAll(procPatch.faceCells(), i) {
            myGlobalIDsToSend[i] = localToGlobal_[procPatch.faceCells()[i]];
        }
        
        UOPstream toNeighbour(procPatch.neighbProcNo(), pBufs);
        toNeighbour << myGlobalIDsToSend;
    }
}

pBufs.finishedSends(); 

```

## Flux and Turbulence Integration

A critical advantage of using OpenFOAM as the physics host is the ability to leverage its massive library of turbulence models without modifying the core linear algebra. The solver natively instantiates the transport models and recalculates the face fluxes dynamically, ensuring mass conservation and accurate eddy viscosity integration within the coupled Jacobian.

**File:** `createFields.H`

```cpp
#include "createPhi.H"

mesh.setFluxRequired(p.name());

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

```

## Advanced Preconditioning via the Schur Complement

The viability of the coupled solver hinges on the FieldSplit preconditioner, which algebraically separates the monolithic matrix back into its physical constituents exclusively for the preconditioning phase. The algorithm isolates the pressure block by calculating the Schur complement of the velocity field.

Mathematically, the exact Schur complement $\mathbf{S}$ is defined by eliminating the velocity degrees of freedom from the block system.

$$\mathbf{S} = - (\nabla \cdot) \mathbf{A}_{uu}^{-1} \nabla$$

Because computing the exact inverse of the momentum block $\mathbf{A}_{uu}^{-1}$ is computationally prohibitive, we rely on sparse approximations of the Schur complement. This approximation provides a highly accurate representation of the pressure velocity coupling without requiring the full inversion.

The momentum block is preconditioned using an Algebraic Multigrid algorithm, which efficiently damps high frequency errors across the computational grid. The approximated Schur complement, representing the elliptic nature of the pressure field, is similarly treated with aggressive multigrid cycles.

## Performance Benchmarking

[Work in Progress]

Scientific transparency is the foundational pillar of this project. The rigorous profiling of the coupled architecture is currently under active development. To ensure absolute validity, no synthetic data or estimated speedups are presented at this stage. The validation campaign is strictly focused on isolating the exact memory bandwidth utilization and floating point efficiency of the preconditioning structures compared to traditional segregated methods.

The upcoming performance analysis will detail several key benchmarking phases. First, it will cover the execution of the standard motorbike test case scaled to 17 million hexahedral cells, directly comparing the OpenFOAM v2512 SIMPLE algorithm against the custom monolithic implementation. This will be followed by an evaluation of convergence metrics and failure states when utilizing standard incomplete LU factorization on the highly indefinite coupled block matrix.

Furthermore, the analysis will provide a rigorous profiling of the FieldSplit preconditioning strategy paired with the BoomerAMG solver and Schur complement approximation, detailing exact iteration counts and wall clock times. Finally, a hardware telemetry analysis on the AMD Ryzen 9 7900X3D will explicitly measure L3 cache hit rates. This phase will verify if the reduced dimensionality of the Schur complement matrix allows it to reside entirely within the 96MB 3D V-Cache of the primary core complex die, effectively bypassing the DDR5 memory bandwidth wall during aggressive Multigrid cycles.

