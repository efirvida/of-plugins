# Changelog

## [Unreleased] — 2025-03-27

### Bug Fixes (OpenFOAM v2506 compatibility)

#### 1. Replace deprecated `Pstream::scatterList`/`gatherList` with `OPstream`/`IPstream`

**File:** `precice-openfoam-adapter/Interface.C`

**Problem:** In OpenFOAM v2506, `Pstream::scatterList` and
`Pstream::gatherList` are deprecated and silently broken:

- `scatterList` returns empty buffers on non-master ranks, leaving
  `dataBuffer_` with stale values from the previous coupling iteration.
  These stale Force values get interpreted as Displacement, causing a
  ~34 million× amplification that results in **SIGFPE** (floating-point
  exception) when the mesh deforms.
- `gatherList` deadlocks after the second coupling read, because the
  deprecated implementation no longer matches the expected MPI
  communication pattern.

**Fix:** Replace all 3 `scatterList` and 4 `gatherList` calls with
explicit point-to-point communication using `OPstream`/`IPstream` with
`Pstream::commsTypes::scheduled`:

- **gather:** ranks → master via `OPstream`, master ← ranks via `IPstream`
- **scatter:** master → ranks via `OPstream`, ranks ← master via `IPstream`

**Affected functions:**
- `gatherRegisterScatterIDs()` — counts gather+scatter, vertex gather, vertex ID scatter
- `configureMesh()` — triangle ID gather
- `readCouplingData()` — data buffer scatter
- `writeCouplingData()` — data buffer gather

---

#### 2. Prevent MPI deadlock in `fvm::laplacian` with distributed cyclicAMI patches

**File:** `solidBodyDisplacementLaplacianZone/solidBodyDisplacementLaplacianZoneFvMotionSolver.C`

**Problem:** After mesh motion, cyclicAMI boundary values become stale.
When the Laplacian motion solver builds its matrix via `fvm::laplacian()`
with a `"limited corrected"` snGrad scheme, the non-orthogonal correction
internally calls `fvc::grad()`, which triggers AMI interpolation. If the
cyclicAMI boundary data is stale at that point, the AMI weight calculation
triggers inconsistent MPI communication patterns across ranks, causing an
**MPI deadlock** (all 20 ranks freeze inside `fvm::laplacian`).

This bug only manifests starting at the second preCICE time window, after
the first mesh motion has occurred.

**Fix:** Call `cellDisplacement_.correctBoundaryConditions()` immediately
before `fvm::laplacian()` in `solve()`. This ensures all boundary
conditions (including cyclicAMI) are evaluated synchronously across all
ranks before the matrix assembly begins.

---

#### 3. Prevent SIGSEGV on overset meshes during motion solver init

**Files:**
- `dynamicOversetMotionSolverFvMesh/dynamicOversetMotionSolverFvMesh.C`
- `dynamicOversetMotionSolverFvMesh/dynamicOversetMotionSolverFvMesh.H`
- `solidBodyDisplacementLaplacianZone/solidBodyDisplacementLaplacianZoneFvMotionSolver.C`

**Problem:** When `solidBodyDisplacementLaplacianZone` is used with the
overset mesh wrapper `dynamicOversetZoneDisplacementFvMesh`, the motion
solver's `inverseDistanceDiffusivity` constructor eagerly computes
`wallDist`. The wall-distance boundary evaluation triggers
`oversetFvPatchField::initEvaluate()`, which constructs the
`cellCellStencil` (via `cellCellStencilObject::New()`). This happens
during the `dynamicMotionSolverFvMesh` base-class constructor — before
`oversetFvMeshBase` is constructed and before the derived-class vtable
is in place — causing a **SIGSEGV** inside `inverseDistance::update()`
(`ListPolicy::deallocate<Map<long>>`). The crash is consistent and
always occurs on the same MPI rank for a given mesh decomposition.

**Fix (two parts):**
1. **Deferred init in `dynamicOversetZoneDisplacementFvMesh`:** Pass
   `doInit=false` to `dynamicMotionSolverFvMesh` and call `init()`
   after `oversetFvMeshBase` is constructed, ensuring the full vtable
   (`lduAddr`, `interfaces`, overset `solve`/`interpolate` overrides)
   is available when the motion solver triggers overset operations.
2. **Lazy diffusivity in `solidBodyDisplacementLaplacianZone`:** Defer
   `motionDiffusivity` creation from the member-initializer list to the
   first call to `diffusivity()` (lazy init already supported by the
   accessor). This avoids triggering `wallDist` → overset stencil
   construction during the motion solver constructor.

---

#### 4. Prevent heap corruption on overset meshes during motion solve

**Files:**
- `solidBodyDisplacementLaplacianZone/solidBodyDisplacementLaplacianZoneFvMotionSolver.C`
- `dynamicOversetMotionSolverFvMesh/dynamicOversetMotionSolverFvMesh.C`

**Problem:** After the SIGSEGV init fix (Bug 3), the overset simulation
started successfully but crashed on the 3rd time step with **"corrupted
double-linked list"** (glibc abort) inside `inverseDistance::update()`
during `fvMesh::movePoints()`.

Root cause: the motion solver's `solve()` method called
`cellDisplacement_.correctBoundaryConditions()`, which triggers
`oversetFvPatchField::initEvaluate()`. When `active_` is false (normal
operation), the overset patch performs field interpolation via the
stencil's `mapDistribute::distribute()`. This introduces MPI
communication through UCX that corrupts glibc heap metadata. The
corruption is detected later when `inverseDistance::update()` allocates
or frees memory.

The standard `displacementLaplacianFvMotionSolver` uses
`cellDisplacement_.boundaryFieldRef().updateCoeffs()` instead, which
does not trigger overset field interpolation.

**Fix:**
1. **Runtime overset detection in `solve()`:** Check at runtime whether
   the `cellDisplacement` field has `overset` patches. If so, use
   `updateCoeffs()` (matching the standard solver). Otherwise, use
   `correctBoundaryConditions()` (preserving the cyclicAMI deadlock
   fix from Bug 2).
2. **Remove redundant destructor call:** Removed the
   `oversetFvMeshBase::clearOut()` call from
   `dynamicOversetZoneDisplacementFvMesh`'s destructor (upstream
   `dynamicOversetFvMesh` uses an empty destructor; `lduPtr_` is
   automatically cleaned up by autoPtr).
