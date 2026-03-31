# Changelog

## [Unreleased] — 2026-03-31

### Code Quality & Maintainability

#### 1. Refactor solidBodyDisplacementLaplacianZone

**Files:**
- `solidBodyDisplacementLaplacianZone/solidBodyDisplacementLaplacianZoneFvMotionSolver.C`
- `solidBodyDisplacementLaplacianZone/solidBodyDisplacementLaplacianZoneFvMotionSolver.H`

**Changes:**
- Extract duplicated constructor logic (~150 lines each) into `initZoneAndPointIDs()`
- Replace deprecated `lookup()` with `get<word>()` (3 occurrences)
- Fix debug block crash: `diffusivityPtr_()` when `nullptr` → read config string instead
- Add `holeCellThreshold_` named constant instead of magic number `0.5`
- Rebuild `pointIDs_` after topology change in `updateMesh()`
- Add explicit `pointMesh.H` include (was relying on transitive include)
- Rename `SBMFPtr_` → `solidBodyMotionPtr_` (follows naming convention)
- Add null check for `solidBodyMotionPtr_` after construction
- Fix `IOobject::MUST_READ` → `NO_READ` (contradicted `typeHeaderOk()` check)
- Use `dimLength` instead of `pointDisplacement_.dimensions()` (avoids UB in member-init)
- Add destructor cleanup for `diffusivityPtr_`

---

#### 2. Fix dynamicOversetZoneDisplacementFvMesh robustness

**Files:**
- `dynamicOversetZoneDisplacementFvMesh/dynamicOversetZoneDisplacementFvMesh.C`
- `dynamicOversetZoneDisplacementFvMesh/dynamicOversetZoneDisplacementFvMesh.H`

**Changes:**
- Fix stale header guard: `dynamicOversetMotionSolverFvMesh_H` → `dynamicOversetZoneDisplacementFvMesh_H`
- Add `IOstreams.H` include for `IOstreamOption` type
- Check `init()` return value and throw `FatalError` on failure
- Document `const_cast` safety for stencil `movePoints()` call
- Clarify critical ordering in `update()` (motion before overset update)
- Simplify `writeObject()` boolean logic (explicit `ok1 && ok2`)
- Document why destructor is empty (RAII handles cleanup)

---

#### 3. Fix fsiOmega field registration and logging

**File:** `fsiOmega/preciceOmega.C`

**Changes:**
- Register created field with Time registry via `checkIn()` (preCICE adapter can now find it)
- Guard `Info<<` output with `debug` flag (was flooding logs every timestep)
- Validate `fieldName_` not empty in `read()`
- Fix copy constructor: share field pointer (shallow copy) instead of creating disconnected instance

---

#### 4. Modernize precice-openfoam-adapter memory management

**Files:**
- `precice-openfoam-adapter/Adapter.C`, `Adapter.H`
- `precice-openfoam-adapter/Interface.C`, `Interface.H`
- `precice-openfoam-adapter/CouplingDataUser.C`, `CouplingDataUser.H`
- `precice-openfoam-adapter/FSI/Displacement.C`, `Displacement.H`
- `precice-openfoam-adapter/FSI/FSI.C`, `FSI.H`
- `precice-openfoam-adapter/FSI/ModuleFSI.C`
- `precice-openfoam-adapter/FSI/ForceBase.C`

**Changes:**
- Replace raw `new`/`delete` with `std::unique_ptr` for:
  - `meshPoints_`, `meshOldPoints_` in Adapter
  - `interpolationObjects_` in Displacement
  - `couplingDataReaders_`, `couplingDataWriters_` in Interface
- Simplify Interface destructor to `= default` (RAII handles cleanup)
- Pass strings by `const reference` to avoid copies in Interface and CouplingDataUser
- Use `vector::insert` instead of manual push_back loops for MPI gather buffers
- Remove stale TODO referencing deleted Stress class in ForceBase.C

---

#### 5. Remove unused FSI modules (Stress, DisplacementDelta)

**Deleted files:**
- `precice-openfoam-adapter/FSI/Stress.C`, `Stress.H`
- `precice-openfoam-adapter/FSI/DisplacementDelta.C`, `DisplacementDelta.H`

**Changes:**
- Remove `new Stress` / `new DisplacementDelta` from FSI.C writer/reader registration
- Remove includes from FSI.H
- Remove `#include` directives from ModuleFSI.C
- Update `docs/config.md` and `docs/README.md` to remove references to deleted data types

---

#### 6. Add AGENTS.md for coding agents

**File:** `AGENTS.md` (new)

**Content:**
- Build/clean commands for all plugins
- Testing approach (manual, no automated suite)
- OpenFOAM code style guidelines
- Naming conventions, types, error handling patterns
- MPI/parallel considerations
- Overset mesh specifics
- Known issues

---

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
- `dynamicOversetZoneDisplacementFvMesh/dynamicOversetZoneDisplacementFvMesh.C`
- `dynamicOversetZoneDisplacementFvMesh/dynamicOversetZoneDisplacementFvMesh.H`
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
- `dynamicOversetZoneDisplacementFvMesh/dynamicOversetZoneDisplacementFvMesh.C`

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

---

#### 5. Prevent Laplacian singularity on overset hole cells at time-window transitions

**File:** `solidBodyDisplacementLaplacianZone/solidBodyDisplacementLaplacianZoneFvMotionSolver.C`

**Problem:** Starting from time-window 2, the `smoothSolver` for the
`cellDisplacement_` Laplacian produced overflow values (~1.7e+16) in a

single iteration, despite starting from a zeroed internal field.

Root cause (two interacting issues):

1. **Stale initial guess from preCICE checkpoint.** The preCICE checkpoint
   is written immediately after `curPoints()` returns. At that point,
   `cellDisplacement_.primitiveField()` holds the Laplacian solution from
   the previous window (~8.6e-10 m). When the checkpoint is restored at
   the start of the next window, the field is restored to this non-zero
   state. The subsequent call to `cellDisplacement_.boundaryFieldRef().updateCoeffs()`
   used this non-zero internal field as donor values for the overset BC,
   producing a non-trivial initial guess for the Laplacian.

2. **Near-zero diagonal on hole cells.** Interior overset hole cells are
   not constrained by any patch BC. On a moved mesh, some hole cells are
   entirely surrounded by other hole cells, giving them a near-zero
   Laplacian diagonal. With a non-zero initial guess `x`, the Gauss-Seidel
   smoother computes the correction as `x / diag ≈ x / 0 = ∞`, producing
   overflow that propagates through the mesh motion to all mesh points.

**Fix (two parts):**

1. **Zero before `updateCoeffs`:** Reset `cellDisplacement_.primitiveFieldRef()`
   to zero *before* calling `updateCoeffs()`, so the overset patch always
   sees a zero donor field, matching the behavior of window 1.

2. **Pin hole cells in matrix:** After assembling `TEqn`, add `VGREAT` to
   the diagonal of all hole cells (`cellMask < 0.5`). This pins their
   solution to zero without affecting the physical cells. The `smoothSolver`
   then has a well-conditioned system regardless of mesh position.

3. **Zero `cellDisplacement_` after `curPoints()` consumes it:** Reset the
   field to zero at the end of `curPoints()` so the preCICE checkpoint
   always stores a clean zero field. This is a belt-and-suspenders measure
   on top of the zeroing in `solve()`.

---

#### 6. Zero `cellDisplacement_` after mesh motion in `curPoints()`

**File:** `solidBodyDisplacementLaplacianZone/solidBodyDisplacementLaplacianZoneFvMotionSolver.C`

**Problem:** `cellDisplacement_` is a scratch field whose only purpose is
to carry BC data from the preCICE adapter into the Laplacian solver and
then propagate it to `pointDisplacement_` via `interpolate()`. After
`curPoints()` returns, the field has no physical meaning. However, it was
left with the Laplacian solution from the current window in its internal
field (`primitiveField()`).

The preCICE adapter writes its checkpoint immediately after
`fvMesh::movePoints()` returns (which calls `curPoints()`). If
`cellDisplacement_` is non-zero at that moment, the checkpoint stores the
Laplacian solution. On checkpoint restore, this non-zero field is reloaded
into OpenFOAM, creating a spurious initial condition for the next
time-window's Laplacian solve.

**Fix:** Add `cellDisplacement_.primitiveFieldRef() = vector::zero` and
`cellDisplacement_.correctBoundaryConditions()` at the end of the
`else` branch in `curPoints()`, after the mesh points have been computed.
This ensures every checkpoint stores `cellDisplacement_ = 0`.

---

#### 7. Refactor preCICE adapter MPI logic into template helpers

**Files:**
- `precice-openfoam-adapter/Utilities.H`
- `precice-openfoam-adapter/Interface.C`

**Problem:** The manual point-to-point MPI communication logic introduced for v2506 compatibility was repetitive across multiple methods in `Interface.C`, making the code harder to maintain and read.

**Fix:** Abstracted the gather, scatter, and broadcast logic into template helper functions in `preciceAdapter` namespace within `Utilities.H`. This reduces code duplication and ensures a consistent implementation of the v2506-compatible MPI patterns.

---

#### 8. Standardize overset wrapper naming

**Files:**
- `dynamicOversetZoneDisplacementFvMesh/` (renamed from `dynamicOversetMotionSolverFvMesh/`)
- `README.md`

**Problem:** The overset wrapper directory and file names did not match the class name `dynamicOversetZoneDisplacementFvMesh`, leading to confusion.

**Fix:** Renamed the directory and source files to match the class name. Updated `Make/files` and `README.md` to reflect these changes.

---

#### 9. Modernize mathematical constants in `fsiOmega`

**File:** `fsiOmega/preciceOmega.C`

**Problem:** Usage of `M_PI` is not fully OpenFOAM-compliant and might rely on non-standard system headers.

**Fix:** Replaced `M_PI` with `Foam::constant::mathematical::pi` and included `mathematicalConstants.H`.

---

#### 10. Fix zonal motion propagation to neighbouring cells

**Files:**
- `solidBodyDisplacementLaplacianZone/solidBodyDisplacementLaplacianZoneFvMotionSolver.H`
- `solidBodyDisplacementLaplacianZone/solidBodyDisplacementLaplacianZoneFvMotionSolver.C`

**Problem:** In the zonal overset workflow, the rigid-body transform for the selected
`cellZone`/`cellSet` was applied directly to zone points, but the Laplacian system did
not explicitly constrain the corresponding zone cells to the same rigid displacement.
This could leave an inconsistent transition at the zone boundary and produce abrupt
deformation of neighbouring cells.

**Fix:**
1. Store selected zone `cellIDs_` and reference cell centres `cellCentres0_`.
2. In `solve()`, compute the rigid displacement of selected cell centres and pin those
   rows in `TEqn` using `VGREAT` diagonal/source contributions.
3. In `curPoints()`, subtract the rigid component on zone points when composing
   `newPoints + pointDisplacement_` to avoid double counting rigid motion.
4. Refresh `cellCentres0_` after topology changes in `updateMesh()`.

---

#### 11. Fix preCICE adapter build compatibility with OpenFOAM v2506

**Files:**
- `precice-openfoam-adapter/Interface.C`
- `precice-openfoam-adapter/Utilities.C`

**Problem:** The adapter failed to compile on OpenFOAM v2506 due to two API
compatibility issues introduced by MPI/data-container changes:

1. `Foam::List` construction from STL iterators used unsupported/private
   constructors.
2. `Pstream::broadcast` was called with `std::vector<double>`, which has no
   OpenFOAM stream operators.

Additionally, `wmkdepend` reported a parse warning in `Utilities.C` because
the file ended without a trailing newline.

**Fix:**
1. Replaced iterator-based `Foam::List` construction with explicit
   `setSize()` + element copy loops for rank-local gather/scatter buffers.
2. Reworked global-data broadcast paths to use temporary `Foam::List<double>`
   containers for MPI serialization and copied data back to `std::vector`.
3. Added trailing newline in `Utilities.C` to silence `wmkdepend` parse
   warnings.
