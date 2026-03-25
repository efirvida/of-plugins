# Changelog

## [Unreleased] — 2025-03-25

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
