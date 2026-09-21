# Changelog

## [Unreleased] — 2026-04-02

### New Features

#### 1. Add `boundaryDecay` motionDiffusivity manipulator for AMI meshes

**Files:**
- `solidBodyDisplacementLaplacianZone/boundaryDecayDiffusivity.H` (new)
- `solidBodyDisplacementLaplacianZone/boundaryDecayDiffusivity.C` (new)
- `solidBodyDisplacementLaplacianZone/Make/files` (updated)

**Problem:** When using `solidBodyDisplacementLaplacianZone` with AMI
(Arbitrary Mesh Interface) meshes, non-linear diffusivity functions
(e.g. `quadratic inverseDistance`) propagate mesh deformation all the way
to the AMI boundary.  Because the diffusivity is non-linear, adjacent
cells on opposite sides of the AMI interface deform by different
magnitudes, breaking interface connectivity.

**Solution:** New `motionDiffusivity` manipulator that multiplies any base
diffusivity by a smooth decay factor that transitions from 0 at the
selected boundary patches to 1 at a user-specified distance.  This
guarantees zero diffusivity — and therefore zero mesh deformation — at
the AMI interface.

The decay function is a quintic smooth-step (C² continuous):
`f(xi) = 6*xi^5 - 15*xi^4 + 10*xi^3`, where `xi = min(d / decayDistance, 1)`
and `d` is the cell distance to the nearest selected patch (computed via
`patchWave` / meshWave).  Both first and second derivatives are zero at
the endpoints, ensuring smooth diffusivity gradients near the decay boundary.

**Usage:**
```
diffusivity boundaryDecay 0.05 (AMI.*) quadratic inverseDistance (blade);
```

Where:
- `0.05` — decay distance (metres): distance from AMI where D ramps 0 → 1
- `(AMI.*)` — patch name regex matching the AMI boundaries
- `quadratic inverseDistance (blade)` — any existing motionDiffusivity chain

#### 2. Add `displacementDecay` post-solve clamping for AMI meshes

**Files:**
- `solidBodyDisplacementLaplacianZone/solidBodyDisplacementLaplacianZoneFvMotionSolver.C`
- `solidBodyDisplacementLaplacianZone/solidBodyDisplacementLaplacianZoneFvMotionSolver.H`

**Problem:** After the Laplacian solve, the cyclicAMI boundary condition
interpolates non-zero displacement from the coupled side back onto AMI
points during `correctBoundaryConditions()`, overriding the `boundaryDecay`
protection.  Residual displacement reaches the AMI interface, degrading
the AMI weights (minimum weights dropping to 0) and eventually causing
floating-point exceptions.

**Solution:** New optional `displacementDecay` sub-dictionary in
`solidBodyDisplacementLaplacianZoneCoeffs`:

```
displacementDecay
{
    distance  0.05;      // decay length (m)
    patches   (AMI.*);   // regex matching the AMI patches
}
```

The clamp is applied **after** `correctBoundaryConditions()` in
`curPoints()` and after the solve in `solve()`, so the cyclicAMI BC cannot
reintroduce non-zero displacement at the interface:
1. Points/cells **outside** the rotation zone: displacement set to exactly
   zero (no gradual decay), preventing any residual displacement from
   reaching the outer domain through the Laplacian solve
2. Points/cells **inside** the zone: quintic smooth-step (C²) decay over
   `distance`, measured from the selected patches with `patchWave`
3. AMI patch boundary fields are explicitly zeroed

The feature is skipped entirely when `moveAllCells` is active (no zone
selected).

#### 3. Vendor `turbinesFoam` as a new plugin

**Files:**
- `turbinesFoam/` (new — complete upstream tree, git history stripped)
- `Allwmake`, `Allclean` (updated)
- `README.md` (updated)

**Problem:** The actuator line library
[turbinesFoam/turbinesFoam](https://github.com/turbinesFoam/turbinesFoam)
was needed as a base for further research development, but the upstream
repository uses git submodules and its own git history.

**Solution:** Vendored the upstream codebase as an independent plugin
directory:
- Full upstream tree imported (`src/`, `tutorials/`, `tests/`) with its
  original GPL-3.0 license preserved
- Original git history removed (independent fork); Docker packaging
  (`Dockerfile` and Docker install instructions) dropped
- Upstream style checker `foamStyleCheck` removed (community vera++ tool,
  not an official OpenFOAM formatter) — the repository remains submodule-free
- Root `Allwmake` builds `libturbinesFoam.so`; root `Allclean` invokes
  the upstream-named `Allwclean` script

#### 4. Refocus README on FSI with OpenFOAM

**Files:**
- `README.md`

**Problem:** The README presented the repository as a generic plugin
collection.  The plugins actually form a complete FSI pipeline for
rotating machinery (coupling → mesh motion → rotation → loading).

**Solution:** Rewrote the introduction around the FSI workflow, adding a
"How the plugins fit together" section mapping each plugin to its role in
the coupled simulation, badges for OpenFOAM/preCICE/licence, and clearer
plugin descriptions.  Technical usage sections unchanged.

#### 5. Add blade actuator surface model (ASM) to turbinesFoam

**Files:**
- `turbinesFoam/src/fvOptions/actuatorLineSource/actuatorLineElement/actuatorLineElement.{H,C}` (updated)
- `turbinesFoam/src/fvOptions/actuatorLineSource/actuatorLineElement/actuatorSurfaceElement.{H,C}` (new)
- `turbinesFoam/src/fvOptions/actuatorLineSource/actuatorLineSource.C` (updated)
- `turbinesFoam/src/Make/files` (updated)
- `turbinesFoam/tutorials/axialFlowTurbineASM/` (new, incl. `compareALMvsASM.py`)
- `turbinesFoam/tests/test_asm.py`, `tests/test_aftal_asm.py`, `tests/axialFlowTurbineASMSource/` (new)
- `turbinesFoam/README.md`, `README.md` (updated)

**Problem:** The actuator line model (ALM) samples the inflow at a single
collocation point per radial station and ties the projection width to the
chord (`epsilon = max(0.25*c, mesh)`), so it cannot resolve chordwise flow
features and does not couple to meshes finer than the chord.

**Solution:** Optional blade actuator surface element
(`actuatorSurfaceElement`, Yang & Sotiropoulos, arXiv:1702.02108v4, Sec.
2.1), selected per line/blade with `elementType actuatorSurfaceElement;`
(+ optional `nChordwise`, default 5), that reuses the full inherited BEM
force chain but:
- averages the inflow over `nChordwise` equal chord strips (midpoint rule),
- distributes the force uniformly across the strips
  (`forceVector_/nChordwise` per strip, own Gaussian each), and
- uses a mesh-only projection width `2*cbrt(V)*meshFactor`.

The change completes the vestigial `actuatorLineElement` run-time
selection table (`New` + `elementType` default) so the ALM path stays
byte-identical when the key is absent.  Includes a HAWT ASM tutorial
(copy of `axialFlowTurbineAL` with the ASM blade config and truncated
`endTime`) and an ALM-vs-ASM comparison script, plus integration tests.
Documented as HAWT-only; CFTAL/VAWT usage untested.

---

#### 6. Add NREL Phase VI validation package

**Files:**
- `turbinesFoam/validation/phaseVI/` (new: YAML config, generator, committed
  `case/` skeleton, geometry/polars/experiment data with `PROVENANCE.md`,
  runners, Slurm scripts, README)
- `turbinesFoam/tests/test_phasevi_case.py`, `tests/test_phasevi_data.py`,
  `tests/test_phasevi_compare.py`, `tests/conftest.py` (new/updated)
- `README.md`, `CHANGELOG.md` (updated)

**Problem:** The fork's ALM and ASM models had no shared, reproducible rotor
validation case: earlier comparisons relied on scaffold data of unverified
provenance (flipped `Cm` signs, suspect `Cd`), there was no way to reproduce
the per-speed measured tip-speed ratio, no convergence/window contract, and no
gate against launching costly runs with a mirrored rotation/pitch convention.

**Fix:** Added the uniform-inflow NREL Phase VI package:
1. `config/case.yaml` single source of truth for the 18-block hex mesh
   (D/32 6,674,304 / D/48 22,525,776 / D/64 53,394,432 cells), fixed time step
   with a tip-displacement assertion, per-speed measured TSR and the ALM/ASM
   `fvOptions` twins; `tools/generate_case.py --check` detects stale renders
   without writing.
2. Committed geometry (TP-500-29955 Table A-1), polars (Tables A-3..A-8,
   TP-442-7817) and WDH experimental rows with per-directory `PROVENANCE.md`
   (DOI, source sha256, sheet/table, row rule, units, no PDF/XLS committed).
3. `scripts/runPhaseVI.sh` / `mesh.sh` / `stage0.sh` / `check_environment.sh`
   plus the Stage 0 development-queue job and the prepared-only 24 h
   production array (restart from `latestTime`, refuses to submit without the
   long-queue authorization and a passing 7 m/s sign gate).
4. `scripts/comparePhaseVI.py` merges `postProcessing/turbines/0/turbine.csv`
   and the element-level CSVs with the measured rows, averages revolutions
   4-12, flags > 1 % drift, interpolates `c_ref_n`/`c_ref_t` at
   30/47/63/80/95 % span, states the air density from `WTBARO`/`WTATEMP`,
   prints the fixed F1 metric definitions and the modelling limitations, and
   produces the sign gate that blocks production.
5. Pure-Python pytest contracts for the case, data and comparison tooling
   (cell counts, mesh blocks, twin diff, TSR, tip constraint, anchors,
   provenance, F1 metric formulas, fail-loud exit codes, sign gate); the
   OpenFOAM-driven upstream test modules are skipped automatically when no
   OpenFOAM environment is loaded.

---

#### 7. Add Stage 3 tooling: IDDES variant, `nChordwise`/ranks overrides, prepared arrays

**Files:**
- `turbinesFoam/validation/phaseVI/config/case.yaml`
- `turbinesFoam/validation/phaseVI/tools/case_config.py`
- `turbinesFoam/validation/phaseVI/tools/generate_case.py`
- `turbinesFoam/validation/phaseVI/scripts/runPhaseVI.sh`
- `turbinesFoam/validation/phaseVI/scripts/slurm/stage3.slurm` (new)
- `turbinesFoam/validation/phaseVI/scripts/slurm/stage3-d64.slurm` (new)
- `turbinesFoam/tests/test_phasevi_case.py`
- `turbinesFoam/validation/phaseVI/README.md`, `CHANGELOG.md` (updated)

**Problem:** The NREL Phase VI package declared the Stage 3 plan in
`config/case.yaml` (`iddes_model: kOmegaSSTIDDES`, the `stage3` block with
`n_chordwise: [1, 3, 5]` and the `ultra` mesh) but the tooling ignored it:
`generate_case.py` hardcoded `simulationType RAS; RASModel kOmegaSST;`, the
runner had no `nChordwise` or rank-count override, and Stage 3 had no Slurm
job — only the 48-rank single-node production array existed.

**Fix:**
1. `config/case.yaml` gains an `iddes:` block (`delta: IDDESDelta`,
   `delta_Cw: 0.15`, `delta_t: 0.0025`, `div_phi_U: Gauss linear`,
   `wall_dist_n_required: true`). `IDDESDelta` is mandatory: OpenFOAM.com
   v2506 `kOmegaSSTIDDES::setDelta()` aborts for `cubeRootVol` and any other
   `LESdelta` model (verified against the installed v2506 sources; the SAL
   scaffold's `delta IDDESDelta` confirms the syntax).
2. `generate_case.py` gains `--solver urans|iddes` (renders
   `simulationType LES; LESModel kOmegaSSTIDDES; delta IDDESDelta` plus the
   DES `deltaT 0.0025`, keeping the tip-displacement assertion),
   `--n-chordwise N` (ASM twin only) and `--ranks N`
   (`decomposeParDict numberOfSubdomains`). `case_config.py` validates the
   `iddes:` block, exposes `solver_delta_t()`/`check_solver()` and threads the
   solver through `kinematics()` and the `--select` CLI.
3. `runPhaseVI.sh` gains `--solver`, `--nchordwise` (rejected for ALM) and
   `--ranks`; run ids gain `-iddes`/`-nc<N>` suffixes, `run.json` records
   `solver`, `n_chordwise` and `ranks`, and inside Slurm a `--ranks` that
   differs from `SLURM_NTASKS` fails with exit 2. `--submit` refuses the new
   Stage 3 flags (Stage 3 has its own arrays).
4. Two prepared-only arrays, guarded by `PHASEVI_LONG_QUEUE_AUTHORIZED=1` and
   resolved from `$SLURM_SUBMIT_DIR` (never `BASH_SOURCE`):
   `slurm/stage3.slurm` (4 nodes / 192 ranks, 96 h, 6 tasks: IDDES ALM/ASM on
   coarse+fine and ASM `nChordwise` 1/3 on fine, all at 7 m/s) and
   `slurm/stage3-d64.slurm` (8 nodes / 384 ranks, 96 h, ALM/ASM ultra URANS at
   7 m/s, `--restart` to allow chaining). At 48 ranks one D/64 run would need
   ~194 h; at 384 ranks it fits in a single 96 h task.
5. Tests cover the IDDES render (`LESModel kOmegaSSTIDDES`, `deltaT 0.0025`,
   scheme/wallDist overrides, tip constraint per mesh), the `--n-chordwise 1`
   override, the `--ranks 192` render, the CLI end-to-end render into a temp
   directory, and the invalid-`delta` rejection. `--check` on the committed
   default case is unchanged and still passes.

---

### Bug Fixes (FSI Physics)

#### 1. Remove ghost time directories from implicit coupling sub-iterations

**Files:**
- `precice-openfoam-adapter/Adapter.C`
- `precice-openfoam-adapter/Adapter.H`

**Problem:** During implicit coupling, function objects (forces, propellerInfo,
etc.) write `functionObjectProperties` to the checkpoint time directory during
sub-iterations.  When the coupling window completes and time advances, these
directories remain as ghosts containing only
`uniform/functionObjects/functionObjectProperties` but no field data (U, p, …).
This breaks post-processing tools like `reconstructPar`, which expect all
time directories to contain field files.

**Fix:** After a coupling window completes, check whether the previous
checkpoint time directory is a ghost (no top-level field files) and remove
it with `rmDir()`.  Real write directories are left untouched.

#### 2. Fix FSI motion: apply deformation BEFORE rotation

**Files:**
- `solidBodyDisplacementLaplacianZone/solidBodyDisplacementLaplacianZoneFvMotionSolver.C`
- `solidBodyDisplacementLaplacianZone/solidBodyDisplacementLaplacianZoneFvMotionSolver.H`

**Problem:** The previous implementation computed:
```
P_final = T(P_0) - (T(P_0) - P_0) + U_lap = P_0 + U_lap
```
This cancelled the rigid rotation, effectively ignoring it for the deformation component.

**Physical Model:** For rotor FSI, the correct motion is:
```
U_total(x) = R(t)·(x + u_fsi(x)) - x
```
where the rotation must be applied AFTER the structural deformation.

**Fix:**
1. **curPoints():** Now computes `T(P_0 + U_deform)` — first deform, then rotate.
2. **solve():** Removed rigid displacement pinning. Added zone boundary constraint
   to prevent deformation diffusion outside the rotation zone.
3. **Documentation:** Updated class docs with the correct FSI motion equation.

#### 3. Upgrade `boundaryDecay` smooth-step to quintic (C²)

**File:** `solidBodyDisplacementLaplacianZone/boundaryDecayDiffusivity.H`

**Problem:** The cubic Hermite smooth-step (`3x²-2x³`) is only C¹ continuous —
the second derivative is discontinuous at ξ=0 and ξ=1 (jumps of ±6).  This
produces abrupt diffusivity gradients at the decay boundary, visible as a
sharp ring in the deformation field.  Can cause numerical stiffness in the
Laplacian mesh motion solver.

**Fix:** Replace with the quintic smooth-step (`6x⁵-15x⁴+10x³`), which has
`f'=0` and `f''=0` at both endpoints (C² continuous).  The transition is
significantly smoother, reducing the risk of mesh quality degradation.

#### 4. Fix zone boundary detection using face topology

**File:** `solidBodyDisplacementLaplacianZone/solidBodyDisplacementLaplacianZoneFvMotionSolver.C`

**Problem:** The zone boundary detection loop used `cellCells` indexing, which
does not account for boundary faces.  When a zone-boundary cell touches a
physical boundary (e.g. `propellerTip`), the `cellCells` array can be shorter
than the cell's face count, causing out-of-bounds access.

**Fix:** Replace `cellCells` lookup with explicit face owner/neighbour traversal.
Skip boundary faces (`!isInternalFace`), and use owner/neighbour to find the
correct adjacent cell.  Also removed a spurious `cellDisplacement_ = zero`
reset in `curPoints()` that was clearing the solved displacement field.

#### 5. Sync processor-boundary point positions in `curPoints()`

**File:** `solidBodyDisplacementLaplacianZone/solidBodyDisplacementLaplacianZoneFvMotionSolver.C`

**Problem:** After computing zone-rotated + elastically-deformed point
positions, processor-shared points can have tiny inconsistencies between
neighbouring ranks (~0.01–0.03% face-area mismatch), causing
`FOAM FATAL ERROR: face area does not match neighbour by ...%` on
processor patches during `calcGeometry()`.

**Fix:** When running in parallel, average the computed `curPoints` across
processor boundaries with `syncTools::syncPointList()` (sum plus count
normalisation) before returning them.

---

### Code Quality & Maintainability

#### 2. Refactor solidBodyDisplacementLaplacianZone

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

#### 3. Fix dynamicOversetZoneDisplacementFvMesh robustness

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

#### 4. Fix fsiOmega field registration and logging

**File:** `fsiOmega/preciceOmega.C`

**Changes:**
- Register created field with Time registry via `checkIn()` (preCICE adapter can now find it)
- Guard field discovery/creation `Info<<` messages with `debug` flag
- Keep omega/rpm logging in `value()`/`integrate()` unconditional
  (deliberate, for FSI debugging)
- Validate `fieldName_` not empty in `read()`
- Fix copy constructor: share field pointer (shallow copy) instead of creating disconnected instance

---

#### 5. Improve fsiDiagLog formatting and scalar data handling

**File:** `precice-openfoam-adapter/Interface.C`

**Changes:**
- Log scalar data (e.g. omega) as single value instead of meaningless norm/sum statistics
- Remove redundant `[FSI-DIAG]` prefix (already tagged by `adapterInfo`)
- Add braces to single-line control flow blocks (style consistency)
- Use `const` qualifiers for local variables

---

#### 6. Modernize precice-openfoam-adapter memory management

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

#### 7. Remove unused FSI modules (Stress, DisplacementDelta)

**Deleted files:**
- `precice-openfoam-adapter/FSI/Stress.C`, `Stress.H`
- `precice-openfoam-adapter/FSI/DisplacementDelta.C`, `DisplacementDelta.H`

**Changes:**
- Remove `new Stress` / `new DisplacementDelta` from FSI.C writer/reader registration
- Remove includes from FSI.H
- Remove `#include` directives from ModuleFSI.C
- Update `docs/config.md` and `docs/README.md` to remove references to deleted data types

---

#### 8. Add AGENTS.md for coding agents

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

#### 2. Replace deprecated `Pstream::scatterList`/`gatherList` with `OPstream`/`IPstream`

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

#### 3. Prevent MPI deadlock in `fvm::laplacian` with distributed cyclicAMI patches

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

#### 4. Prevent SIGSEGV on overset meshes during motion solver init

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

#### 5. Prevent heap corruption on overset meshes during motion solve

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

#### 6. Prevent Laplacian singularity on overset hole cells at time-window transitions

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

#### 7. Zero `cellDisplacement_` after mesh motion in `curPoints()`

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

#### 8. Refactor preCICE adapter MPI logic into template helpers

**Files:**
- `precice-openfoam-adapter/Utilities.H`
- `precice-openfoam-adapter/Interface.C`

**Problem:** The manual point-to-point MPI communication logic introduced for v2506 compatibility was repetitive across multiple methods in `Interface.C`, making the code harder to maintain and read.

**Fix:** Abstracted the gather, scatter, and broadcast logic into template helper functions in `preciceAdapter` namespace within `Utilities.H`. This reduces code duplication and ensures a consistent implementation of the v2506-compatible MPI patterns.

---

#### 9. Standardize overset wrapper naming

**Files:**
- `dynamicOversetZoneDisplacementFvMesh/` (renamed from `dynamicOversetMotionSolverFvMesh/`)
- `README.md`

**Problem:** The overset wrapper directory and file names did not match the class name `dynamicOversetZoneDisplacementFvMesh`, leading to confusion.

**Fix:** Renamed the directory and source files to match the class name. Updated `Make/files` and `README.md` to reflect these changes.

---

#### 10. Modernize mathematical constants in `fsiOmega`

**File:** `fsiOmega/preciceOmega.C`

**Problem:** Usage of `M_PI` is not fully OpenFOAM-compliant and might rely on non-standard system headers.

**Fix:** Replaced `M_PI` with `Foam::constant::mathematical::pi` and included `mathematicalConstants.H`.

---

#### 11. Fix zonal motion propagation to neighbouring cells

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

#### 12. Fix preCICE adapter build compatibility with OpenFOAM v2506

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

---

#### 13. Add nacelle/hub actuator surface model, geometry pipeline and validation package

**Files:**
- `turbinesFoam/src/fvOptions/nacelleSurface/nacelleSurfaceSampler.{H,C,I.H}` (new)
- `turbinesFoam/src/fvOptions/nacelleSurface/nacelleSurfaceSource.{H,C,I.H}` (new)
- `turbinesFoam/src/fvOptions/axialFlowTurbineALSource/axialFlowTurbineALSource.{H,C}` (updated)
- `turbinesFoam/src/fvOptions/turbineALSource/turbineALSource.{H,C}` (updated)
- `turbinesFoam/src/Make/files`, `turbinesFoam/src/Make/options` (updated)
- `turbinesFoam/geometry/` (new: `src/nacelle.geo`, `src/makeGeometry.py`,
  `stl/nacelle.stl`, `metadata/nacelle.json`, `README.md`, `PROVENANCE.md`)
- `turbinesFoam/validation/nacelle-asn/` (new: `config/case.yaml`,
  `tools/generate_case.py`, `case/`, `data/reference/`,
  `scripts/compareNacelle.py`, `scripts/runNacelle.sh`,
  `scripts/slurm/{stage0,production}.slurm`, `README.md`)
- `turbinesFoam/tests/{test_nacelle.py,test_nacelle_data.py,test_nacelle_case.py,test_nacelle_compare.py}` (new)
- `turbinesFoam/README.md`, `README.md` (updated)

**Problem:** `axialFlowTurbineALSource` declared a `nacelle {}` subdict and a
`nacelle_` `autoPtr` but `createNacelle()` was an empty stub, so any case with a
`nacelle {}` block dereferenced null in every `addSup` overload. There was no
nacelle force model, no deterministic geometry pipeline, and the rotor azimuth
was an incremental accumulator that silently reset on `startFrom latestTime`
restarts — an idempotency gap that blocks FSI.

**Fix:**
1. New `nacelleSurfaceSource` (a `cellSetOption` `fv::option`) owning the
   reusable `nacelleSurfaceSampler`: `triSurface` read, per-triangle centroid
   positions/normals/areas, the direct-forcing normal traction (Eq. 19), the
   Schultz–Grunow friction model with a constant `cf` override (Eqs. 21–23),
   the smoothed 4-point cosine kernel (Eqs. 7/8/18) for both interpolation and
   force distribution, MPI replicated-node/local-cell handling and CSV output.
   The path is additive: with no `nacelle {}` the ALM default is byte-identical.
   The source is usable standalone (the rotor-less validation case) or composed
   through the implemented `createNacelle()`, which closes the null-dereference.
2. Deterministic geometry pipeline: gmsh `.geo` + `makeGeometry.py` produce the
   committed binary nacelle STL (2616 triangles, paper target ~2652) plus
   metadata and provenance, with a non-destructive `--check` that verifies
   byte-identical regeneration in both the plain and OpenFOAM environments.
3. Time-derived, restart-safe kinematics: `angleDeg_` is an integral of the
   omega law (closed form for the `tsrAmplitude` oscillation) persisted as the
   `angleDeg.<name>` registry field, with an optional `omega.<name>` override
   for the future preCICE seam; the constant-TSR `angle_deg` CSV column is
   byte-identical to the accumulator.
4. Validation package `turbinesFoam/validation/nacelle-asn/`: the paper's
   periodic-nacelle benchmark (Re = 1000, 30R x 20R x 20R, cyclic/free-slip,
   WALE headline with a URANS k-omega SST fallback), coarse/medium grids
   rendered from a YAML single source of truth, metric-station profile probes,
   digitized wall-resolved-LES reference profiles with provenance, a compare
   tool (`⟨u⟩(z)`, resolved `k(z)`, `CD = |F_drag|/(0.5 ρ U∞² πR²)`; per-station
   and per-grid acceptance with the permeable-disk `CD = 0.48` kept as a datum)
   and staged Slurm jobs (stage0 dev-partition execution; production prepared
   only, never launched before authorization).

Reference: Yang, X. and Sotiropoulos, F., *A new class of actuator surface
models for wind turbines*, arXiv:1702.02108v4 (2018), Sec. 2.2 and 4.1.
