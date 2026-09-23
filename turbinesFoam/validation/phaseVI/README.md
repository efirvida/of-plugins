# NREL Phase VI uniform-inflow validation

Validation package for the NREL/NASA-Ames Phase VI rotor (10.058 m diameter,
standard 5.029 m tip, two blades, 3° pitch) coupled to three vendored
`turbinesFoam` models: the actuator line (ALM), the no-mesh actuator surface
(ASM) and the mesh-backed actuator surface (`fvOptions.ASM-MESH`, the imported
`phaseVI_blade` surface). The three variants differ only in the blade
subdictionary element/surface keys; see "Model variants and formulation" for
the formulation being validated, the sub-grid caveat and the pre-registered
hypothesis.

The case is a uniform-inflow free box (no atmospheric boundary layer, no Mann
turbulence, no external inflow preprocessing):

| Axis | Extent | Fallback (`--domain squat`) |
|---|---|---|
| x | −5D … 15D | −5D … 10D |
| y | ±4D | ±3D |
| z | hub ± 2.5D (hub 12.192 m) | hub ± 3D |

Boundary conditions: uniform velocity inlet, pressure outlet, symmetry far
field. Solver: `pimpleFoam` + `kOmegaSST` URANS, `adjustTimeStep off` with a
fixed time step (D/32 0.008 s, D/48 0.005 s, D/64 0.004 s) chosen so the blade
tip moves less than the hub-adjacent cell per step. The optional Stage 3
variant (`--solver iddes`) switches to LES `kOmegaSSTIDDES` with a fixed
0.0025 s step (see "IDDES variant"). The full design and
requirements live in `openspec/changes/phase-vi-validation/` of the
`of-plugins` repository; the data sources and their sha256 are recorded in the
per-directory `PROVENANCE.md` files.

## Layout

```
config/case.yaml          single source of truth (YAML)
tools/case_config.py      load/validate YAML; derived mesh/kinematics helpers
tools/generate_case.py    renderer + non-destructive --check
tools/element_data.py     blade/hub element rows, -(twist + pitch) convention
tools/stage_blade_stl.py  sha256-gated STL staging for asm-mesh run dirs
case/                     committed generated skeleton (ALM/ASM/ASM-MESH twins)
data/                     geometry, polars, experiment, s809 + PROVENANCE per dir
scripts/                  runners, mesh/stage tooling, comparison
  slurm/stage0.slurm      authorized development-queue job
  slurm/production.slurm  prepared-only long-queue array (do not submit)
  slurm/stage3.slurm      prepared-only IDDES + nChordwise array (do not submit)
  slurm/stage3-d64.slurm  prepared-only D/64 sensitivity array (do not submit)
  slurm/asm-mesh.slurm    prepared-only ASM-mesh D/32 gate + D/48 array (do not submit)
runs/                     rendered run directories (gitignored)
results/                  comparison output (gitignored)
```

## Environment

Load the repository OpenFOAM toolchain and set the wmake identity variables
(see the root `AGENTS.md`; the SDumont `openfoam/v2506_*` module does not
export them):

```sh
module load openfoam/v2506_openmpi-4.1.4_gnu
export WM_ARCH=linux64 WM_COMPILER=Gcc WM_COMPILE_OPTION=Opt
export WM_PRECISION_OPTION=DP WM_LABEL_SIZE=64
export WM_OPTIONS=linux64GccDPInt64Opt
. "$WM_PROJECT_DIR/etc/bashrc"   # the module exports Int32 paths; this puts the Int64 binaries on PATH
export FOAM_USER_LIBBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$WM_PROJECT_DIR/platforms/$WM_OPTIONS/lib:$WM_PROJECT_DIR/platforms/$WM_OPTIONS/lib/sys-openmpi:$LD_LIBRARY_PATH
```

`scripts/check_environment.sh` accepts OpenFOAM **v2506 and v2412**, checks
`blockMesh`/`topoSet`/`checkMesh`/`pimpleFoam`/`decomposePar`/`mpirun`,
`libturbinesFoam.so` in `FOAM_USER_LIBBIN`, and then runs
`generate_case.py --check` (exit 3 environment, exit 4 stale case).

## Case generation

```sh
python3 tools/generate_case.py            # render the committed case/
python3 tools/generate_case.py --check    # exit 1 on missing/stale, never writes
```

`--mesh coarse|fine|ultra`, `--speed <7|10|13|15|20|25>`,
`--domain long|squat`, `--sequence H|S`, `--profile production|smoke`,
`--solver urans|iddes`, `--n-chordwise N`, `--ranks N`, `--case-dir DIR`,
`--end-revs FLOAT`, `--start-from startTime|latestTime`.
Rendered files carry a "Generated from config/case.yaml" banner; `case/` is
the only committed copy and hand edits are detected by `--check`.

The measured per-speed TSR is selected through `tipSpeedRatio` (7 m/s → 5.408,
10 → 3.798, 13 → 2.920, 15 → 2.530, 20 → 1.898, 25 → 1.521; Sequence S 7 → 5.407).
The mesh ladder renders 18 blocks, 100 % hexahedra:

| Resolution | Cell size | Cells | Use |
|---|---|---|---|
| coarse (D/32) | 0.3143 m | 6,674,304 | Stage 0 / development |
| fine (D/48) | 0.2095 m | 22,525,776 | headline |
| ultra (D/64) | 0.1572 m | 53,394,432 | optional sensitivity |

The proposal's "≈40 M" D/64 estimate is superseded by 53,394,432 (the design
refines z like x/y); recorded in `design.md` §5 and `config/case.yaml`.

## Running

```sh
scripts/runPhaseVI.sh -m alm|asm|asm-mesh -u <speed> [-mesh coarse|fine|ultra]
                      [--domain long|squat] [-s H|S] [--solver urans|iddes]
                      [--nchordwise N] [--ranks N] [--stage0] [--restart]
                      [--run] [--submit]
```

The runner validates the (model, speed, mesh, sequence, solver) combination
through `tools/case_config.py`, renders the case into
`runs/<model>-U<speed>-<mesh>` (for example `runs/asm-U7-coarse-s0`),
installs the selected `fvOptions` twin as `system/fvOptions`, hardlinks the
shared `runs/mesh-<mesh>/constant/polyMesh`, and writes `run.json` (git
commit, OpenFOAM version, config sha256, variant, solver, `n_chordwise`,
`ranks`, start time). The committed `case/` skeleton is never modified.

- `--restart` sets `startFrom latestTime` when a written time exists (otherwise
  it warns and keeps `startTime`); the production and Stage 3 arrays use it for
  requeues.
- `--stage0` caps `endTime` at 0.25 revolutions (spec bound ≤ 0.3) so a 48-rank
  job fits the 20-minute `sequana_cpu_dev` partition.
- `--solver iddes` renders the Stage 3 LES `kOmegaSSTIDDES` variant (see
  below) and appends `-iddes` to the run id.
- `--nchordwise N` overrides the ASM chordwise strip count (ASM only; rejected
  for ALM) and appends `-ncN`; the configured default is 5.
- `--ranks N` overrides `decomposition.number_of_subdomains` in the rendered
  `decomposeParDict` and for `mpirun -np`; inside a Slurm allocation the runner
  fails (exit 2) unless `N == SLURM_NTASKS`.
- `--run` executes `decomposePar -force` and
  `mpirun -np <ranks> pimpleFoam -parallel` (inside a Slurm allocation;
  48 subdomains by default).
- `--submit` with `--stage0` submits `slurm/stage0.slurm` to the development
  queue. Without `--stage0` it requires `PHASEVI_LONG_QUEUE_AUTHORIZED=1`
  **and** a passing `results/U7-H/sign_gate.json`; otherwise it exits 5. It
  refuses to carry `--solver iddes`, `--nchordwise` or `--ranks` (Stage 3 has
  its own prepared arrays).
- `-m asm-mesh` selects the mesh-backed surface twin (`fvOptions.ASM-MESH`),
  produces the run id `asm-mesh-U<speed>-<mesh>`, and stages the committed
  `geometry/stl/phaseVI_blade.stl` into `constant/triSurface/` through
  `tools/stage_blade_stl.py` before the mesh link (sha256-checked against the
  geometry metadata; a missing source or a hash mismatch aborts with exit 3 and
  leaves no stale copy). The staged hash is recorded in `run.json` as
  `staged_stl_sha256`. `--nchordwise` is accepted (ASM family) and `--submit`
  is refused (exit 2) pointing at the prepared-only
  `scripts/slurm/asm-mesh.slurm`.

Exit codes: **2** unsupported input (including a `--ranks`/`SLURM_NTASKS`
mismatch), **3** environment/mesh/solver failure, **4** stale generated case,
**5** authorization or sign-gate failure.

### IDDES variant

`--solver iddes` switches `constant/turbulenceProperties` to
`simulationType LES` with `LESModel kOmegaSSTIDDES` and `delta IDDESDelta`
(`IDDESDeltaCoeffs { Cw 0.15; }`), the only LES delta `kOmegaSSTIDDES` accepts
in OpenFOAM.com v2506. The `iddes:` block of `config/case.yaml` also carries
the DES-appropriate `div(phi,U) Gauss linear` scheme, `wallDist { nRequired
true; }` (`IDDESDelta` reads the wall-normal vectors) and the fixed DES time
step **0.0025 s** on every mesh — still below the tip-displacement bound
(37.856 m/s × 0.0025 s = 0.0946 m < 0.1572 m, the D/64 hub-adjacent cell).
The tip-displacement assertion runs against the selected time step, so the
`--solver iddes` render is validated like the URANS one.

### Meshes and Stage 0

```sh
scripts/mesh.sh coarse [--domain long|squat]   # blockMesh + topoSet + checkMesh
scripts/mesh.sh fine
scripts/stage0.sh                              # print the stage plan only
scripts/stage0.sh --submit                     # sbatch the dev-queue job
scripts/slurm/stage0.slurm                     # 48 ranks, 20 min, sequana_cpu_dev
```

`mesh.sh` accepts a mesh only when `checkMesh` reports "Mesh OK", 100 %
hexahedra, non-orthogonality ≤ 1e-10 and the cell count inside the band from
`config/case.yaml` (6.60–6.80 M, 22.40–22.70 M, 53.0–53.8 M).

### Staged plan

| Stage | Content | Queue | Status |
|---|---|---|---|
| 0 | `checkMesh` D/32 + D/48; ALM+ASM 7 m/s D/32, ≤ 0.3 rev | `sequana_cpu_dev` | authorized |
| 1 | 7 m/s URANS ALM+ASM, D/32 first then D/48 | long | prepared |
| 2 | {10, 13, 15, 25} m/s ALM+ASM + Sequence S 7 m/s repeat | long | prepared |
| 3 | IDDES ALM+ASM D/32+D/48, ASM `nChordwise` {1, 3} D/48, D/64 ALM+ASM sensitivity | long | prepared |
| asm-mesh | Mesh-backed ASM 7 m/s: D/32 performance gate, then D/48 headline | long | prepared |

Stages 1–3 and the ASM-mesh array are prepared but **must not be submitted**
until the long-queue authorization is granted. 20 m/s stays renderable but is
deliberately not staged, so no blanket job array can pick it up.
`slurm/production.slurm` is a **prepared-only** array: one task per (model,
speed, mesh, sequence), 48 ranks, ≤ 96 h per task, restart from `latestTime`.
It refuses to run without `PHASEVI_LONG_QUEUE_AUTHORIZED=1`, and no script
submits it automatically.

### Stage 3 arrays (prepared only)

Two prepared-only arrays cover Stage 3; both refuse to run without
`PHASEVI_LONG_QUEUE_AUTHORIZED=1` and resolve the package from
`$SLURM_SUBMIT_DIR` (Slurm spools the script, so `BASH_SOURCE` is not usable —
same fix as `stage0.slurm`/`production.slurm`). Submit them from the package
directory with `sbatch scripts/slurm/stage3.slurm` (or `stage3-d64.slurm`).

| Job | Nodes / ranks | Wall time | Array |
|---|---|---|---|
| `slurm/stage3.slurm` | 4 / 192 | 96 h | 6 tasks: IDDES ALM+ASM on coarse and fine, plus ASM URANS `nChordwise` 1 and 3 on fine (7 m/s, `--restart`) |
| `slurm/stage3-d64.slurm` | 8 / 384 | 96 h | 2 tasks: ALM and ASM URANS on ultra (D/64) at 7 m/s, `--restart` so a task can be chained if the wall time is ever exhausted |

Every task passes `--ranks <ntasks>`, so the rendered `decomposeParDict` and
`mpirun -np` match the allocation (the runner fails otherwise). The multi-node
sizes are required by the mesh/step combination, not by throughput alone: the
D/48 IDDES runs and the D/64 URANS runs are roughly an order of magnitude
heavier than the 48-rank Stage 1 jobs. At 48 ranks a single D/64 run would need
about 194 h, so the D/64 array runs at 384 ranks (8×), which fits comfortably
in one 96 h task.

### ASM-mesh array (prepared only)

`slurm/asm-mesh.slurm` follows the `production.slurm` shape: 48 ranks on one
node, `--time=24:00:00` (the delta's at-most-24 h bound, stricter than the
untouched 96 h `production.slurm`), `--array=0-1`, both tasks `--restart --run`.
Task 0 is `asm-mesh:7:coarse:H` — the **D/32 performance measurement gate** —
and task 1 is `asm-mesh:7:fine:H`, the headline run. The script refuses to run
without `PHASEVI_LONG_QUEUE_AUTHORIZED=1` (exit 5) and resolves the package from
`$SLURM_SUBMIT_DIR`, so it is **prepared only**: the D/32 task is executed and
reviewed (proceed / harden the candidate query / restrict the campaign) before
any D/48 preparation, and no automated step of this change submits it. The
queued `phaseVI-prod`, `phaseVI-stage3` and `phaseVI-stage3-d64` arrays are
read-only baselines.

#### D/32 performance measurement gate (prepared, not run)

The measurement is **prepared but not executed** by this change: it runs only
under explicit authorization, and no `sbatch` is issued here. When authorized,
task 0 is submitted alone and its output is reviewed before any D/48 job is
prepared.

The per-`addSup` instrumentation added by W4 (`logDistribution`, default
`true`) is the measurement source: the master rank prints one line per
`distribute()`,

```
Blade surface distribution 'turbine.blade1.surface': nodes <N>, candidates <C>,
mean <C/N>, max <max>, seconds <s>
```

and appends the matching row to
`postProcessing/bladeSurface/<owner>.surface_distribution.csv`
(`time,nodes,candidates,mean_candidates,max_candidates,seconds`). `candidates`
is the number of candidate cell entries visited in the call — the bounded-query
evidence, far below the naive `nodes × N_local cells` scan — and `seconds` is
the `addSup` wall time. The counters are observational only: the distributed
load is unchanged.

**Gate decision procedure** (recorded before any D/48 preparation):

| Outcome | When | Action |
|---|---|---|
| **proceed** | the measured per-step cost is within the practical budget | keep the prepared D/48 task; no code change |
| **harden** | the measured cost exceeds the budget | harden the candidate query (W4.4: mesh-tree/stencil) and re-measure task 0 |
| **restrict** | the campaign scope must shrink | restrict the authorized mesh/speed ladder; the scope choice stays with the user |

The gate is a **performance** gate only. The sub-grid caveat stands (D/32 cells
0.314 m vs chord 0.218–0.744 m; no resolved chordwise-physics claim) and the
kernel/width confound framing is unchanged.

## Model variants and formulation

Three `fvOptions` twins are rendered from the same configuration and differ
only in the blade subdictionary element/surface keys:

| Variant | Blade keys | Model |
|---|---|---|
| `fvOptions.ALM` | `elementType actuatorLineElement;` | actuator line (baseline) |
| `fvOptions.ASM` | `elementType actuatorSurfaceElement; nChordwise 5;` | no-mesh actuator surface |
| `fvOptions.ASM-MESH` | the surface keys plus `surfaceGeometry "constant/triSurface/phaseVI_blade.stl";` (and `kernel gaussian;` only for the ablation) | mesh-backed actuator surface |

**Formulation.** The mesh-backed variant implements the paper's chord-line blade
ASM (Yang & Sotiropoulos, arXiv:1702.02108v4) **extended to an imported
surface**: the blade load is mapped from the line elements onto the nodes of the
committed `phaseVI_blade` triangulation (element→patch association by radial
station, patch-area share per element) and distributed over the background cells
with the selected kernel (Eq. 18). It is explicitly **not a literal equation
port**: the paper's surface construction is replaced by the imported wetted
mesh, and the surface does not sample per-node inflow or recompute BEM loads
(the element BEM chain, chord-averaged inflow, dynamic stall, added mass and end
effects are unchanged). The surface is therefore a **distribution-only** model.

**Sub-grid caveat.** At D/32–D/64 the background cells are 0.314/0.210/0.157 m
while the blade chord is 0.218–0.744 m, so the imported surface is **sub-grid**:
the three-way comparison tests **model form** (how the load is distributed),
not resolved chordwise physics.

**Kernel/width confound and the Gaussian ablation.** The paper-cosine surface
(`kernel cosine`, default; support `2.5·h_i`) and the no-mesh ASM's mesh-only
Gaussian width (`kernel gaussian`; `ε_i = 2·cbrt(V_i)·meshFactor`) differ in
both kernel and width. `kernel gaussian` on the mesh-backed variant is the
**kernel-matched ablation** that narrows this confound (it does not eliminate
it: the no-mesh ASM uses one `ε` per element while the surface uses one per
node). It is an ablation, not a model, and is rendered only when selected.

**Pre-registered hypothesis.** Before the campaign, the expected direction and
approximate size of the geometry effect on the spanwise load are recorded here,
so the comparison tests a stated hypothesis rather than a post-hoc "the mesh is
better" claim:

- The imported-surface load is expected to shift the spanwise `c_ref_n` toward
  the blade **root transition and tip** first, where the chord-line element
  representation is coarsest, with the mid-span least affected.
- The expected size of the effect is of the order of the spanwise band
  (max(0.15, 20 % of measured)) or smaller, and **not** a resolved-chordwise
  redistribution; no agreement better than the documented band is promised a
  priori.
- The geometry effect is separated from the kernel/width confound by the
  Gaussian ablation above.

**MEXICO naming resolution.** The committed blade is the Phase VI
`phaseVI_blade` component; one STL serves both identical Phase VI blades
(azimuth is runtime). The retired `blade0/1/2` reservation encoded an un-sourced
per-blade assumption and is not re-used: a future MEXICO blade is a separate
`mexico_blade` component owned by the `mexico-validation` change, not by this
package.

**Prepared-only.** All prepared stages and arrays, including the ASM-mesh array,
are never submitted by this change; see the staged plan above.

## Comparison

```sh
scripts/comparePhaseVI.py --alm-dir runs/alm-U7-coarse-s0 \
                          --asm-dir runs/asm-U7-coarse-s0 \
                          --asm-mesh-dir runs/asm-mesh-U7-coarse \
                          --speed 7 --sequence H --sign-gate
```

Reads `postProcessing/turbines/0/turbine.csv` and
`postProcessing/actuatorLineElements/0/*.csv`, merges them with the committed
experimental CSVs, and writes `turbine_comparison.csv`,
`spanwise_comparison.csv`, `metrics.json`, `report.txt` and (with
`--sign-gate`) `sign_gate.json` into `--out` (default
`results/U<speed>-<sequence>/`). It averages revolutions 4–12 after discarding
revolutions 0–4, flags power/thrust drift above 1 %, and interpolates
`c_ref_n`/`c_ref_t` to 30/47/63/80/95 % span.

With `--asm-mesh-dir` the comparison is **three-way**: it additionally reads the
window-filtered `postProcessing/bladeSurface/*.csv` station table (the
per-node `*_nodes.csv` output is not a station table), merges station rows by
station id, and maps each station's `root_dist` to r/R with the existing
element formula (no coefficient definition is redefined). Turbine metrics come
from the same `turbine.csv` (the surface moment is already in the AFTAL
torque). A missing or incomplete ASM-mesh directory fails loudly
(`EXIT_MISSING_INPUT`), never falling back to a two-model merge. `metrics.json`
records the staged STL hash (`staged_stl_sha256`) and the surface kernel
(`surface_kernel`, `cosine` default) for the fair-comparison check.

Metric definitions (F1-corrected, fixed against the NREL reports; printed in
every output):

```
q_dyn   = 1/2 * rho * A * Uinf^2          R = rotorRadius, A = pi R^2
Q       = ct * q_dyn * R                  shaft torque (ct = torque coefficient)
P       = cp * q_dyn * Uinf               cp = ct * TSR
T_blade = sum_blades cd_blade * q_dyn     blade-only thrust (headline)
T_rotor = cd * q_dyn                      hub included; labelled secondary
rho     = WTBARO / (287.058 * (WTATEMP + 273.15))   from the measured row
```

`--thrust-scope blade|rotor` selects the headline thrust (default blade) and
`--match-eaeroth-span` restricts the blade sum to `r/R ≥ 0.25` to match the
EAEROTH pressure integration. `--rho auto|<kg/m³>` selects the measured-row ideal-gas density (default) or a
user value.

Bands (fixed in code and reported with every result): turbine
power/torque/thrust ± 15 %; spanwise max(0.15, 20 % of measured). Claims are
limited to trend, stall onset and agreement within these bands; no sub-5 %
agreement is promised a priori.

Spanwise CM is **excluded**: the element CSV has no `cm` column (follow-up
change adds `momentCoefficient_` to the element writer). The measured CM
values are still written to `metrics.json` for reference.

Exit codes: **1** missing/incomplete input or failed sign gate, **2** averaging
window shorter than four revolutions without `--allow-short-window`, **3**
configured TSR does not match the run, **0** comparison produced (band and
drift flags are reported, not hidden).

The **sign gate** (`--sign-gate`, run on the Stage 0 D/32 7 m/s runs) requires
the run TSR to match the configured value, `cp`/`cd` > 0, and positive
`c_ref_n`/`c_ref_t` at all five stations. Production submission refuses without
a passing `results/U7-H/sign_gate.json`.

## Modelling limitations

- Sub-cell ASM chord strips: at every affordable mesh (D/32–D/64) the
  chordwise strips are sub-grid, so the comparison tests chord-averaged inflow
  and the ALM-vs-ASM load shift, not a chord-resolved surface.
- **Sub-grid imported surface**: the mesh-backed ASM is sub-grid too (cells
  0.314/0.210/0.157 m versus chord 0.218–0.744 m); the three-way comparison
  tests model form, not resolved chordwise physics.
- **Kernel/width confound**: the paper-cosine surface and the no-mesh Gaussian
  ASM differ in kernel and width; the `kernel gaussian` mesh-backed run is the
  kernel-matched ablation, not a model, and narrows but does not eliminate the
  confound.
- **Distribution-only**: the surface redistributes the element load; it does
  not sample per-node inflow or recompute BEM loads.
- URANS `kOmegaSST` cannot capture deep-stall unsteadiness or hysteresis; the
  separated high-speed points (13–25 m/s) are trend and stall-onset evidence
  only.
- The tower is not modelled (rotor + hub only); EAEROTH is a blade pressure
  integration and the measured loads include tower shadow.
- Blockage is not simulated (uniform free box versus the NASA-Ames test
  section); documented, not corrected.

## Data and tests

`data/geometry/` (TP-500-29955 Table A-1), `data/polars/` (TP-500-29955 A-3..A-8
and TP-442-7817 B1..B4), `data/experiment/` (WDH `wt_loads_statistics.xls`,
sheet `ldsmean`, 0° yaw) and `data/s809/` (Somers, NREL/SR-440-6918, Table 2)
each carry a `PROVENANCE.md` with DOI/URL, source sha256, sheet/table, row
selection, units and extraction date. No report PDF or workbook is committed.
The imported blade surface itself is the committed `phaseVI_blade` component in
`turbinesFoam/geometry/` (binary STL + §4.6 metadata + `PROVENANCE.md`), built
by the pure-Python loft from `data/geometry/phaseVI_blade.csv` and the committed
S809 coordinates; `runPhaseVI.sh -m asm-mesh` stages it, hash-checked, into each
run directory.

```sh
python3 -m pytest turbinesFoam/tests/test_phasevi_case.py \
                    turbinesFoam/tests/test_phasevi_data.py \
                    turbinesFoam/tests/test_phasevi_compare.py \
                    turbinesFoam/tests/test_blade_stage.py -q
python3 -m pytest turbinesFoam/tests/ -q   # solver-driven suites auto-skip without OpenFOAM
```

The Phase VI tests are pure Python (no OpenFOAM): the case/data contracts plus
the comparison tooling (F1 metric formulas, `root_dist`→r/R mapping, drift
flag, fail-loud exit codes, sign gate, three-way dry merge on synthetic
fixtures including the real `bladeSurface` CSV schema) and the STL staging
helper (hash match, mismatch/missing-source abort, run-directory layout). The
upstream turbinesFoam tutorial tests are skipped when `WM_PROJECT_VERSION` is
unset.

## Attribution and license

The case generator and tooling in this package are adapted from
[`mttbrbr/single-actuator-line`](https://github.com/mttbrbr/single-actuator-line)
(main `8284be8c...`), GPL-3.0-or-later. The adaptation retains that license:
this package is distributed under **GPL-3.0-or-later**, consistent with
`turbinesFoam` (see `turbinesFoam/LICENSE`).

Data: rotor geometry and S809 polars are derived from the public U.S. Government
reports NREL/TP-500-29955 (DOI 10.2172/15000240), NREL/TP-442-7817 and — for the
committed S809 blade coordinates (`data/s809/`) — Somers, NREL/SR-440-6918; the
measured load statistics come from the public NWTC WDH dataset
(DOI 10.21947/WDH-DAP/1910052). No report PDF or workbook is redistributed —
only derived numeric tables with per-directory `PROVENANCE.md` (source URL,
sha256, sheet/table, row selection, units, extraction date).

The mesh-backed surface model (`bladeSurfaceSampler` / `bladeSurfaceSource`) and
the `phaseVI_blade` geometry component are **fork extensions**, not part of
upstream `turbinesFoam`. This package is validated against the pinned `of-plugins`
commit `26a1f46d79ff1481bba8d7fe2866516c367de232` (*feat(turbinesFoam): add
actuator surface element type*), which carries the blade actuator-surface patch
(`elementType actuatorSurfaceElement;`); later S2 commits add the imported
surface distributor, the geometry component and the ASM-mesh variant.
