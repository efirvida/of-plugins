# NREL Phase VI uniform-inflow validation

Validation package for the NREL/NASA-Ames Phase VI rotor (10.058 m diameter,
standard 5.029 m tip, two blades, 3° pitch) coupled to the vendored
`turbinesFoam` actuator line (ALM) and actuator surface (ASM) models.

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
tip moves less than the hub-adjacent cell per step. The full design and
requirements live in `openspec/changes/phase-vi-validation/` of the
`of-plugins` repository; the data sources and their sha256 are recorded in the
per-directory `PROVENANCE.md` files.

## Layout

```
config/case.yaml          single source of truth (YAML)
tools/case_config.py      load/validate YAML; derived mesh/kinematics helpers
tools/generate_case.py    renderer + non-destructive --check
tools/element_data.py     blade/hub element rows, -(twist + pitch) convention
case/                     committed generated skeleton (ALM/ASM twins)
data/                     geometry, polars, experiment + PROVENANCE per dir
scripts/                  runners, mesh/stage tooling, comparison
  slurm/stage0.slurm      authorized development-queue job
  slurm/production.slurm  prepared-only long-queue array (do not submit)
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
`--case-dir DIR`, `--end-revs FLOAT`, `--start-from startTime|latestTime`.
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
scripts/runPhaseVI.sh -m alm|asm -u <speed> [-mesh coarse|fine|ultra]
                      [--domain long|squat] [-s H|S] [--stage0] [--restart]
                      [--run] [--submit]
```

The runner validates the (model, speed, mesh, sequence) combination through
`tools/case_config.py`, renders the case into `runs/<model>-U<speed>-<mesh>`
(for example `runs/asm-U7-coarse-s0`), installs the selected `fvOptions` twin
as `system/fvOptions`, hardlinks the shared `runs/mesh-<mesh>/constant/polyMesh`,
and writes `run.json` (git commit, OpenFOAM version, config sha256, variant,
start time). The committed `case/` skeleton is never modified.

- `--restart` sets `startFrom latestTime` when a written time exists (otherwise
  it warns and keeps `startTime`); `production.slurm` uses it for requeues.
- `--stage0` caps `endTime` at 0.25 revolutions (spec bound ≤ 0.3) so a 48-rank
  job fits the 20-minute `sequana_cpu_dev` partition.
- `--run` executes `decomposePar -force` and
  `mpirun -np <decomposeParDict> pimpleFoam -parallel` (inside a Slurm
  allocation; 48 subdomains by default).
- `--submit` with `--stage0` submits `slurm/stage0.slurm` to the development
  queue. Without `--stage0` it requires `PHASEVI_LONG_QUEUE_AUTHORIZED=1`
  **and** a passing `results/U7-H/sign_gate.json`; otherwise it exits 5.

Exit codes: **2** unsupported input, **3** environment/mesh/solver failure,
**4** stale generated case, **5** authorization or sign-gate failure.

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
| 3 | optional IDDES, `nChordwise` {1, 3, 5} at D/48, D/64 sensitivity | long | prepared |

Stages 1–3 are prepared but **must not be submitted** until the long-queue
authorization is granted. 20 m/s stays renderable but is deliberately not
staged, so no blanket job array can pick it up. `slurm/production.slurm` is a
**prepared-only** array: one task per (model, speed, mesh, sequence), 48 ranks,
≤ 24 h per task, restart from `latestTime`. It refuses to run without
`PHASEVI_LONG_QUEUE_AUTHORIZED=1`, and no script submits it automatically.

## Comparison

```sh
scripts/comparePhaseVI.py --alm-dir runs/alm-U7-coarse-s0 \
                          --asm-dir runs/asm-U7-coarse-s0 \
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
- URANS `kOmegaSST` cannot capture deep-stall unsteadiness or hysteresis; the
  separated high-speed points (13–25 m/s) are trend and stall-onset evidence
  only.
- The tower is not modelled (rotor + hub only); EAEROTH is a blade pressure
  integration and the measured loads include tower shadow.
- Blockage is not simulated (uniform free box versus the NASA-Ames test
  section); documented, not corrected.

## Data and tests

`data/geometry/` (TP-500-29955 Table A-1), `data/polars/` (TP-500-29955 A-3..A-8
and TP-442-7817 B1..B4) and `data/experiment/` (WDH `wt_loads_statistics.xls`,
sheet `ldsmean`, 0° yaw) each carry a `PROVENANCE.md` with DOI/URL, source
sha256, sheet/table, row selection, units and extraction date. No report PDF
or workbook is committed.

```sh
python3 -m pytest turbinesFoam/tests/test_phasevi_case.py \
                    turbinesFoam/tests/test_phasevi_data.py \
                    turbinesFoam/tests/test_phasevi_compare.py -q
python3 -m pytest turbinesFoam/tests/ -q   # solver-driven suites auto-skip without OpenFOAM
```

The Phase VI tests are pure Python (no OpenFOAM): the case/data contracts plus
the comparison tooling (F1 metric formulas, `root_dist`→r/R mapping, drift
flag, fail-loud exit codes, sign gate, dry merge on synthetic fixtures). The
upstream turbinesFoam tutorial tests are skipped when `WM_PROJECT_VERSION` is
unset.

## Attribution and license

The case generator and tooling in this package are adapted from
[`mttbrbr/single-actuator-line`](https://github.com/mttbrbr/single-actuator-line)
(main `8284be8c...`), GPL-3.0-or-later. The adaptation retains that license:
this package is distributed under **GPL-3.0-or-later**, consistent with
`turbinesFoam` (see `turbinesFoam/LICENSE`).

Data: rotor geometry and S809 polars are derived from the public U.S. Government
reports NREL/TP-500-29955 (DOI 10.2172/15000240) and NREL/TP-442-7817; the
measured load statistics come from the public NWTC WDH dataset
(DOI 10.21947/WDH-DAP/1910052). No report PDF or workbook is redistributed —
only derived numeric tables with per-directory `PROVENANCE.md` (source URL,
sha256, sheet/table, row selection, units, extraction date).
