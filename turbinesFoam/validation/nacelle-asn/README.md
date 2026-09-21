# Periodic-nacelle actuator-surface validation

Validation package for the standalone `nacelleSurfaceSource` against the
periodic-nacelle benchmark of Yang & Sotiropoulos, *A new class of actuator
surface models for wind turbines*,
[arXiv:1702.02108v4](https://arxiv.org/abs/1702.02108v4) (2018), Sec. 4.1
("Flow over periodically placed nacelles").

The case is rotor-less: a hemisphere + cylinder nacelle at Re = 1000 in a
streamwise-periodic 30R x 20R x 20R cell, with the body represented **only** by
the actuator surface model (no body-fitted mesh) — the source distributes the
surface force onto the background grid.

## Setup

| Item | Value |
|---|---|
| Domain | 30R x 20R x 20R: `x` in [-10R, 20R], `y,z` in ±10R; nacelle downstream end at `x = 0`, nose at `x = -7R` |
| Boundaries | streamwise `cyclic` (inlet ↔ outlet); crosswise `symmetry` (free-slip) |
| Geometry | `../../geometry/stl/nacelle.stl` (pipeline artefact, 2616 triangles; paper target ~2652, −1.4 %) |
| Reynolds | Re = `U∞ R / ν` = 1000 (`U∞ = 1 m/s`, `R = 1 m`, `ν = 1e-3 m²/s`) |
| Solver | `pimpleFoam`, fixed `deltaT = 0.1 R/U`, `adjustTimeStep off` |
| Closure | **LES WALE** headline; URANS k-ω SST fallback (`--solver urans`) |
| Run length | wash-out `2 T_ft`, time-averaging `5 T_ft` (`T_ft = Lx/U∞ = 30 s`) |
| Sampling | `writeInterval 1 s` (10 Δt) → 150 writes in the averaging window |
| Profile probes | vertical lines through the wake at every metric station, probed **every** time step (`writeControl timeStep`) → the 0.1 s series in `postProcessing/profiles/` |

Grids (uniform `blockMesh`; the paper's 502x348x348 wall-resolved-LES reference
is declared in the config only and is never rendered as a run grid):

| Grid | Nx x Ny x Nz | ~cells |
|---|---|---|
| coarse | 153 x 80 x 80 | 0.98 M |
| medium | 115 x 151 x 151 | 2.62 M |
| reference (not generated) | 502 x 348 x 348 | 60.8 M |

**Grid spacings are derived from the cell counts**, not the other way round: the
counts are the paper's grids and the config declares only `cells`. Over the
30R x 20R x 20R domain they give `Δx = 30R/153 ≈ R/5.1` and
`Δy = Δz = 20R/80 = R/4` (coarse), `Δx ≈ R/3.8` and
`Δy = Δz ≈ R/7.6` (medium). Exploration-time notes quoting `Δx = R/2.5`
(coarse) / `R/3.75` (medium) are inconsistent with the rendered counts and are
superseded by `generate_case.py`'s `cell_size_R()`: no rendered file or script
asserts those Δ values.

## Profile sampling and stations

The rendered case adds a `probes` function object (`profileSamples`) that writes
the velocity time series on a vertical line through the wake at every metric
station (`1R, 3R, 5R, 7R` plus the `10R` stretch), at cell centres nearest the
declared line (at most half a cell away; the offset is recorded in
`metrics.json`). `scripts/compareNacelle.py` consumes that series to form the
time-averaged `⟨u⟩(z)` and the **resolved** TKE `k(z) = 1/2 ⟨u'_i u'_i⟩` over the
configured window `[2 T_ft, 7 T_ft]`.

Station reconciliation (`10R`): the paper's panels are at **odd** multiples of
`R` (1R..19R) — there is no `10R` panel. A declared station with an exact panel
(`1R/3R/5R/7R`) is compared against it; the `10R` stretch is **bracketed** by
the `9R` and `11R` panels and its acceptance uses the worse of the two. No `10R`
reference profile is invented, and `metrics.json` records the mapping and the
`±1R` offsets.

## Staged run plan

- **stage0** (`sequana_cpu_dev`): `blockMesh` + `checkMesh` + a short stability
  run — mesh generation and sanity now.
- **production** (`sequana_cpu_long`): full wash-out + averaging — **prepared
  only; never submitted before explicit authorization**.

Nothing is launched on the long queue in S1.

## Metrics

Time-averaged `⟨u⟩(z)` and `k(z)` vertical profiles at `x = 1R, 3R, 5R, 7R`
(+ a `10R` stretch) through the axis, and the drag coefficient
`CD = |F_drag| / (0.5 ρ U∞² πR²)` (frontal area `πR²`). `k` is the **resolved**
turbulent kinetic energy computed from the velocity time series (the case's `k`
field is only used by the URANS fallback).

Acceptance is per station and per grid: the **medium** grid must agree at all
stations; the **coarse** grid only at `1R` and the far wake (the paper reports
coarse deficits too large at 3R–7R). The paper's permeable-disk `CD = 0.48` is a
**comparison datum, not a target**.

The bands are package choices (the paper compares graphically) and are
configurable: RMS deviation over the profile `rms(⟨u⟩) ≤ 0.10 u/U∞` and
`rms(k) ≤ 0.02 k/U∞²` (`--band-u` / `--band-k`). A claimed station outside its
band exits non-zero; unclaimed stations are reported but do not fail the run.

## Running and comparing

```sh
# render / check the committed skeleton
python3 tools/generate_case.py --check

# prepare and run one variant (inside an allocation for --run)
scripts/runNacelle.sh --mesh coarse [--solver wale|urans] [--ranks N] --run
scripts/runNacelle.sh --mesh coarse --stage0 --run --ranks 48   # short stability

# staged Slurm jobs (see the scripts for the environment recipe)
sbatch scripts/slurm/stage0.slurm         # dev partition, authorized
# sbatch scripts/slurm/production.slurm   # long queue: PREPARED ONLY

# compare a finished run against the digitized reference
python3 scripts/compareNacelle.py --run-dir runs/nacelle-coarse --mesh coarse
```

`runNacelle.sh` refuses to run on a non-development partition without
`NACELLE_LONG_QUEUE_AUTHORIZED=1`, removes stale `log.*` before running (OpenFOAM
skips a step whose log already exists), and never submits the production job.
`compareNacelle.py` writes `metrics.json`, `profiles.csv` and `report.txt` under
`results/<run-id>/`.

## Closure adaptation and limitations

- The paper's dynamic SGS model (Eq. 25) is **not available in standard
  OpenFOAM**; the headline closure is **WALE** (`cubeRootVol` delta) and the
  adaptation is recorded here and in `config/case.yaml`.
- The nacelle is represented only by the actuator surface: no body-fitted mesh
  and no resolved boundary layer.
- The reference is the paper's 502x348x348 wall-resolved LES; it is not
  reproduced by this package.
- The probes sit at the cell centres nearest the declared station lines (at most
  half a cell away, recorded per station) so no probe lands on a cell face or a
  processor boundary; the `10R` stretch is compared against the `9R`/`11R`
  bracket, not an exact panel.
- URANS k-ω SST is a cheaper smoke **fallback** with a weaker claim; it is not
  the headline.

## Layout

```
config/case.yaml          single source of truth (YAML)
tools/generate_case.py    renderer + non-destructive --check
case/                     committed generated skeleton (never edited by hand)
scripts/compareNacelle.py compare a run against the digitized reference
scripts/runNacelle.sh     prepare / run / submit one variant (queue gates)
scripts/slurm/            staged jobs (stage0 executes; production prepared only)
data/reference/           digitized wall-resolved-LES profiles + PROVENANCE
runs/                     rendered run directories (gitignored)
results/                  comparison output (gitignored)
```

## Case generation

```sh
python3 tools/generate_case.py            # render the committed case/
python3 tools/generate_case.py --check    # exit 1 on missing/stale, never writes
```

`--mesh coarse|medium`, `--solver wale|urans`, `--ranks N`, `--case-dir DIR`,
`--end-time T` (the stage0 stability bound) and `--start-from latestTime`
(restart). Rendered files carry a "Generated from config/case.yaml" banner;
`case/` is the only committed copy and hand edits are detected by `--check`.

## Attribution

Geometry and benchmark: Yang & Sotiropoulos, *A new class of actuator surface
models for wind turbines*, [arXiv:1702.02108v4](https://arxiv.org/abs/1702.02108v4)
(2018), Sec. 4.1 and Fig. 3.
