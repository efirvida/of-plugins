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
| Sampling | `writeInterval 1 s` (10 Δt) → 150 samples in the averaging window |

Grids (uniform `blockMesh`; the paper's 502x348x348 wall-resolved-LES reference
is declared in the config only and is never rendered as a run grid):

| Grid | Nx x Ny x Nz | ~cells |
|---|---|---|
| coarse | 153 x 80 x 80 | 0.98 M |
| medium | 115 x 151 x 151 | 2.62 M |
| reference (not generated) | 502 x 348 x 348 | 60.8 M |

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

## Closure adaptation and limitations

- The paper's dynamic SGS model (Eq. 25) is **not available in standard
  OpenFOAM**; the headline closure is **WALE** (`cubeRootVol` delta) and the
  adaptation is recorded here and in `config/case.yaml`.
- The nacelle is represented only by the actuator surface: no body-fitted mesh
  and no resolved boundary layer.
- The reference is the paper's 502x348x348 wall-resolved LES; it is not
  reproduced by this package.
- URANS k-ω SST is a cheaper smoke **fallback** with a weaker claim; it is not
  the headline.

## Layout

```
config/case.yaml          single source of truth (YAML)
tools/generate_case.py    renderer + non-destructive --check
case/                     committed generated skeleton (never edited by hand)
runs/                     rendered run directories (gitignored)
results/                  comparison output (gitignored)
```

## Case generation

```sh
python3 tools/generate_case.py            # render the committed case/
python3 tools/generate_case.py --check    # exit 1 on missing/stale, never writes
```

`--mesh coarse|medium`, `--solver wale|urans`, `--ranks N`, `--case-dir DIR`.
Rendered files carry a "Generated from config/case.yaml" banner; `case/` is the
only committed copy and hand edits are detected by `--check`.

## Attribution

Geometry and benchmark: Yang & Sotiropoulos, *A new class of actuator surface
models for wind turbines*, [arXiv:1702.02108v4](https://arxiv.org/abs/1702.02108v4)
(2018), Sec. 4.1 and Fig. 3.
