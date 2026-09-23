#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Adapted from mttbrbr/single-actuator-line (main 8284be8c...), GPL-3.0-or-later; see README.md "Attribution and license".
"""Shared configuration, validation and kinematics helpers for the NREL
Phase VI uniform-inflow validation case.

The YAML file `config/case.yaml` is the single source of truth: every
generated file under `case/` (and every run directory) is rendered from it by
`tools/generate_case.py`. This module owns the derived quantities so that the
generator, the runners and the tests agree on the arithmetic.

Adapted from `mttbrbr/single-actuator-line` (main 8284be8c..., GPL-3.0-or-later).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config" / "case.yaml"
DEFAULT_CASE_DIR = ROOT / "case"
BLADE_CSV = ROOT / "data" / "geometry" / "phaseVI_blade.csv"
POLARS_DIR = ROOT / "data" / "polars"
EXPERIMENT_DIR = ROOT / "data" / "experiment"

# Rejected scaffold value: a fixed 0.03 s cannot satisfy the tip-displacement
# constraint at the measured rotor speed.
SCAFFOLD_DELTA_T = 0.03

# Solver variants: URANS `kOmegaSST` (headline) and the Stage 3 IDDES
# confirmation. The name selects the turbulence properties, the scheme
# overrides and the fixed time step.
SOLVER_CHOICES = ("urans", "iddes")

# Model variants: the actuator line, the no-mesh actuator surface and the
# mesh-backed actuator surface (ASM over the imported blade STL).
MODEL_CHOICES = ("alm", "asm", "asm-mesh")

# Blade root cutout station; everything inboard of it is modelled as cylinder.
ROOT_CUTOUT_RADIUS = 0.5083
ROOT_CYLINDER_END = 1.2575


def load_config(path: Path | str = DEFAULT_CONFIG) -> dict[str, Any]:
    """Load and validate the YAML configuration."""
    source = Path(path).resolve()
    with source.open(encoding="utf-8") as stream:
        data = yaml.safe_load(stream)
    if not isinstance(data, dict):
        raise ValueError(f"{source} does not contain a YAML mapping")
    validate_config(data)
    return data


def _speeds(cfg: dict[str, Any]) -> dict[str, dict[str, Any]]:
    speeds = cfg["kinematics"]["speeds"]
    return {str(int(key)): value for key, value in speeds.items()}


def validate_config(cfg: dict[str, Any]) -> None:
    """Validate the structural invariants of the case configuration."""
    turbine = cfg["turbine"]
    diameter = float(turbine["diameter"])
    radius = float(turbine["radius"])
    if not math.isclose(diameter, 2.0 * radius, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError("turbine.diameter must be exactly twice turbine.radius")
    if int(turbine["n_blades"]) != 2:
        raise ValueError("the Phase VI rotor has two blades")
    if float(turbine["pitch_deg"]) != 3.0:
        raise ValueError("the Phase VI reference pitch is 3 degrees")

    hub = float(turbine["hub_height"])
    for domain_key in ("domain_D", "domain_squat_D"):
        domain = cfg[domain_key]
        if domain["x"][0] >= domain["x"][1] or domain["y"][0] >= domain["y"][1]:
            raise ValueError(f"{domain_key} x/y extents must be increasing")
        if domain["z_offset"][0] >= domain["z_offset"][1]:
            raise ValueError(f"{domain_key} z_offset must be increasing")
        # The free box has no ground plane: z may extend below 0 m.
        # The rotor disk must stay clear of every far boundary.
        if float(domain["x"][0]) > -0.5 or float(domain["x"][1]) < 0.5:
            raise ValueError(f"{domain_key} x does not contain the rotor plane")
        if float(domain["y"][0]) > -0.5 or float(domain["y"][1]) < 0.5:
            raise ValueError(f"{domain_key} y does not contain the rotor disk")
        if (
            float(domain["z_offset"][0]) > -0.5
            or float(domain["z_offset"][1]) < 0.5
        ):
            raise ValueError(f"{domain_key} z does not contain the rotor disk")

    breaks = cfg["mesh"]["breaks_D"]
    for axis in ("x", "y", "z"):
        values = [float(value) for value in breaks[axis]]
        if any(b <= a for a, b in zip(values, values[1:])):
            raise ValueError(f"mesh.breaks_D.{axis} must be strictly increasing")

    for name, resolution in cfg["mesh"]["resolutions"].items():
        for axis in ("x", "y", "z"):
            spec = resolution[axis]
            n_breaks = len(breaks[axis])
            if len(spec["cells"]) != n_breaks - 1:
                raise ValueError(f"mesh.resolutions.{name}.{axis}.cells has wrong length")
            if len(spec["grading"]) != n_breaks - 1:
                raise ValueError(f"mesh.resolutions.{name}.{axis}.grading has wrong length")
            if any(int(value) <= 0 for value in spec["cells"]):
                raise ValueError(f"mesh.resolutions.{name}.{axis}.cells must be positive")
            if any(float(value) <= 0 for value in spec["grading"]):
                raise ValueError(f"mesh.resolutions.{name}.{axis}.grading must be positive")
        delta_t = float(resolution["delta_t"])
        if math.isclose(delta_t, SCAFFOLD_DELTA_T, rel_tol=0.0, abs_tol=1e-12):
            raise ValueError(
                "the scaffold 0.03 s time step is rejected: the tip displacement "
                "per step would exceed the hub-adjacent cell"
            )
        low, high = (float(value) for value in resolution["target_cells"])
        total = cell_count(cfg, name)
        if not low <= total <= high:
            raise ValueError(
                f"mesh.resolutions.{name} renders {total} cells, outside the "
                f"configured band [{low:.3g}, {high:.3g}]"
            )

    speeds = _speeds(cfg)
    if sorted(speeds) != ["10", "13", "15", "20", "25", "7"]:
        raise ValueError("kinematics.speeds must define exactly 7/10/13/15/20/25 m/s")
    for key, entry in speeds.items():
        for field in ("tsr", "rotor_speed_rpm", "row"):
            if field not in entry:
                raise ValueError(f"kinematics.speeds.{key} is missing '{field}'")
        if float(entry["tsr"]) <= 0.0:
            raise ValueError(f"kinematics.speeds.{key}.tsr must be positive")
        if float(entry["rotor_speed_rpm"]) <= 0.0:
            raise ValueError(
                f"kinematics.speeds.{key}.rotor_speed_rpm must be positive"
            )

    _validate_stages(cfg)

    solver = cfg["solver"]
    for field in ("end_revolutions", "discard_revolutions"):
        if int(solver[field]) <= 0:
            raise ValueError(f"solver.{field} must be positive")
    if int(solver["discard_revolutions"]) >= int(solver["end_revolutions"]):
        raise ValueError("solver.discard_revolutions must be smaller than end_revolutions")
    if str(solver["iddes_model"]) != "kOmegaSSTIDDES":
        raise ValueError("solver.iddes_model must be kOmegaSSTIDDES")
    if int(cfg["decomposition"]["number_of_subdomains"]) <= 0:
        raise ValueError("decomposition.number_of_subdomains must be positive")

    iddes = cfg.get("iddes")
    if not isinstance(iddes, dict):
        raise ValueError("iddes must define the Stage 3 IDDES variant")
    # kOmegaSSTIDDES::setDelta() (OpenFOAM.com v2506) aborts unless the LES
    # delta is an IDDESDelta-based model; cubeRootVol and friends are invalid.
    if str(iddes.get("delta")) != "IDDESDelta":
        raise ValueError(
            "iddes.delta must be IDDESDelta: kOmegaSSTIDDES rejects any other "
            "LESdelta model at construction time"
        )
    if float(iddes.get("delta_t", 0.0)) <= 0.0:
        raise ValueError("iddes.delta_t must be positive")
    if not isinstance(iddes.get("wall_dist_n_required"), bool):
        raise ValueError("iddes.wall_dist_n_required must be a boolean")
    if not str(iddes.get("div_phi_U", "")).strip():
        raise ValueError("iddes.div_phi_U must not be empty")


def check_solver(solver: str) -> str:
    """Validate and return a solver variant name."""
    if solver not in SOLVER_CHOICES:
        raise KeyError(
            f"unsupported solver {solver!r}; choose from {list(SOLVER_CHOICES)}"
        )
    return solver


def speed_keys(cfg: dict[str, Any]) -> list[str]:
    return sorted(_speeds(cfg), key=float)


def _validate_stages(cfg: dict[str, Any]) -> None:
    """Validate the staged run plan (F3): only Stage 0 may execute now."""
    stages = cfg.get("stages")
    if not isinstance(stages, dict) or not stages:
        raise ValueError("stages must define the staged run plan")
    executing = [name for name, stage in stages.items() if stage.get("executes")]
    if executing != ["stage0"]:
        raise ValueError("only stage0 may be marked as executing")
    if float(stages["stage0"]["max_revolutions"]) > 0.3:
        raise ValueError("stage0.max_revolutions must not exceed 0.3")
    for name, stage in stages.items():
        assert isinstance(stage, dict)
        speeds = [float(value) for value in stage.get("speeds", [])]
        if 20.0 in speeds:
            raise ValueError(
                f"stages.{name}: 20 m/s is renderable but must not be staged"
            )
        for mesh in stage.get("meshes", []):
            if mesh not in cfg["mesh"]["resolutions"]:
                raise ValueError(f"stages.{name}: unknown mesh {mesh!r}")


def stages(cfg: dict[str, Any]) -> dict[str, Any]:
    return cfg["stages"]


def tsr_from_rpm(rpm: float, speed: float, radius: float) -> float:
    """Measured tip-speed ratio from the row RPM: TSR = RPM*2*pi*R/(60*U)."""
    return float(rpm) * 2.0 * math.pi * float(radius) / (60.0 * float(speed))


def assert_tsr_matches_experiment(
    cfg: dict[str, Any],
    experiment_csv: Path | None = None,
    tolerance: float = 0.002,
) -> None:
    """Check every configured TSR against the committed measured rows.

    The configuration carries the frozen measured values; the committed
    experiment CSV carries the row RPM. This asserts they agree within the
    documented +-0.002 tolerance.
    """
    path = Path(experiment_csv) if experiment_csv is not None else (
        EXPERIMENT_DIR / "sequence_H_performance.csv"
    )
    if not path.exists():
        raise FileNotFoundError(f"experiment CSV not found: {path}")
    radius = float(cfg["turbine"]["radius"])
    measured: dict[str, float] = {}
    with path.open(encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            measured[f"{float(row['wind_speed_m_s']):g}"] = float(row["rpm"])
    for key in speed_keys(cfg):
        if key not in measured:
            raise ValueError(f"{path}: no measured row for {key} m/s")
        configured = float(select_speed(cfg, key)["tsr"])
        derived = tsr_from_rpm(measured[key], float(key), radius)
        if abs(configured - derived) > tolerance:
            raise ValueError(
                f"{key} m/s: configured TSR {configured} differs from the "
                f"measured row TSR {derived:.6f} by more than {tolerance}"
            )
    sequence_s = cfg.get("sequence_s")
    if sequence_s:
        speed = float(sequence_s["speed"])
        path_s = EXPERIMENT_DIR / "sequence_S_performance.csv"
        if path_s.exists():
            with path_s.open(encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            match = [
                row for row in rows
                if math.isclose(
                    float(row["wind_speed_m_s"]), speed, abs_tol=1e-9
                )
            ]
            if not match:
                raise ValueError(f"{path_s}: no measured row for {speed:g} m/s")
            derived = tsr_from_rpm(
                float(match[0]["rpm"]), speed, radius
            )
            configured = float(sequence_s["tsr"])
            if abs(configured - derived) > tolerance:
                raise ValueError(
                    f"Sequence S {speed:g} m/s: configured TSR {configured} "
                    f"differs from the measured row TSR {derived:.6f} by more "
                    f"than {tolerance}"
                )


def select_speed(
    cfg: dict[str, Any], speed: str | float | int, sequence: str = "H"
) -> dict[str, Any]:
    """Return the kinematics entry for one supported wind speed.

    `sequence="S"` selects the Sequence S repeat (7 m/s, s0700000); only a
    speed with a configured `sequence_s` row is accepted.
    """
    key = f"{float(speed):g}"
    if sequence.upper() == "S":
        row = cfg.get("sequence_s")
        if not row or not math.isclose(float(row["speed"]), float(speed)):
            raise KeyError(
                f"no Sequence S row for {speed!r} m/s; configured: "
                f"{cfg.get('sequence_s', {}).get('speed')}"
            )
        return {
            "speed": float(row["speed"]),
            "tsr": float(row["tsr"]),
            "rotor_speed_rpm": float(row["rotor_speed_rpm"]),
            "row": row["row"],
        }
    speeds = _speeds(cfg)
    if key not in speeds:
        raise KeyError(
            f"unsupported wind speed {speed!r}; configured: {speed_keys(cfg)}"
        )
    entry = dict(speeds[key])
    entry["speed"] = float(key)
    entry["tsr"] = float(entry["tsr"])
    entry["rotor_speed_rpm"] = float(entry["rotor_speed_rpm"])
    return entry


def cell_sizes(length: float, cells: int, grading: float) -> list[float]:
    """OpenFOAM `simpleGrading` cell sizes for one block direction.

    `simpleGrading` fixes the ratio of the last cell to the first cell; the
    sizes form a geometric progression of `cells` terms.
    """
    if cells < 1:
        raise ValueError("cells must be positive")
    if cells == 1:
        return [length]
    q = float(grading) ** (1.0 / (cells - 1))
    first = length * (q - 1.0) / (q**cells - 1.0)
    return [first * q**i for i in range(cells)]


def axis_extents_D(cfg: dict[str, Any], axis: str, domain: str = "long") -> list[float]:
    """Physical extents (m) of the mesh blocks for one axis.

    x and y breaks are absolute D coordinates; z breaks are hub offsets. The
    squat fallback rescales the block breaks linearly into the tighter domain
    so the 18-block topology is preserved.
    """
    turbine = cfg["turbine"]
    diameter = float(turbine["diameter"])
    hub = float(turbine["hub_height"])
    breaks = [float(value) for value in cfg["mesh"]["breaks_D"][axis]]

    if domain == "squat":
        target = cfg["domain_squat_D"]
    else:
        target = cfg["domain_D"]

    if axis == "z":
        low, high = (float(value) for value in target["z_offset"])
        offsets = breaks
        return [
            hub
            + diameter
            * (
                low
                + (offset - offsets[0])
                * (high - low)
                / (offsets[-1] - offsets[0])
            )
            for offset in offsets
        ]

    low, high = (float(value) for value in target[axis])
    return [
        diameter
        * (
            low
            + (value - breaks[0]) * (high - low) / (breaks[-1] - breaks[0])
        )
        for value in breaks
    ]


def cell_count(cfg: dict[str, Any], mesh: str, domain: str = "long") -> int:
    """Total hexahedral cell count of a rendered mesh."""
    del domain  # cell counts are identical for both domains
    resolution = _resolution(cfg, mesh)
    total = 1
    for axis in ("x", "y", "z"):
        total *= sum(int(value) for value in resolution[axis]["cells"])
    return total


def _resolution(cfg: dict[str, Any], mesh: str) -> dict[str, Any]:
    resolutions = cfg["mesh"]["resolutions"]
    if mesh not in resolutions:
        raise KeyError(
            f"unsupported mesh resolution {mesh!r}; choose from "
            f"{sorted(resolutions)}"
        )
    return resolutions[mesh]


def hub_cell_size(cfg: dict[str, Any], mesh: str) -> float:
    """The finest (hub-adjacent) z cell of a rendered mesh, in metres."""
    resolution = _resolution(cfg, mesh)
    extents = axis_extents_D(cfg, "z")
    z_cells = [int(value) for value in resolution["z"]["cells"]]
    z_grading = [float(value) for value in resolution["z"]["grading"]]
    hub_adjacent = []
    for index, cells in enumerate(z_cells):
        sizes = cell_sizes(extents[index + 1] - extents[index], cells, z_grading[index])
        # The hub split sits at the end of the lower block and the start of the
        # upper block; both cells are hub-adjacent.
        hub_adjacent.append(sizes[-1] if index == 0 else sizes[0])
    return max(hub_adjacent)


def tip_speed(cfg: dict[str, Any], speed: str | float | int) -> float:
    """Blade tip speed Omega*R at the measured tip-speed ratio, in m/s."""
    entry = select_speed(cfg, speed)
    return entry["tsr"] * entry["speed"]


def solver_delta_t(cfg: dict[str, Any], mesh: str, solver: str = "urans") -> float:
    """Fixed time step (s) of a solver variant.

    URANS uses the per-mesh `delta_t` (D/32 0.008, D/48 0.005, D/64 0.004);
    the IDDES variant uses its own fixed `iddes.delta_t` on every mesh.
    """
    check_solver(solver)
    if solver == "iddes":
        return float(cfg["iddes"]["delta_t"])
    return float(_resolution(cfg, mesh)["delta_t"])


def tip_displacement(
    cfg: dict[str, Any], speed: str | float | int, mesh: str, solver: str = "urans"
) -> float:
    """Tip displacement per time step, in metres."""
    return tip_speed(cfg, speed) * solver_delta_t(cfg, mesh, solver)


def assert_tip_constraint(
    cfg: dict[str, Any], speed: str | float | int, mesh: str, solver: str = "urans"
) -> None:
    """Fail when the tip displacement per step reaches the hub-adjacent cell."""
    displacement = tip_displacement(cfg, speed, mesh, solver)
    hub_cell = hub_cell_size(cfg, mesh)
    if displacement >= hub_cell:
        raise ValueError(
            f"Omega*R*dt = {displacement:.6g} m at {speed} m/s on the {mesh} mesh "
            f"is not below the hub-adjacent cell {hub_cell:.6g} m"
        )


def kinematics(
    cfg: dict[str, Any],
    speed: str | float | int,
    mesh: str,
    sequence: str = "H",
    solver: str = "urans",
) -> dict[str, float]:
    """All derived per-speed/per-mesh quantities used by the renderers."""
    check_solver(solver)
    entry = select_speed(cfg, speed, sequence)
    omega = entry["tsr"] * entry["speed"] / float(cfg["turbine"]["radius"])
    t_rev = 2.0 * math.pi / omega
    delta_t = solver_delta_t(cfg, mesh, solver)
    assert_tip_constraint(cfg, speed, mesh, solver)
    return {
        "speed": entry["speed"],
        "tsr": entry["tsr"],
        "rpm": entry["rotor_speed_rpm"],
        "omega": omega,
        "t_rev": t_rev,
        "delta_t": delta_t,
        "end_time": float(cfg["solver"]["end_revolutions"]) * t_rev,
        "write_interval": float(cfg["solver"]["write_interval_rev"]) * t_rev,
        "turbulence_intensity": float(cfg["inflow"]["turbulence_intensity"]),
        "mixing_length": float(cfg["inflow"]["mixing_length_D"])
        * float(cfg["turbine"]["diameter"]),
    }


def inflow_fields(cfg: dict[str, Any], speed: str | float | int) -> dict[str, float]:
    """k and omega for the chosen inflow turbulence intensity."""
    velocity = float(speed)
    intensity = float(cfg["inflow"]["turbulence_intensity"])
    length = float(cfg["inflow"]["mixing_length_D"]) * float(cfg["turbine"]["diameter"])
    k = 1.5 * (intensity * velocity) ** 2
    c_mu = 0.09
    omega = math.sqrt(k) / (c_mu**0.25 * length)
    return {"k": k, "omega": omega}


def read_blade(path: Path | None = None) -> list[dict[str, float]]:
    """Read the committed Table A-1 blade geometry (comments start with '#')."""
    source = Path(path) if path is not None else BLADE_CSV
    rows: list[dict[str, float]] = []
    with source.open(encoding="utf-8") as stream:
        reader = csv.reader(line for line in stream if not line.startswith("#"))
        for row in reader:
            if not row:
                continue
            if len(row) != 4:
                raise ValueError(f"{source}: expected 4 columns, got {len(row)}")
            rows.append(
                {
                    "radius": float(row[0]),
                    "chord": float(row[1]),
                    "twist": float(row[2]),
                    "mount": float(row[3]),
                }
            )
    if len(rows) < 2:
        raise ValueError(f"{source}: at least two radial stations are required")
    if any(b["radius"] <= a["radius"] for a, b in zip(rows, rows[1:])):
        raise ValueError(f"{source}: radial stations must be strictly increasing")
    return rows


def turbine_origin(cfg: dict[str, Any]) -> tuple[float, float, float]:
    """Rotor centre in mesh coordinates: (x_D, y_D) * D at the hub height."""
    turbine = cfg["turbine"]
    diameter = float(turbine["diameter"])
    position = turbine["position_D"]
    return (
        float(position[0]) * diameter,
        float(position[1]) * diameter,
        float(turbine["hub_height"]),
    )


def foam_header(object_name: str, klass: str = "dictionary") -> str:
    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Website:  www.openfoam.com                      |
|   \\\\  /    A nd           |                                                 |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {klass};
    object      {object_name};
}}
// Generated from config/case.yaml. Do not edit by hand.

"""


def foam_vector(values: list[float] | tuple[float, ...]) -> str:
    return "(" + " ".join(f"{float(value):.9g}" for value in values) + ")"


def render_cell_list(values: list[list[float]], indent: str = "            ") -> str:
    return "\n".join(
        indent + "(" + " ".join(f"{float(value):.9g}" for value in row) + ")"
        for row in values
    )


def _main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--select",
        nargs="+",
        metavar="SPEED MESH MODEL [SEQUENCE]",
        help="validate a (speed, mesh, model) triple and print its kinematics; "
             "an optional fourth value selects the sequence (H or S)",
    )
    parser.add_argument(
        "--solver",
        choices=SOLVER_CHOICES,
        default="urans",
        help="solver variant validated by --select (default: urans)",
    )
    parser.add_argument(
        "--target-cells",
        metavar="MESH",
        help="print the accepted cell-count band of a mesh as 'LOW HIGH'",
    )
    parser.add_argument(
        "--stage",
        metavar="NAME",
        help="print one stage of the staged run plan as JSON",
    )
    args = parser.parse_args(argv)
    try:
        cfg = load_config(args.config)
        if args.select:
            if len(args.select) == 3:
                speed, mesh, model = args.select
                sequence = "H"
            elif len(args.select) == 4:
                speed, mesh, model, sequence = args.select
            else:
                raise KeyError("--select takes SPEED MESH MODEL [SEQUENCE]")
            if model not in MODEL_CHOICES:
                raise KeyError(
                    f"unsupported model {model!r}; choose from "
                    f"{list(MODEL_CHOICES)}"
                )
            values = kinematics(cfg, speed, mesh, sequence, args.solver)
            print(json.dumps({"model": model, **values}, indent=2, sort_keys=True))
        elif args.target_cells:
            resolution = _resolution(cfg, args.target_cells)
            low, high = (int(value) for value in resolution["target_cells"])
            print(f"{low} {high}")
        elif args.stage:
            stage = stages(cfg).get(args.stage)
            if stage is None:
                raise KeyError(
                    f"unknown stage {args.stage!r}; configured: "
                    f"{sorted(stages(cfg))}"
                )
            print(json.dumps(stage, indent=2, sort_keys=True))
        return 0
    except (KeyError, ValueError) as exc:
        print(f"case configuration error: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(_main())
