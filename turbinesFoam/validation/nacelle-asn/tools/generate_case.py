#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Render the periodic-nacelle actuator-surface validation case from `config/case.yaml`.

Every generated file carries a "Generated from config/case.yaml" banner and
`--check` re-renders in memory and fails when a file on disk is missing or
stale, without writing anything.

Usage:
    generate_case.py [--mesh coarse|medium] [--solver wale|urans] [--ranks N]
                     [--case-dir DIR] [--check]

The package mirrors `validation/phaseVI/`: `config/case.yaml` is the single
source of truth and `case/` is the committed skeleton rendered from it. The
nacelle is represented only by the standalone `nacelleSurfaceSource` (no
body-fitted mesh); the surface comes from the shared geometry pipeline
(`turbinesFoam/geometry/stl/nacelle.stl`).
"""

from __future__ import annotations

import argparse
import difflib
import math
import os
import sys
from pathlib import Path
from typing import Any

import yaml

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config" / "case.yaml"
DEFAULT_CASE_DIR = ROOT / "case"

# Run grids only: the paper's 502x348x348 wall-resolved-LES reference is
# declared in the config for provenance and must never be rendered as a run.
MESH_CHOICES = ("coarse", "medium")
SOLVER_CHOICES = ("wale", "urans")
FIELDS = ("U", "p", "k", "omega", "nut")
SIDE_PATCHES = ("bottom", "top", "sideMinus", "sidePlus")


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


def load_config(path: Path | str = DEFAULT_CONFIG) -> dict[str, Any]:
    """Load and validate the YAML configuration."""
    source = Path(path).resolve()
    with source.open(encoding="utf-8") as stream:
        data = yaml.safe_load(stream)
    if not isinstance(data, dict):
        raise ValueError(f"{source} does not contain a YAML mapping")
    validate_config(data)
    return data


def validate_config(cfg: dict[str, Any]) -> None:
    """Validate the structural invariants of the case configuration."""
    if not str(cfg["project"].get("name", "")).strip():
        raise ValueError("project.name must not be empty")

    nacelle = cfg["nacelle"]
    radius = float(nacelle["radius"])
    ratio = float(nacelle["cylinder_ratio"])
    if radius <= 0.0:
        raise ValueError("nacelle.radius must be positive")
    if ratio <= 0.0:
        raise ValueError("nacelle.cylinder_ratio must be positive")
    if int(nacelle["n_triangles"]) <= 0:
        raise ValueError("nacelle.n_triangles must be positive")
    stl = (ROOT / str(nacelle["stl"])).resolve()
    if not stl.is_file():
        raise ValueError(f"nacelle.stl does not resolve to a committed STL: {stl}")

    domain = cfg["domain_R"]
    spans: dict[str, float] = {}
    for axis in ("x", "y", "z"):
        low, high = (float(value) for value in domain[axis])
        if low >= high:
            raise ValueError(f"domain_R.{axis} extents must be increasing")
        spans[axis] = high - low
    if not math.isclose(spans["x"], 30.0 * radius, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError("domain_R.x must span 30R")
    for axis in ("y", "z"):
        if not math.isclose(spans[axis], 20.0 * radius, rel_tol=0.0, abs_tol=1e-9):
            raise ValueError(f"domain_R.{axis} must span 20R")

    nose, downstream_end = nacelle_bounds_R(cfg)
    # The nacelle downstream end sits at the streamwise origin (x = 0) and the
    # whole body must stay inside the periodic cell.
    if float(domain["x"][0]) >= nose or float(domain["x"][1]) <= downstream_end:
        raise ValueError(
            "domain_R.x must contain the whole nacelle, whose downstream end "
            "sits at the streamwise origin"
        )
    for axis in ("y", "z"):
        low, high = (float(value) for value in domain[axis])
        if low >= -radius or high <= radius:
            raise ValueError(f"domain_R.{axis} must contain the nacelle cross-section")

    boundaries = cfg["boundaries"]
    if str(boundaries["streamwise"]) != "cyclic":
        raise ValueError("boundaries.streamwise must be cyclic")
    if str(boundaries["crosswise"]) != "symmetry":
        raise ValueError("boundaries.crosswise must be symmetry (free-slip)")

    mesh = cfg["mesh"]
    if "reference" in mesh["resolutions"]:
        raise ValueError("the 502x348x348 reference grid must not be a run resolution")
    for name, resolution in mesh["resolutions"].items():
        cells = [int(value) for value in resolution["cells"]]
        if len(cells) != 3 or any(value <= 0 for value in cells):
            raise ValueError(f"mesh.resolutions.{name}.cells must be three positives")
    reference = mesh.get("reference")
    if not isinstance(reference, dict) or not reference.get("declared_only"):
        raise ValueError(
            "mesh.reference must be declared_only: the wall-resolved-LES "
            "resolution is never generated as a run grid"
        )

    box = mesh["cell_set_R"]
    for axis in ("x", "y", "z"):
        low, high = (float(value) for value in box[axis])
        if low >= high:
            raise ValueError(f"mesh.cell_set_R.{axis} extents must be increasing")
    if float(box["x"][0]) > nose or float(box["x"][1]) < downstream_end:
        raise ValueError("mesh.cell_set_R.x must cover the nacelle")
    for axis in ("y", "z"):
        low, high = (float(value) for value in box[axis])
        if low > -radius or high < radius:
            raise ValueError(f"mesh.cell_set_R.{axis} must cover the nacelle")

    inflow = cfg["inflow"]
    velocity = float(inflow["velocity"])
    nu = float(inflow["kinematic_viscosity"])
    if velocity <= 0.0 or nu <= 0.0:
        raise ValueError("inflow.velocity and inflow.kinematic_viscosity must be positive")
    if float(inflow["turbulence_intensity"]) < 0.0:
        raise ValueError("inflow.turbulence_intensity must not be negative")
    if float(inflow["mixing_length_R"]) <= 0.0:
        raise ValueError("inflow.mixing_length_R must be positive")

    solver = cfg["solver"]
    if str(solver["application"]) != "pimpleFoam":
        raise ValueError("solver.application must be pimpleFoam")
    if str(solver["headline"]) != "wale":
        raise ValueError("solver.headline must be wale (the documented adaptation)")
    if str(solver["fallback"]) != "urans":
        raise ValueError("solver.fallback must be urans (the weaker smoke path)")
    if float(solver["delta_t_R_U"]) <= 0.0:
        raise ValueError("solver.delta_t_R_U must be positive")
    if bool(solver["adjust_time_step"]):
        raise ValueError("the case runs a fixed time step (adjustTimeStep off)")

    run = cfg["run"]
    if not math.isclose(
        float(run["flow_through_time_s"]),
        spans["x"] / velocity,
        rel_tol=0.0,
        abs_tol=1e-9,
    ):
        raise ValueError("run.flow_through_time_s must equal Lx / U_inf")
    for key in ("wash_out_T", "average_T", "write_interval_s"):
        if float(run[key]) <= 0.0:
            raise ValueError(f"run.{key} must be positive")
    if int(run["purge_write"]) != 0:
        raise ValueError(
            "run.purge_write must stay 0: the averaging window keeps every sample"
        )

    metrics = cfg["metrics"]
    stations = [float(value) for value in metrics["stations_R"]]
    stretch = float(metrics["stretch_R"])
    if not stations or any(value <= 0.0 for value in stations):
        raise ValueError("metrics.stations_R must be positive")
    if stations != sorted(stations) or len(set(stations)) != len(stations):
        raise ValueError("metrics.stations_R must be strictly increasing")
    if stretch <= stations[-1]:
        raise ValueError("metrics.stretch_R must be downstream of every station")
    if not math.isclose(
        float(metrics["permeable_disk_datum"]), 0.48, rel_tol=0.0, abs_tol=1e-12
    ):
        raise ValueError("metrics.permeable_disk_datum must be the paper's 0.48 datum")

    acceptance = cfg["acceptance"]
    if str(acceptance["medium"]) != "all_stations":
        raise ValueError("acceptance.medium must be all_stations")
    coarse = [float(value) for value in acceptance["coarse"]]
    if stations[0] not in coarse or stretch not in coarse:
        raise ValueError("acceptance.coarse must claim 1R and the far wake")
    if any(station in coarse for station in stations[1:]):
        raise ValueError(
            "acceptance.coarse must not claim the 3R-7R stations (the paper "
            "reports coarse deficits too large there)"
        )

    _validate_stages(cfg)

    if int(cfg["decomposition"]["number_of_subdomains"]) <= 0:
        raise ValueError("decomposition.number_of_subdomains must be positive")


def _validate_stages(cfg: dict[str, Any]) -> None:
    """Validate the staged run plan: only stage0 may execute now."""
    stages = cfg.get("stages")
    if not isinstance(stages, dict) or "stage0" not in stages or "production" not in stages:
        raise ValueError("stages must define stage0 and production")
    executing = [name for name, stage in stages.items() if stage.get("executes")]
    if executing != ["stage0"]:
        raise ValueError("only stage0 may be marked as executing")
    if str(stages["stage0"]["queue"]) == str(stages["production"]["queue"]):
        raise ValueError("stage0 and production must target different queues")
    if stages["production"].get("executes"):
        raise ValueError("production must stay prepared-only (executes false)")


def check_solver(solver: str) -> str:
    """Validate and return a solver variant name."""
    if solver not in SOLVER_CHOICES:
        raise KeyError(
            f"unsupported solver {solver!r}; choose from {list(SOLVER_CHOICES)}"
        )
    return solver


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


def resolution(cfg: dict[str, Any], mesh: str) -> dict[str, Any]:
    resolutions = cfg["mesh"]["resolutions"]
    if mesh not in resolutions:
        raise KeyError(
            f"unsupported mesh resolution {mesh!r}; choose from {sorted(resolutions)}"
        )
    return resolutions[mesh]


def cell_counts(cfg: dict[str, Any], mesh: str) -> tuple[int, int, int]:
    """(Nx, Ny, Nz) of a run grid."""
    cells = [int(value) for value in resolution(cfg, mesh)["cells"]]
    return cells[0], cells[1], cells[2]


def cell_count(cfg: dict[str, Any], mesh: str) -> int:
    nx, ny, nz = cell_counts(cfg, mesh)
    return nx * ny * nz


def domain_extents_R(cfg: dict[str, Any]) -> dict[str, tuple[float, float]]:
    """Physical (R-scaled) extents of the periodic cell."""
    return {
        axis: (float(cfg["domain_R"][axis][0]), float(cfg["domain_R"][axis][1]))
        for axis in ("x", "y", "z")
    }


def cell_size_R(cfg: dict[str, Any], mesh: str, axis: str) -> float:
    """Uniform cell size in R units for one axis of a run grid."""
    low, high = domain_extents_R(cfg)[axis]
    counts = dict(zip(("x", "y", "z"), cell_counts(cfg, mesh)))
    return (high - low) / counts[axis]


def nacelle_bounds_R(cfg: dict[str, Any]) -> tuple[float, float]:
    """(nose, downstream end) streamwise positions of the nacelle, in R."""
    radius = float(cfg["nacelle"]["radius"])
    ratio = float(cfg["nacelle"]["cylinder_ratio"])
    return -radius * (1.0 + ratio), 0.0


def cell_set_box(cfg: dict[str, Any]) -> tuple[tuple[float, ...], tuple[float, ...]]:
    box = cfg["mesh"]["cell_set_R"]
    minimum = tuple(float(box[axis][0]) for axis in ("x", "y", "z"))
    maximum = tuple(float(box[axis][1]) for axis in ("x", "y", "z"))
    return minimum, maximum


def delta_t(cfg: dict[str, Any]) -> float:
    """Fixed time step (s): 0.1 R/U."""
    return (
        float(cfg["solver"]["delta_t_R_U"])
        * float(cfg["nacelle"]["radius"])
        / float(cfg["inflow"]["velocity"])
    )


def flow_through_time(cfg: dict[str, Any]) -> float:
    """T_ft = Lx / U_inf (s)."""
    low, high = domain_extents_R(cfg)["x"]
    return (high - low) / float(cfg["inflow"]["velocity"])


def run_times(cfg: dict[str, Any]) -> dict[str, float]:
    """Derived run schedule: wash-out, averaging window and write cadence."""
    run = cfg["run"]
    t_ft = flow_through_time(cfg)
    wash_out = float(run["wash_out_T"]) * t_ft
    average = float(run["average_T"]) * t_ft
    write_interval = float(run["write_interval_s"])
    return {
        "t_ft": t_ft,
        "wash_out_end": wash_out,
        "average_window": average,
        "end_time": wash_out + average,
        "delta_t": delta_t(cfg),
        "write_interval": write_interval,
        "samples": int(round(average / write_interval)),
    }


def inflow_fields(cfg: dict[str, Any]) -> dict[str, float]:
    """k and omega for the URANS fallback inflow (ignored by WALE)."""
    velocity = float(cfg["inflow"]["velocity"])
    intensity = float(cfg["inflow"]["turbulence_intensity"])
    length = float(cfg["inflow"]["mixing_length_R"]) * float(cfg["nacelle"]["radius"])
    k = 1.5 * (intensity * velocity) ** 2
    c_mu = 0.09
    if k <= 0.0:
        return {"k": 0.0, "omega": 0.0}
    return {"k": k, "omega": math.sqrt(k) / (c_mu**0.25 * length)}


def geometry_stl(cfg: dict[str, Any]) -> Path:
    return (ROOT / str(cfg["nacelle"]["stl"])).resolve()


def geometry_path(cfg: dict[str, Any], case_dir: Path) -> str:
    """The `geometry` entry as seen from the rendered case (or run) directory."""
    return os.path.relpath(geometry_stl(cfg), case_dir)


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


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


def render_block_mesh(cfg: dict[str, Any], mesh: str) -> str:
    """Single-block uniform blockMeshDict: cyclic streamwise, symmetry sides."""
    extents = domain_extents_R(cfg)
    nx, ny, nz = cell_counts(cfg, mesh)
    x0, x1 = extents["x"]
    y0, y1 = extents["y"]
    z0, z1 = extents["z"]
    vertices = [
        (x0, y0, z0),
        (x1, y0, z0),
        (x1, y1, z0),
        (x0, y1, z0),
        (x0, y0, z1),
        (x1, y0, z1),
        (x1, y1, z1),
        (x0, y1, z1),
    ]
    patches = [
        ("inlet", "cyclic", "neighbourPatch outlet;", "(0 4 7 3)"),
        ("outlet", "cyclic", "neighbourPatch inlet;", "(1 2 6 5)"),
        ("bottom", "symmetry", None, "(0 3 2 1)"),
        ("top", "symmetry", None, "(4 5 6 7)"),
        ("sideMinus", "symmetry", None, "(0 1 5 4)"),
        ("sidePlus", "symmetry", None, "(3 7 6 2)"),
    ]

    def patch(name: str, patch_type: str, extra: str | None, face: str) -> str:
        block = f"    {name}\n    {{\n        type {patch_type};\n"
        if extra:
            block += f"        {extra}\n"
        return block + f"        faces\n        (\n            {face}\n        );\n    }}"

    return (
        foam_header("blockMeshDict")
        + "scale 1;\n\nvertices\n(\n"
        + "\n".join(f"    {foam_vector(vertex)}" for vertex in vertices)
        + "\n);\n\nblocks\n(\n"
        + f"    hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz}) simpleGrading (1 1 1)\n"
        + ");\n\nedges ();\n\nboundary\n(\n"
        + "\n".join(
            patch(name, patch_type, extra, face)
            for name, patch_type, extra, face in patches
        )
        + "\n);\n\nmergePatchPairs ();\n"
    )


def render_toposet(cfg: dict[str, Any]) -> str:
    """One cellSet box enclosing the nacelle plus the kernel margin."""
    minimum, maximum = cell_set_box(cfg)
    return foam_header("topoSetDict") + f"""actions
(
    {{
        name nacelleCells;
        type cellSet;
        action new;
        source boxToCell;
        box {foam_vector(minimum)} {foam_vector(maximum)};
    }}
);
"""


def render_control_dict(cfg: dict[str, Any]) -> str:
    times = run_times(cfg)
    return foam_header("controlDict") + f"""application {cfg['solver']['application']};
startFrom startTime;
startTime 0;
stopAt endTime;
endTime {times['end_time']:.8g};
deltaT {times['delta_t']:.8g};
writeControl runTime;
writeInterval {times['write_interval']:.8g};
purgeWrite {int(cfg['run']['purge_write'])};
writeFormat ascii;
writePrecision 10;
writeCompression off;
timeFormat general;
timePrecision 8;
runTimeModifiable true;
adjustTimeStep off;

libs
(
    "libturbinesFoam.so"
);
"""


def render_fv_schemes(cfg: dict[str, Any], solver: str = "wale") -> str:
    """Central-differencing LES schemes; URANS switches to upwind + backward.

    The WALE headline follows the paper's second-order central differencing
    with CrankNicolson 1.0; the URANS fallback is the conventional bounded
    upwind / backward-Euler smoke configuration.
    """
    check_solver(solver)
    if solver == "urans":
        ddt = "backward"
        div_u = "bounded Gauss linearUpwind grad(U)"
        extra = (
            "    div(phi,k)      bounded Gauss upwind;\n"
            "    div(phi,omega)  bounded Gauss upwind;\n"
        )
    else:
        ddt = "CrankNicolson 1"
        div_u = "Gauss linear"
        extra = ""
    return foam_header("fvSchemes") + f"""ddtSchemes
{{
    default         {ddt};
}}

gradSchemes
{{
    default         Gauss linear;
    grad(U)         cellLimited Gauss linear 1;
}}

divSchemes
{{
    default         none;
    div(phi,U)      {div_u};
{extra}    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}}

laplacianSchemes
{{
    default         Gauss linear corrected;
}}

interpolationSchemes
{{
    default         linear;
}}

snGradSchemes
{{
    default         corrected;
}}
"""


def render_fv_solution(cfg: dict[str, Any]) -> str:
    """GAMG pressure with a reference cell (the domain is fully periodic)."""
    del cfg
    return foam_header("fvSolution") + """solvers
{
    p
    {
        solver          GAMG;
        tolerance       1e-7;
        relTol          0.05;
        smoother        DICGaussSeidel;
        cacheAgglomeration true;
        agglomerator    faceAreaPair;
        mergeLevels     1;
    }

    pFinal
    {
        $p;
        tolerance       1e-8;
        relTol          0;
    }

    "(U|k|omega)"
    {
        solver          smoothSolver;
        smoother        symGaussSeidel;
        tolerance       1e-7;
        relTol          0.05;
    }

    "(U|k|omega)Final"
    {
        $U;
        tolerance       1e-8;
        relTol          0;
    }
}

PIMPLE
{
    momentumPredictor          yes;
    nOuterCorrectors           1;
    nCorrectors                2;
    nNonOrthogonalCorrectors   0;
    pRefCell                   0;
    pRefValue                  0;
}

relaxationFactors
{
    equations
    {
        U       0.9;
        k       0.7;
        omega   0.7;
    }
}
"""


def render_decompose_par(cfg: dict[str, Any], ranks: int | None = None) -> str:
    nproc = int(ranks) if ranks is not None else int(
        cfg["decomposition"]["number_of_subdomains"]
    )
    return foam_header("decomposeParDict") + f"""numberOfSubdomains {nproc};
method scotch;
"""


def render_transport_properties(cfg: dict[str, Any]) -> str:
    nu = float(cfg["inflow"]["kinematic_viscosity"])
    return foam_header("transportProperties") + f"""transportModel Newtonian;

nu              nu [0 2 -1 0 0 0 0] {nu:.8g};
"""


def render_turbulence_properties(cfg: dict[str, Any], solver: str = "wale") -> str:
    """LES WALE headline; URANS k-omega SST fallback (weaker claim)."""
    check_solver(solver)
    if solver == "urans":
        model = cfg["solver"]["urans_model"]
        return foam_header("turbulenceProperties") + f"""simulationType RAS;

RAS
{{
    RASModel {model};
    turbulence on;
    printCoeffs on;
}}
"""
    model = cfg["solver"]["wale_model"]
    return foam_header("turbulenceProperties") + f"""simulationType LES;

LES
{{
    LESModel {model};
    turbulence on;
    printCoeffs on;
    delta cubeRootVol;
}}
"""


def render_field(cfg: dict[str, Any], field: str) -> str:
    """0.org field dictionaries: cyclic streamwise, symmetry crosswise."""
    velocity = float(cfg["inflow"]["velocity"])
    fields = inflow_fields(cfg)
    if field == "U":
        dims = "[0 1 -1 0 0 0 0]"
        header = "volVectorField"
        internal = f"uniform ({velocity:.8g} 0 0)"
    elif field == "p":
        dims = "[0 2 -2 0 0 0 0]"
        header = "volScalarField"
        internal = "uniform 0"
    elif field == "k":
        dims = "[0 2 -2 0 0 0 0]"
        header = "volScalarField"
        internal = f"uniform {fields['k']:.8g}"
    elif field == "omega":
        dims = "[0 0 -1 0 0 0 0]"
        header = "volScalarField"
        internal = f"uniform {fields['omega']:.8g}"
    elif field == "nut":
        dims = "[0 2 -1 0 0 0 0]"
        header = "volScalarField"
        internal = "uniform 0"
    else:
        raise ValueError(f"unknown field {field!r}")

    cyclic = "inlet\n    {\n        type cyclic;\n    }\n\n    outlet\n    {\n        type cyclic;\n    }"
    sides = "\n".join(
        f"    {name}\n    {{\n        type symmetry;\n    }}"
        for name in SIDE_PATCHES
    )
    return foam_header(field, header) + f"""dimensions      {dims};

internalField   {internal};

boundaryField
{{
    {cyclic}

{sides}
}}
"""


def render_fv_options(cfg: dict[str, Any], case_dir: Path) -> str:
    """Standalone nacelle actuator surface source (no turbine)."""
    nacelle = cfg["nacelle"]
    inflow = cfg["inflow"]
    radius = float(nacelle["radius"])
    reference_area = math.pi * radius**2
    geometry = geometry_path(cfg, case_dir)
    return foam_header("fvOptions") + f"""// Standalone nacelle actuator surface source (rotor-less validation case):
// the nacelle is represented only by this source, which distributes the
// surface force onto the background grid (no body-fitted mesh).
nacelle
{{
    type            nacelleSurfaceSource;
    active          true;

    nacelleSurfaceSourceCoeffs
    {{
        fieldNames          (U);
        selectionMode       cellSet;
        cellSet             nacelleCells;

        geometry            "{geometry}";
        referenceVelocity   {float(inflow['velocity']):.8g};
        rho                 1.0;
        nu                  {float(inflow['kinematic_viscosity']):.8g};
        cfModel             schultzGrunow;
        cf                  -1.0;
        bodyOrigin          (0 0 0);
        bodyAxis            (1 0 0);
        referenceArea       {reference_area:.12g};
        writeForceField     true;
        writePerf           true;
        writeNodePerf       false;
    }}
}}
"""


def outputs(
    cfg: dict[str, Any],
    mesh: str,
    solver: str,
    case_dir: Path,
    ranks: int | None = None,
) -> dict[Path, str]:
    check_solver(solver)
    system = case_dir / "system"
    constant = case_dir / "constant"
    zero = case_dir / "0.org"
    rendered: dict[Path, str] = {
        system / "blockMeshDict": render_block_mesh(cfg, mesh),
        system / "topoSetDict": render_toposet(cfg),
        system / "controlDict": render_control_dict(cfg),
        system / "decomposeParDict": render_decompose_par(cfg, ranks),
        system / "fvSchemes": render_fv_schemes(cfg, solver),
        system / "fvSolution": render_fv_solution(cfg),
        system / "fvOptions": render_fv_options(cfg, case_dir),
        constant / "transportProperties": render_transport_properties(cfg),
        constant / "turbulenceProperties": render_turbulence_properties(cfg, solver),
    }
    for field in FIELDS:
        rendered[zero / field] = render_field(cfg, field)
    return rendered


def check_outputs(rendered: dict[Path, str], case_dir: Path) -> bool:
    clean = True
    for path, expected in rendered.items():
        if not path.exists():
            print(f"missing generated file: {_display(path, case_dir)}", file=sys.stderr)
            clean = False
            continue
        current = path.read_text(encoding="utf-8")
        if current != expected:
            print(f"stale generated file: {_display(path, case_dir)}", file=sys.stderr)
            diff = difflib.unified_diff(
                current.splitlines(),
                expected.splitlines(),
                fromfile="current",
                tofile="expected",
                n=2,
            )
            print("\n".join(list(diff)[:40]), file=sys.stderr)
            clean = False
    return clean


def _display(path: Path, case_dir: Path) -> str:
    try:
        return str(path.relative_to(case_dir))
    except ValueError:
        return str(path)


def _positive_int(value: str) -> int:
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--mesh", choices=MESH_CHOICES, default="coarse")
    parser.add_argument(
        "--solver",
        choices=SOLVER_CHOICES,
        default="wale",
        help="wale renders the LES WALE headline; urans renders the "
             "k-omega SST smoke fallback",
    )
    parser.add_argument(
        "--ranks",
        type=_positive_int,
        default=None,
        metavar="N",
        help="override decomposition.number_of_subdomains in decomposeParDict",
    )
    parser.add_argument("--case-dir", type=Path, default=DEFAULT_CASE_DIR)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail (exit 1) if generated files are missing or stale; never writes",
    )
    args = parser.parse_args(argv)

    try:
        cfg = load_config(args.config)
        case_dir = args.case_dir.resolve()
        rendered = outputs(cfg, args.mesh, args.solver, case_dir, args.ranks)
    except (KeyError, ValueError) as exc:
        print(f"case generation error: {exc}", file=sys.stderr)
        return 2

    if args.check:
        return 0 if check_outputs(rendered, case_dir) else 1

    for path, content in rendered.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
        print(_display(path, ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
