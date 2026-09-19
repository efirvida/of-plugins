#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Adapted from mttbrbr/single-actuator-line (main 8284be8c...), GPL-3.0-or-later; see README.md "Attribution and license".
"""Render the NREL Phase VI uniform-inflow case from `config/case.yaml`.

Every generated file carries a "Generated from config/case.yaml" banner and
`--check` re-renders in memory and fails when a file on disk is missing or
stale, without writing anything.

Usage:
    generate_case.py [--mesh coarse|fine|ultra] [--speed 7|10|13|15|20|25]
                     [--domain long|squat] [--profile production|smoke]
                     [--solver urans|iddes] [--n-chordwise N] [--ranks N]
                     [--case-dir DIR] [--end-revs FLOAT]
                     [--start-from startTime|latestTime] [--check]

Adapted from `mttbrbr/single-actuator-line` (main 8284be8c..., GPL-3.0-or-later).
"""

from __future__ import annotations

import argparse
import difflib
import os
import sys
from pathlib import Path
from typing import Any

from case_config import (
    DEFAULT_CASE_DIR,
    DEFAULT_CONFIG,
    POLARS_DIR,
    ROOT,
    SOLVER_CHOICES,
    axis_extents_D,
    foam_header,
    foam_vector,
    inflow_fields,
    kinematics,
    load_config,
    turbine_origin,
)
from element_data import (
    blade_element_profiles,
    hub_element_rows,
    render_blade_element_data,
)

MESH_CHOICES = ("coarse", "fine", "ultra")
DOMAIN_CHOICES = ("long", "squat")
PROFILE_CHOICES = ("production", "smoke")
START_FROM_CHOICES = ("startTime", "latestTime")

ALM_ELEMENT = "actuatorLineElement"
ASM_ELEMENT = "actuatorSurfaceElement"

PATCH_NAMES = ("inlet", "outlet", "bottom", "top", "sideMinus", "sidePlus")


def render_block_mesh(cfg: dict[str, Any], mesh: str, domain: str, profile: str) -> str:
    """18-block hexahedral blockMeshDict with six far-field patches."""
    resolution = cfg["mesh"]["resolutions"][mesh]
    breaks = {axis: axis_extents_D(cfg, axis, domain) for axis in ("x", "y", "z")}
    cells = {
        axis: [int(value) for value in resolution[axis]["cells"]]
        for axis in ("x", "y", "z")
    }
    grading = {
        axis: [float(value) for value in resolution[axis]["grading"]]
        for axis in ("x", "y", "z")
    }
    if profile == "smoke":
        cells = {
            axis: [max(1, round(value / 4)) for value in cells[axis]]
            for axis in ("x", "y", "z")
        }
    nx, ny, nz = (len(breaks[axis]) for axis in ("x", "y", "z"))

    def vertex(ix: int, iy: int, iz: int) -> int:
        return (iz * ny + iy) * nx + ix

    vertices = []
    for z in breaks["z"]:
        for y in breaks["y"]:
            for x in breaks["x"]:
                vertices.append(f"    ({x:.9g} {y:.9g} {z:.9g})")

    blocks = []
    faces: dict[str, list[str]] = {name: [] for name in PATCH_NAMES}
    for iz in range(nz - 1):
        for iy in range(ny - 1):
            for ix in range(nx - 1):
                v000 = vertex(ix, iy, iz)
                v100 = vertex(ix + 1, iy, iz)
                v110 = vertex(ix + 1, iy + 1, iz)
                v010 = vertex(ix, iy + 1, iz)
                v001 = vertex(ix, iy, iz + 1)
                v101 = vertex(ix + 1, iy, iz + 1)
                v111 = vertex(ix + 1, iy + 1, iz + 1)
                v011 = vertex(ix, iy + 1, iz + 1)
                blocks.append(
                    "    hex "
                    f"({v000} {v100} {v110} {v010} {v001} {v101} {v111} {v011}) "
                    f"({cells['x'][ix]} {cells['y'][iy]} {cells['z'][iz]}) "
                    f"simpleGrading ({grading['x'][ix]:.9g} "
                    f"{grading['y'][iy]:.9g} {grading['z'][iz]:.9g})"
                )
                if ix == 0:
                    faces["inlet"].append(f"            ({v000} {v001} {v011} {v010})")
                if ix == nx - 2:
                    faces["outlet"].append(f"            ({v100} {v110} {v111} {v101})")
                if iz == 0:
                    faces["bottom"].append(f"            ({v000} {v010} {v110} {v100})")
                if iz == nz - 2:
                    faces["top"].append(f"            ({v001} {v101} {v111} {v011})")
                if iy == 0:
                    faces["sideMinus"].append(f"            ({v000} {v100} {v101} {v001})")
                if iy == ny - 2:
                    faces["sidePlus"].append(f"            ({v010} {v011} {v111} {v110})")

    patch_types = {
        "inlet": "patch",
        "outlet": "patch",
        "bottom": "symmetryPlane",
        "top": "symmetryPlane",
        "sideMinus": "symmetryPlane",
        "sidePlus": "symmetryPlane",
    }

    def patch(name: str) -> str:
        return (
            f"    {name}\n"
            "    {\n"
            f"        type {patch_types[name]};\n"
            "        faces\n"
            "        (\n"
            + "\n".join(faces[name])
            + "\n        );\n"
            "    }"
        )

    return (
        foam_header("blockMeshDict")
        + "scale 1;\n\nvertices\n(\n"
        + "\n".join(vertices)
        + "\n);\n\nblocks\n(\n"
        + "\n".join(blocks)
        + "\n);\n\nedges ();\n\nboundary\n(\n"
        + "\n".join(patch(name) for name in PATCH_NAMES)
        + "\n);\n\nmergePatchPairs ();\n"
    )


def render_toposet(cfg: dict[str, Any]) -> str:
    """One cellSet box enclosing the rotor disk and hub."""
    diameter = float(cfg["turbine"]["diameter"])
    origin = turbine_origin(cfg)
    minimum = (
        origin[0] - 0.5 * diameter,
        origin[1] - 1.0 * diameter,
        origin[2] - 1.0 * diameter,
    )
    maximum = (
        origin[0] + 0.5 * diameter,
        origin[1] + 1.0 * diameter,
        origin[2] + 1.0 * diameter,
    )
    return foam_header("topoSetDict") + f"""actions
(
    {{
        name turbine;
        type cellSet;
        action new;
        source boxToCell;
        box {foam_vector(minimum)} {foam_vector(maximum)};
    }}
);
"""


def render_control_dict(
    cfg: dict[str, Any],
    speed: str,
    mesh: str,
    profile: str,
    end_revs: float | None,
    start_from: str,
    sequence: str = "H",
    solver: str = "urans",
) -> str:
    values = kinematics(cfg, speed, mesh, sequence, solver)
    solver_cfg = cfg["solver"]
    if profile == "smoke":
        end_time = 2.0 * values["t_rev"]
        purge_write = 4
    else:
        end_time = (end_revs if end_revs is not None else float(solver_cfg["end_revolutions"])) * values["t_rev"]
        purge_write = int(solver_cfg["purge_write"])
    return foam_header("controlDict") + f"""application {solver_cfg['application']};
startFrom {start_from};
startTime 0;
stopAt endTime;
endTime {end_time:.8g};
deltaT {values['delta_t']:.8g};
writeControl runTime;
writeInterval {values['write_interval']:.8g};
purgeWrite {purge_write};
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


def render_fv_schemes(cfg: dict[str, Any], solver: str = "urans") -> str:
    """Cartesian-mesh schemes; backward ddt, linearUpwind for U.

    The IDDES variant switches U convection to the configured (non-dissipative)
    scheme and asks `wallDist` for the wall-normal vectors up front, because
    `IDDESDelta` reads `wallDist::n()`.
    """
    div_phi_u = "bounded Gauss linearUpwind grad(U)"
    wall_dist = "    method          meshWave;"
    if solver == "iddes":
        iddes = cfg["iddes"]
        div_phi_u = str(iddes["div_phi_U"])
        if bool(iddes["wall_dist_n_required"]):
            wall_dist = "    method          meshWave;\n    nRequired       true;"
    return foam_header("fvSchemes") + f"""ddtSchemes
{{
    default         backward;
}}

gradSchemes
{{
    default         Gauss linear;
    grad(U)         cellLimited Gauss linear 1;
}}

divSchemes
{{
    default         none;
    div(phi,U)      {div_phi_u};
    div(phi,k)      bounded Gauss upwind;
    div(phi,omega)  bounded Gauss upwind;
    div((nuEff*dev2(T(grad(U))))) Gauss linear;
}}

laplacianSchemes
{{
    default         Gauss linear orthogonal;
}}

interpolationSchemes
{{
    default         linear;
}}

snGradSchemes
{{
    default         orthogonal;
}}

wallDist
{{
{wall_dist}
}}
"""


def render_fv_solution(cfg: dict[str, Any]) -> str:
    """GAMG pressure, smoothSolver for the transported fields."""
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


def render_turbulence_properties(cfg: dict[str, Any], solver: str = "urans") -> str:
    """RAS kOmegaSST by default; LES kOmegaSSTIDDES for the Stage 3 variant.

    `kOmegaSSTIDDES::setDelta()` accepts only an IDDESDelta-based delta model
    (validated in `case_config.validate_config`), configured here as
    `delta IDDESDelta` with its `IDDESDeltaCoeffs` sub-dictionary.
    """
    if solver == "iddes":
        iddes = cfg["iddes"]
        model = cfg["solver"]["iddes_model"]
        return foam_header("turbulenceProperties") + f"""simulationType LES;

LES
{{
    LESModel {model};
    turbulence on;
    printCoeffs on;
    delta {iddes['delta']};

    {iddes['delta']}Coeffs
    {{
        Cw {float(iddes['delta_Cw']):.8g};
    }}
}}
"""
    model = cfg["solver"]["turbulence_model"]
    return foam_header("turbulenceProperties") + f"""simulationType RAS;

RAS
{{
    RASModel {model};
    turbulence on;
    printCoeffs on;
}}
"""


def render_field_uniform(cfg: dict[str, Any], speed: str, field: str) -> str:
    """0.org field dictionaries (uniform inflow, symmetry far field)."""
    velocity = float(speed)
    fields = inflow_fields(cfg, speed)
    far = "        type symmetryPlane;"
    if field == "U":
        dims = "[0 1 -1 0 0 0 0]"
        header = "volVectorField"
        internal = f"uniform ({velocity:.8g} 0 0)"
        inlet = f"type fixedValue;\n        value uniform ({velocity:.8g} 0 0);"
        outlet = (
            "type inletOutlet;\n"
            "        inletValue uniform (0 0 0);\n"
            f"        value uniform ({velocity:.8g} 0 0);"
        )
    elif field == "p":
        dims = "[0 2 -2 0 0 0 0]"
        header = "volScalarField"
        internal = "uniform 0"
        inlet = "type zeroGradient;"
        outlet = "type fixedValue;\n        value uniform 0;"
    elif field == "k":
        dims = "[0 2 -2 0 0 0 0]"
        header = "volScalarField"
        value = fields["k"]
        internal = f"uniform {value:.8g}"
        inlet = f"type fixedValue;\n        value uniform {value:.8g};"
        outlet = (
            "type inletOutlet;\n"
            "        inletValue uniform 0;\n"
            f"        value uniform {value:.8g};"
        )
    elif field == "omega":
        dims = "[0 0 -1 0 0 0 0]"
        header = "volScalarField"
        value = fields["omega"]
        internal = f"uniform {value:.8g}"
        inlet = f"type fixedValue;\n        value uniform {value:.8g};"
        outlet = (
            "type inletOutlet;\n"
            "        inletValue uniform 0;\n"
            f"        value uniform {value:.8g};"
        )
    elif field == "nut":
        dims = "[0 2 -1 0 0 0 0]"
        header = "volScalarField"
        internal = "uniform 0"
        inlet = "type fixedValue;\n        value uniform 0;"
        outlet = (
            "type inletOutlet;\n"
            "        inletValue uniform 0;\n"
            "        value uniform 0;"
        )
    else:
        raise ValueError(f"unknown field {field!r}")

    far_block = "\n".join(
        f"    {name}\n    {{\n{far}\n    }}"
        for name in ("bottom", "top", "sideMinus", "sidePlus")
    )
    return foam_header(field, header) + f"""dimensions      {dims};

internalField   {internal};

boundaryField
{{
    inlet
    {{
        {inlet}
    }}

    outlet
    {{
        {outlet}
    }}

{far_block}
}}
"""


def render_fv_options(
    cfg: dict[str, Any],
    speed: str,
    mesh: str,
    case_dir: Path,
    element_type: str,
    sequence: str = "H",
    n_chordwise: int | None = None,
) -> str:
    """ALM/ASM twins rendered from one function.

    The twins differ only in the blade element keys: `elementType` and, for
    the surface element, `nChordwise`. Path targets are found by
    `test_twins_differ_only_in_blade_keys`. `n_chordwise` overrides the
    configured ASM strip count (the ALM has no such key).
    """
    values = kinematics(cfg, speed, mesh, sequence)
    turbine = cfg["turbine"]
    actuator = cfg["actuator"]
    origin = turbine_origin(cfg)
    polars = os.path.relpath(POLARS_DIR / "S809_OSU_Re1M_total.dat", case_dir / "system")
    profiles = " ".join(blade_element_profiles(cfg))
    if n_chordwise is None:
        n_chordwise = int(actuator["n_chordwise"])
    elif int(n_chordwise) <= 0:
        raise ValueError("nChordwise must be a positive integer")
    else:
        n_chordwise = int(n_chordwise)
    element_keys = f"                elementType {element_type};\n"
    if element_type == ASM_ELEMENT:
        element_keys += f"                nChordwise {n_chordwise};\n"
    dynamic = actuator["dynamic_stall"]
    end_effects = actuator["end_effects"]
    hub_rows = "\n".join(
        "                (" + " ".join(f"{value:.9g}" for value in row) + ")"
        for row in hub_element_rows()
    )
    return foam_header("fvOptions") + f"""turbine
{{
    type axialFlowTurbineALSource;
    active on;

    axialFlowTurbineALSourceCoeffs
    {{
        fieldNames (U);
        selectionMode cellSet;
        cellSet turbine;
        origin {foam_vector(origin)};
        axis {foam_vector(turbine['axis'])};
        verticalDirection {foam_vector(turbine['vertical_direction'])};
        freeStreamVelocity ({values['speed']:.8g} 0 0);
        tipSpeedRatio {values['tsr']:.8g};
        rotorRadius {float(turbine['radius']):.8g};
        azimuthalOffset 0;

        dynamicStall
        {{
            active {'on' if dynamic['active'] else 'off'};
            dynamicStallModel {dynamic['model']};
        }}

        endEffects
        {{
            active {'on' if end_effects['active'] else 'off'};
            endEffectsModel {end_effects['model']};
            GlauertCoeffs
            {{
                tipEffects {'on' if end_effects['tip'] else 'off'};
                rootEffects {'on' if end_effects['root'] else 'off'};
            }}
        }}

        blades
        {{
            blade1
            {{
                writePerf true;
                writeElementPerf true;
{element_keys.rstrip()}
                nElements {int(turbine['n_elements'])};
                elementProfiles
                (
                    {profiles}
                );
                elementData
                (
{render_blade_element_data(cfg)}
                );
            }}
            blade2
            {{
                $blade1;
                writePerf false;
                writeElementPerf false;
                azimuthalOffset 180;
            }}
        }}

        hub
        {{
            nElements 4;
            elementProfiles (cylinder);
            elementData
            (
{hub_rows}
            );
        }}

        profileData
        {{
            S809
            {{
                Re 1e6;
                GaussianCoeffs
                {{
                    chordFactor 0.25;
                    dragFactor 1.0;
                    meshFactor {float(turbine['gaussian_mesh_factor']):.8g};
                }}
                data (#include "{polars}");
            }}
            cylinder {{ data ((-180 0 1.1 0) (180 0 1.1 0)); }}
        }}
    }}
}}
"""


def outputs(
    cfg: dict[str, Any],
    mesh: str,
    speed: str,
    domain: str,
    profile: str,
    case_dir: Path,
    end_revs: float | None,
    start_from: str,
    sequence: str = "H",
    solver: str = "urans",
    n_chordwise: int | None = None,
    ranks: int | None = None,
) -> dict[Path, str]:
    system = case_dir / "system"
    constant = case_dir / "constant"
    zero = case_dir / "0.org"
    rendered: dict[Path, str] = {
        system / "blockMeshDict": render_block_mesh(cfg, mesh, domain, profile),
        system / "topoSetDict": render_toposet(cfg),
        system / "controlDict": render_control_dict(
            cfg, speed, mesh, profile, end_revs, start_from, sequence, solver
        ),
        system / "decomposeParDict": render_decompose_par(cfg, ranks),
        system / "fvSchemes": render_fv_schemes(cfg, solver),
        system / "fvSolution": render_fv_solution(cfg),
        system / "fvOptions.ALM": render_fv_options(
            cfg, speed, mesh, case_dir, ALM_ELEMENT, sequence, n_chordwise
        ),
        system / "fvOptions.ASM": render_fv_options(
            cfg, speed, mesh, case_dir, ASM_ELEMENT, sequence, n_chordwise
        ),
        constant / "transportProperties": render_transport_properties(cfg),
        constant / "turbulenceProperties": render_turbulence_properties(cfg, solver),
    }
    for field in ("U", "p", "k", "omega", "nut"):
        rendered[zero / field] = render_field_uniform(cfg, speed, field)
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
    parser.add_argument("--speed", default="7")
    parser.add_argument("--sequence", choices=("H", "S"), default="H",
                        help="Sequence H headline or the Sequence S 7 m/s repeat")
    parser.add_argument("--domain", choices=DOMAIN_CHOICES, default="long")
    parser.add_argument("--profile", choices=PROFILE_CHOICES, default="production")
    parser.add_argument(
        "--solver",
        choices=SOLVER_CHOICES,
        default="urans",
        help="urans renders RAS kOmegaSST; iddes renders the Stage 3 LES "
             "kOmegaSSTIDDES variant with its own fixed time step",
    )
    parser.add_argument(
        "--n-chordwise",
        type=_positive_int,
        default=None,
        metavar="N",
        help="override the configured ASM chordwise strip count (ASM twin only)",
    )
    parser.add_argument(
        "--ranks",
        type=_positive_int,
        default=None,
        metavar="N",
        help="override decomposition.number_of_subdomains in decomposeParDict",
    )
    parser.add_argument("--case-dir", type=Path, default=DEFAULT_CASE_DIR)
    parser.add_argument("--end-revs", type=float, default=None)
    parser.add_argument("--start-from", choices=START_FROM_CHOICES, default="startTime")
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail (exit 1) if generated files are missing or stale; never writes",
    )
    args = parser.parse_args(argv)

    try:
        cfg = load_config(args.config)
        case_dir = args.case_dir.resolve()
        rendered = outputs(
            cfg,
            args.mesh,
            args.speed,
            args.domain,
            args.profile,
            case_dir,
            args.end_revs,
            args.start_from,
            args.sequence,
            solver=args.solver,
            n_chordwise=args.n_chordwise,
            ranks=args.ranks,
        )
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
