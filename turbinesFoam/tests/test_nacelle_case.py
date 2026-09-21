# SPDX-License-Identifier: GPL-3.0-or-later
"""Pure-Python contracts for the periodic-nacelle actuator-surface validation case.

These tests never require OpenFOAM: they exercise the YAML single source of
truth, the derived mesh/run arithmetic, the rendered dictionaries in memory and
the non-destructive `--check` gate (spec `nacelle-validation-case`).

The tool module is loaded under a private name because `validation/phaseVI/`
ships a module with the same file name (`tools/generate_case.py`).
"""

from __future__ import annotations

import importlib.util
import re
import subprocess
import sys
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).resolve().parent
TURBINESFOAM = TESTS_DIR.parent
PACKAGE = TURBINESFOAM / "validation" / "nacelle-asn"
CONFIG = PACKAGE / "config" / "case.yaml"
CASE_DIR = PACKAGE / "case"
GENERATOR = PACKAGE / "tools" / "generate_case.py"
NACELLE_STL = TURBINESFOAM / "geometry" / "stl" / "nacelle.stl"

MODULE_NAME = "nacelle_generate_case"
_spec = importlib.util.spec_from_file_location(MODULE_NAME, GENERATOR)
generate_case = importlib.util.module_from_spec(_spec)
sys.modules[MODULE_NAME] = generate_case
_spec.loader.exec_module(generate_case)

EXPECTED_CELLS = {"coarse": (153, 80, 80), "medium": (115, 151, 151)}
REFERENCE_CELLS = (502, 348, 348)
FIELDS = ("U", "p", "k", "omega", "nut")
PATCH_NAMES = ("inlet", "outlet", "bottom", "top", "sideMinus", "sidePlus")


@pytest.fixture(scope="module")
def cfg():
    return generate_case.load_config(CONFIG)


def _render_into(case_dir: Path, *extra: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(GENERATOR), "--case-dir", str(case_dir), *extra],
        capture_output=True,
        text=True,
    )


def _check(case_dir: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(GENERATOR), "--check", "--case-dir", str(case_dir)],
        capture_output=True,
        text=True,
    )


def test_config_schema(cfg):
    assert cfg["project"]["openfoam_versions"] == ["v2506", "v2412"]
    nacelle = cfg["nacelle"]
    assert nacelle["radius"] == 1.0
    assert nacelle["cylinder_ratio"] == 6.0
    assert (PACKAGE / nacelle["stl"]).resolve() == NACELLE_STL
    assert NACELLE_STL.is_file()

    domain = cfg["domain_R"]
    assert domain["x"] == [-10.0, 20.0]
    assert domain["y"] == [-10.0, 10.0]
    assert domain["z"] == [-10.0, 10.0]
    assert cfg["boundaries"] == {"streamwise": "cyclic", "crosswise": "symmetry"}

    inflow = cfg["inflow"]
    re = inflow["velocity"] * nacelle["radius"] / inflow["kinematic_viscosity"]
    assert re == pytest.approx(1000.0)

    solver = cfg["solver"]
    assert solver["application"] == "pimpleFoam"
    assert solver["headline"] == "wale"
    assert solver["fallback"] == "urans"
    assert solver["adjust_time_step"] is False

    # The nacelle downstream end sits at the streamwise origin; the whole body
    # stays inside the periodic cell.
    nose, downstream_end = generate_case.nacelle_bounds_R(cfg)
    assert nose == pytest.approx(-7.0)
    assert downstream_end == pytest.approx(0.0)
    assert domain["x"][0] < nose < downstream_end < domain["x"][1]


def test_grids_and_reference_declared_only(cfg):
    for mesh, expected in EXPECTED_CELLS.items():
        assert generate_case.cell_counts(cfg, mesh) == expected
        assert generate_case.cell_count(cfg, mesh) == (
            expected[0] * expected[1] * expected[2]
        )
    assert tuple(cfg["mesh"]["reference"]["cells"]) == REFERENCE_CELLS
    assert cfg["mesh"]["reference"]["declared_only"] is True
    assert "reference" not in cfg["mesh"]["resolutions"]
    assert "reference" not in generate_case.MESH_CHOICES

    # Uniform spacings follow from the domain and the counts.
    assert generate_case.cell_size_R(cfg, "coarse", "x") == pytest.approx(30.0 / 153)
    assert generate_case.cell_size_R(cfg, "coarse", "y") == pytest.approx(20.0 / 80)
    assert generate_case.cell_size_R(cfg, "medium", "y") == pytest.approx(20.0 / 151)


def test_block_mesh_dict(cfg):
    text = generate_case.render_block_mesh(cfg, "coarse")
    assert "hex (0 1 2 3 4 5 6 7) (153 80 80) simpleGrading (1 1 1)" in text
    assert text.count("hex (") == 1
    assert "type cyclic;" in text
    assert "neighbourPatch outlet;" in text
    assert "neighbourPatch inlet;" in text
    assert text.count("type symmetry;") == 4
    for name in PATCH_NAMES:
        assert f"    {name}\n" in text
    for vertex in ("(-10 -10 -10)", "(20 10 10)"):
        assert vertex in text

    medium = generate_case.render_block_mesh(cfg, "medium")
    assert "hex (0 1 2 3 4 5 6 7) (115 151 151)" in medium


def test_control_dict_run_schedule(cfg):
    text = generate_case.render_control_dict(cfg)
    assert "application pimpleFoam;" in text
    assert "deltaT 0.1;" in text
    assert "adjustTimeStep off;" in text
    assert "endTime 210;" in text
    assert "writeInterval 1;" in text
    assert "purgeWrite 0;" in text
    assert '"libturbinesFoam.so"' in text

    times = generate_case.run_times(cfg)
    assert times["t_ft"] == pytest.approx(30.0)
    assert times["wash_out_end"] == pytest.approx(60.0)
    assert times["average_window"] == pytest.approx(150.0)
    assert times["end_time"] == pytest.approx(210.0)
    assert times["delta_t"] == pytest.approx(0.1)
    assert times["samples"] == 150


def test_sampling_function_object(cfg):
    """The rendered case probes the metric-station lines every time step.

    `scripts/compareNacelle.py` reads those series for <u>(z) and the resolved
    TKE k(z), so the sampling contract is part of the case, not of the compare
    tool.
    """
    assert cfg["inflow"]["density"] == 1.0
    sampling = cfg["sampling"]
    low, high = sampling["z_R"]
    assert low < 0.0 < high
    assert sampling["dz_R"] > 0.0

    stations = generate_case.sampling_stations_R(cfg)
    assert stations == [1.0, 3.0, 5.0, 7.0, 10.0]

    for mesh in ("coarse", "medium"):
        text = generate_case.render_control_dict(cfg, mesh)
        assert "type            probes;" in text
        assert "libs            (sampling);" in text
        # OpenFOAM's `probes` writes to postProcessing/<function-object name>/,
        # so the dictionary key is the contract with compareNacelle.py; a `name`
        # entry is ignored by the function object and must not come back.
        assert re.search(r"^\s{4}profiles\s*$", text, re.MULTILINE)
        assert not re.search(r"^\s{4}name\s", text, re.MULTILINE)
        assert "fields          (U);" in text
        assert "writeControl    timeStep;" in text
        assert "probeLocations" in text

        locations = generate_case.probe_locations(cfg, mesh)
        assert len(set(locations)) == len(locations)  # no duplicate probe
        station_lines = {round(point[0], 9) for point in locations}
        assert len(station_lines) == len(stations)
        z_span = max(point[2] for point in locations) - min(
            point[2] for point in locations
        )
        # The lines must cover the reference profiles (|z/D| <= 1.32, D = 2R).
        assert z_span > 2.0 * 1.32 * cfg["nacelle"]["radius"]
        for point in locations:
            for axis, value in zip(("x", "y", "z"), point):
                low, high = cfg["domain_R"][axis]
                assert low < value < high

    committed = (CASE_DIR / "system" / "controlDict").read_text(encoding="utf-8")
    assert committed == generate_case.render_control_dict(cfg, "coarse")
    assert "type            nacelleSurfaceSource;" in (
        CASE_DIR / "system" / "fvOptions"
    ).read_text(encoding="utf-8")


def test_turbulence_headline_and_fallback(cfg):
    les = generate_case.render_turbulence_properties(cfg)
    assert "simulationType LES;" in les
    assert "LESModel WALE;" in les
    assert "delta cubeRootVol;" in les
    assert "RASModel" not in les

    ras = generate_case.render_turbulence_properties(cfg, "urans")
    assert "simulationType RAS;" in ras
    assert "RASModel kOmegaSST;" in ras
    assert "LESModel" not in ras

    schemes = generate_case.render_fv_schemes(cfg)
    assert "default         CrankNicolson 1;" in schemes
    assert "div(phi,U)      Gauss linear;" in schemes
    urans_schemes = generate_case.render_fv_schemes(cfg, "urans")
    assert "bounded Gauss linearUpwind grad(U)" in urans_schemes
    assert "bounded Gauss upwind" in urans_schemes

    # The adaptation and the weaker fallback claim are documented.
    readme = (PACKAGE / "README.md").read_text(encoding="utf-8")
    assert "WALE" in readme
    assert "dynamic SGS" in readme
    assert "fallback" in readme
    assert "arXiv:1702.02108v4" in readme


def test_fv_options_standalone_source(cfg):
    text = generate_case.render_fv_options(cfg, CASE_DIR)
    assert "type            nacelleSurfaceSource;" in text
    assert "axialFlowTurbineALSource" not in text
    assert "fieldNames          (U);" in text
    assert "cellSet             nacelleCells;" in text
    assert "cfModel             schultzGrunow;" in text
    assert "referenceArea       3.14159265359;" in text
    assert "writePerf           true;" in text

    # The rendered geometry path resolves to the pipeline STL from the case
    # root (where the solver runs). Run directories are rendered with their own
    # `--case-dir`, so the path is recomputed for their depth.
    match = re.search(r'geometry\s+"([^"]+)";', text)
    assert match is not None
    assert (CASE_DIR / match.group(1)).resolve() == NACELLE_STL

    run_dir = PACKAGE / "runs" / "coarse-s0"
    run_text = generate_case.render_fv_options(cfg, run_dir)
    run_match = re.search(r'geometry\s+"([^"]+)";', run_text)
    assert run_match is not None
    assert (run_dir / run_match.group(1)).resolve() == NACELLE_STL

    committed = (CASE_DIR / "system" / "fvOptions").read_text(encoding="utf-8")
    assert committed == text


def test_toposet_covers_the_nacelle(cfg):
    text = generate_case.render_toposet(cfg)
    assert "name nacelleCells;" in text
    assert "type cellSet;" in text
    assert "box (-8 -2 -2) (1 2 2);" in text


def test_fields_boundary_conditions(cfg):
    for field in FIELDS:
        text = generate_case.render_field(cfg, field)
        assert text.count("type cyclic;") == 2
        assert text.count("type symmetry;") == 4
        for name in PATCH_NAMES:
            assert f"    {name}\n" in text

    u = generate_case.render_field(cfg, "U")
    assert "dimensions      [0 1 -1 0 0 0 0];" in u
    assert "internalField   uniform (1 0 0);" in u


def test_metrics_and_acceptance(cfg):
    metrics = cfg["metrics"]
    assert metrics["stations_R"] == [1.0, 3.0, 5.0, 7.0]
    assert metrics["stretch_R"] == 10.0
    assert metrics["profiles"] == ["u", "k"]
    assert metrics["permeable_disk_datum"] == pytest.approx(0.48)
    assert "F_drag" in metrics["drag_coefficient"]

    acceptance = cfg["acceptance"]
    assert acceptance["medium"] == "all_stations"
    assert acceptance["coarse"] == [1.0, 10.0]
    for station in (3.0, 5.0, 7.0):
        assert station not in acceptance["coarse"]


def test_stage_plan(cfg):
    stages = cfg["stages"]
    executing = [name for name, stage in stages.items() if stage.get("executes")]
    assert executing == ["stage0"]
    assert stages["stage0"]["queue"] == "sequana_cpu_dev"
    assert stages["production"]["queue"] == "sequana_cpu"
    assert stages["production"]["executes"] is False


def test_cli_variant_render(tmp_path):
    case_dir = tmp_path / "case"
    result = _render_into(case_dir, "--mesh", "medium", "--solver", "urans", "--ranks", "96")
    assert result.returncode == 0, result.stderr

    block_mesh = (case_dir / "system" / "blockMeshDict").read_text(encoding="utf-8")
    assert "hex (0 1 2 3 4 5 6 7) (115 151 151)" in block_mesh
    turbulence = (case_dir / "constant" / "turbulenceProperties").read_text(
        encoding="utf-8"
    )
    assert "RASModel kOmegaSST;" in turbulence
    decompose = (case_dir / "system" / "decomposeParDict").read_text(encoding="utf-8")
    assert "numberOfSubdomains 96;" in decompose


def test_generated_case_is_current():
    result = _check(CASE_DIR)
    assert result.returncode == 0, result.stderr


def test_stale_file_detection_is_non_destructive(tmp_path):
    case_dir = tmp_path / "case"
    assert _render_into(case_dir).returncode == 0
    assert _check(case_dir).returncode == 0

    control = case_dir / "system" / "controlDict"
    tampered = control.read_text(encoding="utf-8") + "\n// hand edit\n"
    control.write_text(tampered, encoding="utf-8")

    result = _check(case_dir)
    assert result.returncode != 0
    assert "stale generated file" in result.stderr
    assert "system/controlDict" in result.stderr
    # only the tampered file is reported, and it is left untouched
    assert "stale generated file: 0.org/U" not in result.stderr
    assert control.read_text(encoding="utf-8") == tampered


def test_missing_file_detection_is_non_destructive(tmp_path):
    case_dir = tmp_path / "case"
    assert _render_into(case_dir).returncode == 0

    (case_dir / "0.org" / "U").unlink()
    result = _check(case_dir)
    assert result.returncode != 0
    assert "missing generated file" in result.stderr
    assert not (case_dir / "0.org" / "U").exists()
