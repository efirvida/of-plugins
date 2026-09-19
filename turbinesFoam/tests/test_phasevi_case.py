"""Pure-Python contracts for the NREL Phase VI validation case package.

These tests never require OpenFOAM: they exercise the YAML single source of
truth, the derived mesh arithmetic and the rendered dictionaries in memory
(spec requirements V1-V11, design sections 4-6).
"""

from __future__ import annotations

import copy
import re
import subprocess
import sys
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).resolve().parent
TURBINESFOAM = TESTS_DIR.parent
PACKAGE = TURBINESFOAM / "validation" / "phaseVI"
sys.path.insert(0, str(PACKAGE / "tools"))

import generate_case  # noqa: E402
from case_config import (  # noqa: E402
    assert_tip_constraint,
    assert_tsr_matches_experiment,
    cell_count,
    hub_cell_size,
    kinematics,
    load_config,
    select_speed,
    speed_keys,
    tip_displacement,
    validate_config,
)
from element_data import (  # noqa: E402
    HUB_DIAMETER,
    blade_element_rows,
    hub_element_rows,
)

CONFIG = PACKAGE / "config" / "case.yaml"
EXPECTED_CELLS = {
    "coarse": 6_674_304,
    "fine": 22_525_776,
    "ultra": 53_394_432,
}
EXPECTED_HUB_CELL = {"coarse": 0.3143, "fine": 0.2095, "ultra": 0.1572}
FROZEN_TSR = {
    "7": 5.408,
    "10": 3.798,
    "13": 2.920,
    "15": 2.530,
    "20": 1.898,
    "25": 1.521,
}
BLADE_KEYS = ("elementType", "nChordwise")


@pytest.fixture(scope="module")
def cfg():
    return load_config(CONFIG)


def test_config_schema(cfg):
    assert cfg["project"]["openfoam_versions"] == ["v2506", "v2412"]
    turbine = cfg["turbine"]
    assert turbine["n_blades"] == 2
    assert turbine["pitch_deg"] == 3.0
    assert turbine["diameter"] == 2.0 * turbine["radius"]
    assert turbine["hub_height"] == 12.192
    assert speed_keys(cfg) == ["7", "10", "13", "15", "20", "25"]
    assert cfg["solver"]["application"] == "pimpleFoam"
    assert cfg["solver"]["turbulence_model"] == "kOmegaSST"
    assert cfg["solver"]["adjust_time_step"] is False
    solver = cfg["solver"]
    assert (solver["discard_revolutions"], solver["end_revolutions"]) == (4, 12)


def test_mesh_arithmetic_and_cell_counts(cfg):
    for mesh, expected in EXPECTED_CELLS.items():
        assert cell_count(cfg, mesh) == expected
        low, high = cfg["mesh"]["resolutions"][mesh]["target_cells"]
        assert low <= expected <= high
        assert hub_cell_size(cfg, mesh) == pytest.approx(
            EXPECTED_HUB_CELL[mesh], abs=5e-4
        )


def test_block_mesh_dict_18_hex_blocks(cfg):
    pattern = re.compile(
        r"^\s+hex \(([^)]*)\) \((\d+) (\d+) (\d+)\) simpleGrading \(",
        re.MULTILINE,
    )
    for mesh in EXPECTED_CELLS:
        text = generate_case.render_block_mesh(cfg, mesh, "long", "production")
        blocks = pattern.findall(text)
        assert len(blocks) == 18, f"{mesh}: expected 18 blocks"
        assert "poly" not in text
        total = sum(int(nx) * int(ny) * int(nz) for _, nx, ny, nz in blocks)
        assert total == EXPECTED_CELLS[mesh]
        for name in ("inlet", "outlet", "bottom", "top", "sideMinus", "sidePlus"):
            assert f"    {name}\n" in text
        assert text.count("type symmetryPlane;") == 4
        assert text.count("type patch;") == 2
    squat = generate_case.render_block_mesh(cfg, "coarse", "squat", "production")
    assert len(pattern.findall(squat)) == 18


def test_fallback_domain_extents(cfg):
    from case_config import axis_extents_D

    diameter = cfg["turbine"]["diameter"]
    hub = cfg["turbine"]["hub_height"]
    long_x = axis_extents_D(cfg, "x", "long")
    long_y = axis_extents_D(cfg, "y", "long")
    long_z = axis_extents_D(cfg, "z", "long")
    assert (long_x[0], long_x[-1]) == (-5.0 * diameter, 15.0 * diameter)
    assert (long_y[0], long_y[-1]) == (-4.0 * diameter, 4.0 * diameter)
    assert (long_z[0], long_z[-1]) == (hub - 2.5 * diameter, hub + 2.5 * diameter)

    squat_x = axis_extents_D(cfg, "x", "squat")
    squat_y = axis_extents_D(cfg, "y", "squat")
    squat_z = axis_extents_D(cfg, "z", "squat")
    assert (squat_x[0], squat_x[-1]) == (-5.0 * diameter, 10.0 * diameter)
    assert (squat_y[0], squat_y[-1]) == (-3.0 * diameter, 3.0 * diameter)
    assert (squat_z[0], squat_z[-1]) == (hub - 3.0 * diameter, hub + 3.0 * diameter)
    assert cell_count(cfg, "coarse", "squat") == EXPECTED_CELLS["coarse"]


def test_control_dict_fixed_time_step(cfg):
    for mesh, delta_t in (("coarse", 0.008), ("fine", 0.005), ("ultra", 0.004)):
        text = generate_case.render_control_dict(
            cfg, "7", mesh, "production", None, "startTime"
        )
        assert f"deltaT {delta_t:.8g};" in text
        assert "adjustTimeStep off;" in text
        values = kinematics(cfg, "7", mesh)
        assert f"endTime {values['end_time']:.8g};" in text
        assert values["end_time"] == pytest.approx(
            12.0 * values["t_rev"], rel=1e-12
        )


def test_twins_differ_only_in_blade_keys(cfg):
    case_dir = PACKAGE / "case"
    alm = generate_case.render_fv_options(
        cfg, "7", "coarse", case_dir, generate_case.ALM_ELEMENT
    )
    asm = generate_case.render_fv_options(
        cfg, "7", "coarse", case_dir, generate_case.ASM_ELEMENT
    )
    assert "elementType actuatorLineElement;" in alm
    assert "elementType actuatorSurfaceElement;" in asm
    assert "nChordwise 5;" in asm
    assert "nChordwise" not in alm

    def stripped(text):
        return [
            line
            for line in text.splitlines()
            if not any(key in line for key in BLADE_KEYS)
        ]

    assert stripped(alm) == stripped(asm)

    committed_alm = (case_dir / "system" / "fvOptions.ALM").read_text()
    committed_asm = (case_dir / "system" / "fvOptions.ASM").read_text()
    assert stripped(committed_alm) == stripped(committed_asm)
    assert "elementType actuatorLineElement;" in committed_alm
    assert "nChordwise 5;" in committed_asm


def test_per_speed_measured_tsr(cfg):
    assert_tsr_matches_experiment(cfg)
    for speed, tsr in FROZEN_TSR.items():
        entry = select_speed(cfg, speed)
        assert entry["tsr"] == tsr
        text = generate_case.render_fv_options(
            cfg, speed, "coarse", PACKAGE / "case", generate_case.ALM_ELEMENT
        )
        assert f"tipSpeedRatio {tsr:.8g};" in text
        if speed != "7":
            # the scaffold's fixed 5.3894 must not be used globally
            assert abs(entry["tsr"] - 5.3894) > 0.1
    assert select_speed(cfg, "7", "S")["tsr"] == 5.407


def test_tip_displacement_constraint(cfg):
    for speed in speed_keys(cfg):
        for mesh in EXPECTED_CELLS:
            assert_tip_constraint(cfg, speed, mesh)
            assert tip_displacement(cfg, speed, mesh) < hub_cell_size(cfg, mesh)
    displacement = tip_displacement(cfg, "7", "coarse")
    assert displacement == pytest.approx(0.3028, abs=5e-4)
    assert displacement < EXPECTED_HUB_CELL["coarse"]

    invalid = copy.deepcopy(cfg)
    invalid["mesh"]["resolutions"]["coarse"]["delta_t"] = 0.03
    with pytest.raises(ValueError, match="0.03"):
        validate_config(invalid)


def test_element_data_convention(cfg):
    blade = blade_element_rows(cfg)
    assert len(blade) == 26
    assert blade[0][1] == pytest.approx(0.5083)
    assert blade[-1][1] == pytest.approx(5.029)
    for station, row in zip(
        [(0.5083, 0.0), (1.0085, 6.7), (5.029, -1.815)], (blade[0], blade[3], blade[-1])
    ):
        radius, twist = station
        assert row[1] == pytest.approx(radius)
        assert row[5] == pytest.approx(-(twist + 3.0))
    hub = hub_element_rows()
    assert HUB_DIAMETER == pytest.approx(1.0166)
    assert hub[0][1] == pytest.approx(0.5083)
    assert hub[1][1] == pytest.approx(-0.5083)


def test_stage_matrix(cfg):
    stages = cfg["stages"]
    executing = [name for name, stage in stages.items() if stage.get("executes")]
    assert executing == ["stage0"]
    stage0 = stages["stage0"]
    assert stage0["queue"] == "sequana_cpu_dev"
    assert stage0["max_revolutions"] <= 0.3
    assert stage0["speeds"] == [7]
    assert stage0["meshes"] == ["coarse", "fine"]

    assert stages["stage1"]["meshes"] == ["coarse", "fine"]
    assert stages["stage1"]["speeds"] == [7]
    assert sorted(stages["stage2"]["speeds"]) == [10, 13, 15, 25]
    assert stages["stage2"]["sequence_s"] is True
    stage3 = stages["stage3"]
    assert stage3["optional"] is True
    assert stage3["iddes"] is True
    assert stage3["n_chordwise"] == [1, 3, 5]
    assert stage3["n_chordwise_mesh"] == "fine"
    assert "ultra" in stage3["meshes"]
    for stage in stages.values():
        assert 20.0 not in [float(value) for value in stage.get("speeds", [])]

    invalid = copy.deepcopy(cfg)
    invalid["stages"]["stage2"]["speeds"] = [10, 20]
    with pytest.raises(ValueError, match="20 m/s"):
        validate_config(invalid)


def test_decomposition_rendered(cfg):
    text = generate_case.render_decompose_par(cfg)
    assert "numberOfSubdomains 48;" in text


def test_generated_case_is_current():
    result = subprocess.run(
        [sys.executable, str(PACKAGE / "tools" / "generate_case.py"), "--check"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
