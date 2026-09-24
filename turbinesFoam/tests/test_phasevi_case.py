# SPDX-License-Identifier: GPL-3.0-or-later
"""Pure-Python contracts for the NREL Phase VI validation case package.

These tests never require OpenFOAM: they exercise the YAML single source of
truth, the derived mesh arithmetic and the rendered dictionaries in memory
(spec requirements V1-V11, design sections 4-6).
"""

from __future__ import annotations

import copy
import hashlib
import json
import os
import re
import shutil
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
    solver_delta_t,
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
#: Blade-subdictionary keys the three twins are allowed to differ in. `kernel`
#: is rendered only for the Gaussian ablation (cosine is the default and the
#: key is absent); `surfaceGeometry` is rendered for the mesh-backed variant.
BLADE_KEYS = ("elementType", "nChordwise", "surfaceGeometry", "kernel")


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
    # Rotational augmentation: committed default off, paper constants only.
    augmentation = cfg["actuator"]["rotational_augmentation"]
    assert augmentation["active"] is False
    assert augmentation["model"] == "DuSelig"
    assert (augmentation["a"], augmentation["b"], augmentation["d"]) == (1, 1, 1)
    assert isinstance(cfg["actuator"]["end_effects"]["root"], bool)


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


def block(text: str, header: str) -> str:
    """The `{ ... }` block that starts at `header` (balanced braces)."""
    start = text.index(header)
    depth = 0
    for index in range(start, len(text)):
        if text[index] == "{":
            depth += 1
        elif text[index] == "}":
            depth -= 1
            if depth == 0:
                return text[start : index + 1]
    raise AssertionError(f"unterminated block {header!r}")


def test_twins_differ_only_in_blade_keys(cfg):
    case_dir = PACKAGE / "case"

    def render(element_type, **surface):
        return generate_case.render_fv_options(
            cfg, "7", "coarse", case_dir, element_type, **surface
        )

    alm = render(generate_case.ALM_ELEMENT)
    asm = render(generate_case.ASM_ELEMENT)
    asm_mesh = render(
        generate_case.ASM_ELEMENT,
        surface_geometry=generate_case.SURFACE_GEOMETRY,
    )
    gaussian = render(
        generate_case.ASM_ELEMENT,
        surface_geometry=generate_case.SURFACE_GEOMETRY,
        surface_kernel="gaussian",
    )
    assert "elementType actuatorLineElement;" in alm
    assert "elementType actuatorSurfaceElement;" in asm
    assert "nChordwise 5;" in asm
    assert "nChordwise" not in alm
    assert (
        f'surfaceGeometry "{generate_case.SURFACE_GEOMETRY}";' in asm_mesh
    )
    # Cosine is the default: the ablation key is absent unless selected.
    assert "kernel" not in asm_mesh
    assert "kernel gaussian;" in gaussian
    assert "surfaceGeometry" in block(asm_mesh, "blade1")
    # blade2 is `$blade1;` + azimuthalOffset 180, so the surface keys
    # propagate to both identical blades without a second render.
    assert "$blade1;" in block(asm_mesh, "blade2")

    def stripped(text):
        return [
            line
            for line in text.splitlines()
            if not any(key in line for key in BLADE_KEYS)
        ]

    assert stripped(alm) == stripped(asm) == stripped(asm_mesh) == stripped(gaussian)

    committed = {
        name: (case_dir / "system" / f"fvOptions.{name}").read_text()
        for name in ("ALM", "ASM", "ASM-MESH")
    }
    assert stripped(committed["ALM"]) == stripped(committed["ASM"])
    assert stripped(committed["ASM"]) == stripped(committed["ASM-MESH"])
    assert "elementType actuatorLineElement;" in committed["ALM"]
    assert "nChordwise 5;" in committed["ASM"]
    assert (
        f'surfaceGeometry "{generate_case.SURFACE_GEOMETRY}";'
        in committed["ASM-MESH"]
    )
    assert "surfaceGeometry" in block(committed["ASM-MESH"], "blade1")
    assert "$blade1;" in block(committed["ASM-MESH"], "blade2")

    # The augmentation block is present and identical in every twin and every
    # committed file: it is not a blade key, so stripping the blade keys must
    # leave it untouched and identical.
    augmentation_blocks = {
        block(text, "rotationalAugmentation")
        for text in (alm, asm, asm_mesh, gaussian)
    }
    assert len(augmentation_blocks) == 1
    assert "active off;" in augmentation_blocks.pop()
    committed_augmentation = {
        block(text, "rotationalAugmentation") for text in committed.values()
    }
    assert len(committed_augmentation) == 1


def test_rotational_augmentation_rendered(cfg, tmp_path):
    """The switch renders on/off identically across the twins; default off."""
    constants = {"active": True, "model": "DuSelig", "a": 1, "b": 1, "d": 1}

    def render(element_type, **kwargs):
        return generate_case.render_fv_options(
            cfg,
            "7",
            "coarse",
            PACKAGE / "case",
            element_type,
            rotational_augmentation=constants,
            **kwargs,
        )

    alm = block(render(generate_case.ALM_ELEMENT), "rotationalAugmentation")
    asm = block(render(generate_case.ASM_ELEMENT), "rotationalAugmentation")
    asm_mesh = block(
        render(
            generate_case.ASM_ELEMENT,
            surface_geometry=generate_case.SURFACE_GEOMETRY,
        ),
        "rotationalAugmentation",
    )
    assert alm == asm == asm_mesh
    for key in ("active on;", "model DuSelig;", "a 1;", "b 1;", "d 1;"):
        assert key in alm

    # The block must sit at the rotor-coeffs level (mirroring `dynamicStall`),
    # never inside a blade subdict: only the rotor-level block is forwarded by
    # AFTAL together with `rotorRadius`/`rootRadius`, so a blade-level block
    # makes the element skip the correction (the proxy matrix caught this).
    coeffs = block(
        render(generate_case.ALM_ELEMENT), "axialFlowTurbineALSourceCoeffs"
    )
    assert "rotationalAugmentation" in coeffs
    assert "rotationalAugmentation" not in block(coeffs, "blade1")

    # The committed default renders the same block with the switch off.
    committed = (PACKAGE / "case" / "system" / "fvOptions.ALM").read_text()
    default = block(committed, "rotationalAugmentation")
    assert "active off;" in default
    assert "model DuSelig;" in default
    for key in ("a 1;", "b 1;", "d 1;"):
        assert key in default

    # The CLI toggle reaches every rendered twin.
    case_dir = tmp_path / "case"
    result = subprocess.run(
        [
            sys.executable,
            str(PACKAGE / "tools" / "generate_case.py"),
            "--rotational-augmentation", "on",
            "--case-dir", str(case_dir),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    for name in ("ALM", "ASM", "ASM-MESH"):
        rendered = (case_dir / "system" / f"fvOptions.{name}").read_text()
        assert "active on;" in block(rendered, "rotationalAugmentation")


def test_root_effect_ablation_rendered(cfg, tmp_path):
    """The root-effect ablation is render-time; the committed default keeps on."""
    committed = (PACKAGE / "case" / "system" / "fvOptions.ALM").read_text()
    assert "rootEffects on;" in block(committed, "GlauertCoeffs")

    ablation = generate_case.render_fv_options(
        cfg,
        "7",
        "coarse",
        PACKAGE / "case",
        generate_case.ALM_ELEMENT,
        root_effects=False,
    )
    assert "rootEffects off;" in block(ablation, "GlauertCoeffs")
    assert "rootEffects on;" not in ablation
    # The tip effect is not part of the ablation.
    assert "tipEffects on;" in ablation

    case_dir = tmp_path / "case"
    result = subprocess.run(
        [
            sys.executable,
            str(PACKAGE / "tools" / "generate_case.py"),
            "--root-effects", "off",
            "--case-dir", str(case_dir),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    rendered = (case_dir / "system" / "fvOptions.ALM").read_text()
    assert "rootEffects off;" in block(rendered, "GlauertCoeffs")


def test_asm_mesh_selection(cfg, tmp_path):
    """`asm-mesh` is selectable end-to-end and its runs are prepared-only."""
    selected = subprocess.run(
        [
            sys.executable,
            str(PACKAGE / "tools" / "case_config.py"),
            "--select", "7", "coarse", "asm-mesh",
        ],
        capture_output=True,
        text=True,
    )
    assert selected.returncode == 0, selected.stderr
    payload = json.loads(selected.stdout)
    assert payload["model"] == "asm-mesh"
    assert payload["tsr"] == pytest.approx(5.408, abs=1e-9)

    rejected = subprocess.run(
        [
            sys.executable,
            str(PACKAGE / "tools" / "case_config.py"),
            "--select", "7", "coarse", "bogus",
        ],
        capture_output=True,
        text=True,
    )
    assert rejected.returncode == 2
    assert "['alm', 'asm', 'asm-mesh']" in rejected.stderr

    # The runner refuses --submit for the model before any environment check:
    # the ASM-mesh array is prepared-only and must not be submitted.
    refused = subprocess.run(
        [
            "sh", str(PACKAGE / "scripts" / "runPhaseVI.sh"),
            "-m", "asm-mesh", "-u", "7", "-mesh", "coarse", "--submit",
        ],
        capture_output=True,
        text=True,
    )
    assert refused.returncode == 2
    assert "scripts/slurm/asm-mesh.slurm" in refused.stderr
    assert "prepared-only" in refused.stderr

    # A prepare-only invocation (no --run/--submit) renders the run directory
    # `runs/asm-mesh-U7-coarse` without OpenFOAM. `check_environment.sh` is
    # stubbed in a sandbox copy and the shared mesh is pre-created so the
    # runner stops after preparation, exactly like the W3 gate rehearsal.
    sandbox = tmp_path / "sandbox"
    package = sandbox / "validation" / "phaseVI"
    package.parent.mkdir(parents=True)
    # `runs/`, `results/` and any generated `polyMesh` are gitignored
    # working-tree state (the committed package is small); the sandbox only
    # needs the rendered skeleton, the tooling and the data.
    shutil.copytree(
        PACKAGE,
        package,
        ignore=shutil.ignore_patterns(
            "runs", "results", "polyMesh", "__pycache__"
        ),
    )
    (sandbox / "geometry").symlink_to(TURBINESFOAM / "geometry")
    environment = package / "scripts" / "check_environment.sh"
    environment.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    environment.chmod(0o755)
    mesh = package / "runs" / "mesh-coarse" / "constant" / "polyMesh"
    mesh.mkdir(parents=True)
    (mesh / ".keep").write_text("", encoding="utf-8")

    # The runner calls `python3` from PATH; the test's interpreter is the one
    # that must be used (a `python3` shim, so the test does not depend on the
    # invoking environment's PATH order).
    shim = tmp_path / "bin"
    shim.mkdir()
    (shim / "python3").symlink_to(sys.executable)
    # The sandbox must be isolated from an ambient Slurm allocation: the runner
    # enforces --ranks == SLURM_NTASKS when SLURM_NTASKS is set, so inheriting
    # the batch job's SLURM_* would make this test fail inside a compute node.
    env = {
        key: value
        for key, value in os.environ.items()
        if not key.startswith("SLURM_")
    }
    env["PATH"] = str(shim) + os.pathsep + env.get("PATH", "")

    prepared = subprocess.run(
        [
            "sh", str(package / "scripts" / "runPhaseVI.sh"),
            "-m", "asm-mesh", "-u", "7", "-mesh", "coarse",
        ],
        cwd=package,
        capture_output=True,
        text=True,
        env=env,
    )
    assert prepared.returncode == 0, prepared.stderr
    run_dir = package / "runs" / "asm-mesh-U7-coarse"
    assert run_dir.is_dir()
    committed_stl = TURBINESFOAM / "geometry" / "stl" / "phaseVI_blade.stl"
    staged = run_dir / "constant" / "triSurface" / "phaseVI_blade.stl"
    assert staged.read_bytes() == committed_stl.read_bytes()
    payload = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    assert payload["model"] == "asm-mesh"
    assert payload["staged_stl_sha256"] == hashlib.sha256(
        staged.read_bytes()
    ).hexdigest()
    # The installed twin is the committed ASM-MESH render; only the
    # case-relative polars `#include` differs (the run directory sits
    # elsewhere in the tree), exactly as for the committed skeleton render.
    def without_polars_include(text):
        return re.sub(r'#include "[^"]*"', '#include "POLARS"', text)

    installed = (run_dir / "system" / "fvOptions").read_text(encoding="utf-8")
    committed_twin = (
        PACKAGE / "case" / "system" / "fvOptions.ASM-MESH"
    ).read_text(encoding="utf-8")
    assert without_polars_include(installed) == without_polars_include(committed_twin)


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

    # The ASM-mesh third variant is prepared, never executed by the plan.
    asm_mesh = stages["asm-mesh"]
    assert asm_mesh["queue"] == "sequana_cpu_long"
    assert asm_mesh["models"] == ["asm-mesh"]
    assert asm_mesh["meshes"] == ["coarse", "fine"]
    assert asm_mesh["speeds"] == [7]
    assert asm_mesh["executes"] is False

    for stage in stages.values():
        assert 20.0 not in [float(value) for value in stage.get("speeds", [])]

    invalid = copy.deepcopy(cfg)
    invalid["stages"]["stage2"]["speeds"] = [10, 20]
    with pytest.raises(ValueError, match="20 m/s"):
        validate_config(invalid)


def test_decomposition_rendered(cfg):
    text = generate_case.render_decompose_par(cfg)
    assert "numberOfSubdomains 48;" in text
    override = generate_case.render_decompose_par(cfg, ranks=192)
    assert "numberOfSubdomains 192;" in override
    assert "numberOfSubdomains 48;" not in override


def test_iddes_variant(cfg):
    assert cfg["solver"]["iddes_model"] == "kOmegaSSTIDDES"
    iddes = cfg["iddes"]
    assert iddes["delta"] == "IDDESDelta"
    assert iddes["delta_t"] == 0.0025

    les = generate_case.render_turbulence_properties(cfg, "iddes")
    assert "simulationType LES;" in les
    assert "LESModel kOmegaSSTIDDES;" in les
    assert "delta IDDESDelta;" in les
    assert "IDDESDeltaCoeffs" in les
    assert "Cw 0.15;" in les
    assert "RASModel" not in les

    ras = generate_case.render_turbulence_properties(cfg)
    assert "simulationType RAS;" in ras
    assert "RASModel kOmegaSST;" in ras
    assert "LESModel" not in ras

    schemes = generate_case.render_fv_schemes(cfg, "iddes")
    assert "div(phi,U)      Gauss linear;" in schemes
    assert "nRequired       true;" in schemes
    assert "bounded Gauss linearUpwind grad(U)" not in schemes
    urans_schemes = generate_case.render_fv_schemes(cfg)
    assert "div(phi,U)      bounded Gauss linearUpwind grad(U);" in urans_schemes
    assert "nRequired" not in urans_schemes

    for mesh in EXPECTED_CELLS:
        assert solver_delta_t(cfg, mesh, "iddes") == 0.0025
        control = generate_case.render_control_dict(
            cfg, "7", mesh, "production", None, "startTime", solver="iddes"
        )
        assert "deltaT 0.0025;" in control
        assert "adjustTimeStep off;" in control
        assert_tip_constraint(cfg, "7", mesh, "iddes")
        assert tip_displacement(cfg, "7", mesh, "iddes") < hub_cell_size(cfg, mesh)
    assert tip_displacement(cfg, "7", "fine", "iddes") == pytest.approx(
        5.408 * 7.0 * 0.0025, rel=1e-12
    )
    assert kinematics(cfg, "7", "fine", "H", "iddes")["delta_t"] == 0.0025
    assert kinematics(cfg, "7", "fine")["delta_t"] == 0.005

    invalid = copy.deepcopy(cfg)
    invalid["iddes"]["delta"] = "cubeRootVol"
    with pytest.raises(ValueError, match="IDDESDelta"):
        validate_config(invalid)
    invalid = copy.deepcopy(cfg)
    invalid["iddes"]["delta_t"] = 0.0
    with pytest.raises(ValueError, match="delta_t"):
        validate_config(invalid)


def test_n_chordwise_override(cfg):
    asm = generate_case.render_fv_options(
        cfg, "7", "fine", PACKAGE / "case", generate_case.ASM_ELEMENT,
        n_chordwise=1,
    )
    assert "nChordwise 1;" in asm
    assert "nChordwise 5;" not in asm
    # The override reaches the ASM twin only; the ALM twin has no such key.
    alm = generate_case.render_fv_options(
        cfg, "7", "fine", PACKAGE / "case", generate_case.ALM_ELEMENT,
        n_chordwise=1,
    )
    assert "nChordwise" not in alm
    with pytest.raises(ValueError, match="positive"):
        generate_case.render_fv_options(
            cfg, "7", "fine", PACKAGE / "case", generate_case.ASM_ELEMENT,
            n_chordwise=0,
        )


def test_cli_variant_render(cfg, tmp_path):
    """The CLI renders IDDES + nChordwise/ranks overrides in one pass."""
    case_dir = tmp_path / "case"
    result = subprocess.run(
        [
            sys.executable,
            str(PACKAGE / "tools" / "generate_case.py"),
            "--mesh", "fine",
            "--solver", "iddes",
            "--n-chordwise", "1",
            "--ranks", "192",
            "--case-dir", str(case_dir),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr

    turbulence = (case_dir / "constant" / "turbulenceProperties").read_text()
    assert "simulationType LES;" in turbulence
    assert "LESModel kOmegaSSTIDDES;" in turbulence
    assert "delta IDDESDelta;" in turbulence
    control = (case_dir / "system" / "controlDict").read_text()
    assert "deltaT 0.0025;" in control
    assert "nChordwise 1;" in (case_dir / "system" / "fvOptions.ASM").read_text()
    assert "numberOfSubdomains 192;" in (
        case_dir / "system" / "decomposeParDict"
    ).read_text()


def test_generated_case_is_current(cfg):
    # The committed third twin is exactly the in-memory render (so `--check`
    # covers it and the committed file is the documented renderer output).
    rendered = generate_case.render_fv_options(
        cfg,
        "7",
        "coarse",
        PACKAGE / "case",
        generate_case.ASM_ELEMENT,
        surface_geometry=generate_case.SURFACE_GEOMETRY,
    )
    committed = (PACKAGE / "case" / "system" / "fvOptions.ASM-MESH").read_text()
    assert committed == rendered

    # The other two committed twins are their in-memory renders too.
    for name, element_type in (
        ("ALM", generate_case.ALM_ELEMENT),
        ("ASM", generate_case.ASM_ELEMENT),
    ):
        rendered_twin = generate_case.render_fv_options(
            cfg, "7", "coarse", PACKAGE / "case", element_type
        )
        committed_twin = (
            PACKAGE / "case" / "system" / f"fvOptions.{name}"
        ).read_text()
        assert committed_twin == rendered_twin

    result = subprocess.run(
        [sys.executable, str(PACKAGE / "tools" / "generate_case.py"), "--check"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
