# SPDX-License-Identifier: GPL-3.0-or-later
"""Pure-Python contracts for the ASM-mesh STL staging helper (W3.2/W3.5).

`tools/stage_blade_stl.py` is what makes a surface-model run directory
self-contained: it verifies the committed blade STL against the geometry
metadata sha256, copies it to the case-relative
`constant/triSurface/phaseVI_blade.stl` the `fvOptions.ASM-MESH` twin
references, and aborts (exit 3) on a hash mismatch or a missing source.

The helper is pure Python (no gmsh, no OpenFOAM), so none of these tests needs
a runnable solver; the module is never added to `conftest.SOLVER_DRIVEN`.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).resolve().parent
TURBINESFOAM = TESTS_DIR.parent
PACKAGE = TURBINESFOAM / "validation" / "phaseVI"
TOOLS = PACKAGE / "tools"
GEOMETRY = TURBINESFOAM / "geometry"
STL = GEOMETRY / "stl" / "phaseVI_blade.stl"
METADATA = GEOMETRY / "metadata" / "phaseVI_blade.json"
ASM_MESH_TWIN = PACKAGE / "case" / "system" / "fvOptions.ASM-MESH"
sys.path.insert(0, str(TOOLS))

from stage_blade_stl import (  # noqa: E402
    DEFAULT_METADATA,
    DEFAULT_STL,
    EXIT_STAGING,
    STAGED_RELATIVE,
    StageError,
    expected_stl_sha256,
    sha256_file,
    stage_blade_stl,
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_defaults_point_at_the_committed_geometry():
    """The helper's defaults must be the committed STL and metadata."""
    assert DEFAULT_STL.resolve() == STL.resolve()
    assert DEFAULT_METADATA.resolve() == METADATA.resolve()
    assert STL.is_file() and METADATA.is_file()
    digest = expected_stl_sha256(METADATA)
    assert digest == sha256_file(STL) == sha256(STL)
    assert len(digest) == 64


def test_run_dir_layout_is_the_twin_reference():
    """The staged path is exactly the case-relative `surfaceGeometry` value."""
    assert STAGED_RELATIVE.as_posix() == "constant/triSurface/phaseVI_blade.stl"
    twin = ASM_MESH_TWIN.read_text(encoding="utf-8")
    assert f'surfaceGeometry "{STAGED_RELATIVE.as_posix()}";' in twin


def test_stage_into_run_dir(tmp_path):
    run_dir = tmp_path / "asm-mesh-U7-coarse"
    staged = stage_blade_stl(STL, METADATA, run_dir)
    destination = run_dir / STAGED_RELATIVE
    assert staged == sha256(STL)
    assert destination.is_file()
    assert destination.read_bytes() == STL.read_bytes()
    # Staging is idempotent: a prepared run can be re-staged unchanged.
    assert stage_blade_stl(STL, METADATA, run_dir) == staged


def test_stage_aborts_on_mismatch(tmp_path):
    corrupt = tmp_path / "corrupt.stl"
    shutil.copyfile(STL, corrupt)
    with corrupt.open("ab") as stream:
        stream.write(b"stale")
    run_dir = tmp_path / "run"
    with pytest.raises(StageError, match="sha256 mismatch"):
        stage_blade_stl(corrupt, METADATA, run_dir)
    # The mismatch aborts before any copy: no stale surface is left behind.
    assert not (run_dir / STAGED_RELATIVE).exists()


def test_stage_aborts_on_missing_source(tmp_path):
    run_dir = tmp_path / "run"
    with pytest.raises(StageError, match="not found"):
        stage_blade_stl(tmp_path / "missing.stl", METADATA, run_dir)
    assert not (run_dir / STAGED_RELATIVE).exists()


def test_stage_aborts_on_unusable_metadata(tmp_path):
    run_dir = tmp_path / "run"
    with pytest.raises(StageError, match="cannot read blade metadata"):
        stage_blade_stl(STL, tmp_path / "missing.json", run_dir)
    empty = tmp_path / "empty.json"
    empty.write_text(json.dumps({"component": "phaseVI_blade"}), encoding="utf-8")
    with pytest.raises(StageError, match="does not record stl_sha256"):
        stage_blade_stl(STL, empty, run_dir)
    invalid = tmp_path / "invalid.json"
    invalid.write_text("{not json", encoding="utf-8")
    with pytest.raises(StageError, match="not valid JSON"):
        stage_blade_stl(STL, invalid, run_dir)
    assert not (run_dir / STAGED_RELATIVE).exists()


def test_cli_contract(tmp_path):
    run_dir = tmp_path / "run"
    command = [sys.executable, str(TOOLS / "stage_blade_stl.py"), "--run-dir", str(run_dir)]
    happy = subprocess.run(command, capture_output=True, text=True)
    assert happy.returncode == 0, happy.stderr
    # The runner captures the bare hash from stdout for `run.json`.
    assert happy.stdout == f"{sha256(STL)}\n"
    assert "staged" in happy.stderr and STAGED_RELATIVE.as_posix() in happy.stderr

    corrupt = tmp_path / "corrupt.stl"
    shutil.copyfile(STL, corrupt)
    with corrupt.open("ab") as stream:
        stream.write(b"stale")
    broken = subprocess.run(
        command + ["--stl", str(corrupt)], capture_output=True, text=True
    )
    assert broken.returncode == EXIT_STAGING
    assert broken.stderr.startswith("staging error:")
    # The corrupt run must not overwrite the good copy staged above.
    assert (run_dir / STAGED_RELATIVE).read_bytes() == STL.read_bytes()

    bad_run_dir = tmp_path / "corrupt-run"
    broken = subprocess.run(
        [sys.executable, str(TOOLS / "stage_blade_stl.py"), "--run-dir", str(bad_run_dir),
         "--stl", str(corrupt)],
        capture_output=True, text=True,
    )
    assert broken.returncode == EXIT_STAGING
    assert not (bad_run_dir / STAGED_RELATIVE).exists()
