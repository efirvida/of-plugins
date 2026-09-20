# SPDX-License-Identifier: GPL-3.0-or-later
"""Pure-Python checks for the committed nacelle geometry pipeline (W2).

Covers the committed artefacts under `turbinesFoam/geometry/`: the binary STL
(shape, closure, orientation, analytic area/volume), the metadata JSON
(per-node normals/areas, reserved blade fields, input/stl hashes), the
provenance and README documentation, and the generator's CLI contract.

The regeneration checks shell out to the pinned gmsh through
`makeGeometry.py --check` (rather than importing gmsh) and **skip** with an
explicit reason when gmsh is missing, not runnable, or not the pinned version.
No OpenFOAM is required, so this module is never added to
`conftest.SOLVER_DRIVEN`.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
import struct
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

TESTS_DIR = Path(__file__).resolve().parent
TURBINESFOAM = TESTS_DIR.parent
GEOMETRY = TURBINESFOAM / "geometry"
SRC = GEOMETRY / "src"
GENERATOR = SRC / "makeGeometry.py"
STL = GEOMETRY / "stl" / "nacelle.stl"
METADATA = GEOMETRY / "metadata" / "nacelle.json"
PROVENANCE = GEOMETRY / "PROVENANCE.md"
README = GEOMETRY / "README.md"

sys.path.insert(0, str(SRC))

from makeGeometry import (  # noqa: E402
    GMSH_VERSION,
    geo_parameters,
    input_sha256,
)

# canonical S1 geometry (src/nacelle.geo)
R = 1.0
CYLINDER_RATIO = 6.0
BODY_LENGTH = (CYLINDER_RATIO + 1.0) * R  # nose at -7R, downstream end at 0

# the paper's surface resolution (arXiv:1702.02108v4 Sec. 4.1) and the accepted band
N_TRIANGLES_PAPER = 2652
TRIANGLE_TOLERANCE = 0.05

ANALYTIC_AREA = 15.0 * math.pi * R * R  # hemisphere 2 + lateral 12 + cap 1
ANALYTIC_VOLUME = (2.0 / 3.0) * math.pi * R**3 + math.pi * R * R * CYLINDER_RATIO * R


# --------------------------------------------------------------------------
# independent STL reader: the artefacts must be readable without the generator
# --------------------------------------------------------------------------
def read_binary_stl(path: Path):
    """Return `[(facet_normal, v0, v1, v2), ...]` from a binary STL."""
    data = path.read_bytes()
    assert len(data) >= 84, f"{path}: too short for a binary STL"
    count = struct.unpack_from("<I", data, 80)[0]
    assert len(data) == 84 + 50 * count, (
        f"{path}: not a binary STL "
        f"({len(data)} bytes for {count} triangles, expected {84 + 50 * count})"
    )
    triangles = []
    offset = 84
    for _ in range(count):
        values = struct.unpack_from("<12f", data, offset)
        offset += 50
        triangles.append(
            (values[0:3], values[3:6], values[6:9], values[9:12])
        )
    return triangles


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _norm(a):
    return math.sqrt(_dot(a, a))


def triangle_area(triangle):
    _, v0, v1, v2 = triangle
    return 0.5 * _norm(_cross(_sub(v1, v0), _sub(v2, v0)))


def signed_volume(triangles):
    total = 0.0
    for _, v0, v1, v2 in triangles:
        total += _dot(v0, _cross(v1, v2)) / 6.0
    return total


def load_metadata():
    return json.loads(METADATA.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# --------------------------------------------------------------------------
# committed artefacts (always run, no gmsh, no OpenFOAM)
# --------------------------------------------------------------------------
def test_stl_is_binary_with_the_paper_triangle_count():
    triangles = read_binary_stl(STL)
    assert len(triangles) == load_metadata()["n_triangles"]
    assert abs(len(triangles) - N_TRIANGLES_PAPER) <= TRIANGLE_TOLERANCE * N_TRIANGLES_PAPER


def test_stl_geometry_is_the_hemisphere_plus_six_radius_cylinder():
    triangles = read_binary_stl(STL)
    xs = [vertex[0] for triangle in triangles for vertex in triangle[1:]]
    ys = [vertex[1] for triangle in triangles for vertex in triangle[1:]]
    zs = [vertex[2] for triangle in triangles for vertex in triangle[1:]]

    # the downstream end (flat cap) sits at the streamwise origin, the nose upstream
    assert max(xs) == pytest.approx(0.0, abs=1e-6)
    assert min(xs) == pytest.approx(-BODY_LENGTH, abs=0.01 * R)
    # a radius-R body of revolution: nothing outside R, and the surface reaches R
    radius2 = [y * y + z * z for y, z in zip(ys, zs)]
    assert max(radius2) <= R * R + 1e-6
    assert math.sqrt(max(radius2)) == pytest.approx(R, abs=1e-3)
    for axis in (ys, zs):
        assert max(abs(value) for value in axis) <= R + 1e-6
        assert max(abs(value) for value in axis) == pytest.approx(R, abs=1e-3)


def test_surface_is_closed_and_outward_oriented():
    triangles = read_binary_stl(STL)
    # every edge of a closed surface is shared by exactly two triangles
    edges: dict[tuple, int] = {}
    for _, v0, v1, v2 in triangles:
        for first, second in ((v0, v1), (v1, v2), (v2, v0)):
            key = tuple(sorted((first, second)))
            edges[key] = edges.get(key, 0) + 1
    assert len(edges) == 3 * len(triangles) // 2
    assert max(edges.values()) == 2

    # facet normals agree with the winding, and the winding is outward
    for facet_normal, v0, v1, v2 in triangles:
        assert _dot(_cross(_sub(v1, v0), _sub(v2, v0)), facet_normal) > 0.0
    volume = signed_volume(triangles)
    assert volume > 0.0
    assert volume == pytest.approx(ANALYTIC_VOLUME, rel=0.01)


def test_total_area_matches_the_analytic_nacelle_surface():
    triangles = read_binary_stl(STL)
    area = sum(triangle_area(triangle) for triangle in triangles)
    # an inscribed polyhedron is slightly smaller than the smooth surface
    assert area < ANALYTIC_AREA
    assert area == pytest.approx(ANALYTIC_AREA, rel=0.01)


def test_stl_triangle_order_is_canonical():
    """The committed order must not depend on gmsh's (environment-dependent) order.

    The generator re-serializes the mesh: each triangle cyclically rotated to its
    lexicographically smallest form, triangles sorted.  This keeps `--check`
    byte-identical across environments, so the invariant is asserted here on the
    committed artefact.
    """
    triangles = read_binary_stl(STL)
    ordered = []
    for _, v0, v1, v2 in triangles:
        ordered.append(min((v0, v1, v2), (v1, v2, v0), (v2, v0, v1)))
    assert ordered == sorted(ordered)
    assert len(set(ordered)) == len(ordered)  # no duplicate triangles
    assert STL.read_bytes()[:6] == b"turbin"  # fixed header, not gmsh's


def test_metadata_schema_and_node_records():
    metadata = load_metadata()
    for key in (
        "component",
        "format_version",
        "generator",
        "input_sha256",
        "units",
        "n_triangles",
        "parameters",
        "stl_sha256",
        "nodes",
        "_reserved",
    ):
        assert key in metadata, f"metadata is missing '{key}'"

    assert metadata["component"] == "nacelle"
    assert metadata["format_version"] == 1
    assert metadata["units"] == "m"
    assert metadata["generator"]["script"] == "src/nacelle.geo"
    assert metadata["generator"]["gmsh_version"] == GMSH_VERSION
    assert metadata["_reserved"] == {
        "blade": {"radial_station": None, "chord_fraction": None}
    }

    triangles = read_binary_stl(STL)
    nodes = metadata["nodes"]
    assert metadata["n_triangles"] == len(nodes) == len(triangles)

    for position, node in enumerate(nodes):
        assert node["index"] == position
        assert len(node["normal"]) == 3
        assert _norm(tuple(node["normal"])) == pytest.approx(1.0, abs=1e-5)
        assert node["area"] > 0.0
        # the metadata describes the same triangles as the committed STL
        facet_normal = triangles[position][0]
        assert _dot(tuple(node["normal"]), facet_normal) == pytest.approx(1.0, abs=5e-4)
        assert node["area"] == pytest.approx(triangle_area(triangles[position]), rel=1e-5)


def test_metadata_hashes_trace_inputs_and_artefacts():
    metadata = load_metadata()
    geo_bytes = (SRC / "nacelle.geo").read_bytes()
    parameters = geo_parameters(geo_bytes.decode("utf-8"))

    assert parameters == {
        "R": R,
        "cylinderRatio": CYLINDER_RATIO,
        "meshSizeRatio": metadata["parameters"]["meshSizeRatio"],
    }
    assert metadata["parameters"] == parameters
    assert metadata["input_sha256"] == input_sha256(geo_bytes, parameters)
    assert metadata["stl_sha256"] == sha256(STL)


def test_provenance_records_sources_generator_and_hashes():
    text = PROVENANCE.read_text(encoding="utf-8")
    metadata = load_metadata()

    # source references: the nacelle (generated) and the deferred MEXICO blades
    assert "1702.02108v4" in text
    assert "Sec. 4.1" in text
    assert "Fig. 3" in text
    assert "MEXICO" in text
    assert "cylinderRatio" in text or "6R" in text
    # generator + environment requirements
    assert GMSH_VERSION in text
    assert "makeGeometry.py" in text
    assert "glu" in text
    # input and artefact hashes
    assert metadata["input_sha256"] in text
    assert metadata["stl_sha256"] in text
    assert metadata["n_triangles"] == int(metadata["n_triangles"])
    assert str(metadata["n_triangles"]) in text


def test_readme_documents_usage_and_s2_deferral():
    text = README.read_text(encoding="utf-8")
    for token in (
        "makeGeometry.py",
        "--check",
        "--component",
        "--gmsh",
        "nacelle",
        "blade0",
        "S2",
        "input_sha256",
        "_reserved",
        "PROVENANCE.md",
    ):
        assert token in text, f"geometry/README.md does not document '{token}'"


# --------------------------------------------------------------------------
# generator CLI (no gmsh needed for the argument handling)
# --------------------------------------------------------------------------
def test_generator_rejects_deferred_blade_components():
    completed = subprocess.run(
        [sys.executable, str(GENERATOR), "--component", "blade0"],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert completed.returncode != 0
    assert "S2" in completed.stderr
    assert "blade0" in completed.stderr


def test_generator_rejects_unknown_components():
    completed = subprocess.run(
        [sys.executable, str(GENERATOR), "--component", "bogus"],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert completed.returncode != 0
    assert "bogus" in completed.stderr


# --------------------------------------------------------------------------
# regeneration through the pinned gmsh (skip cleanly when unusable)
# --------------------------------------------------------------------------
GMSH_SKIP = (
    "the pinned gmsh {version} is not runnable in this environment "
    "(geometry/PROVENANCE.md: pip gmsh {version} plus the system GL libraries, "
    "e.g. `module load glu/9.0.2_gnu`); the committed artefacts were not "
    "re-verified"
).format(version=GMSH_VERSION)


def gmsh_candidates():
    """Candidate gmsh executables: $GMSH, PATH, then the pytest interpreter's bin."""
    candidates = []
    if os.environ.get("GMSH"):
        candidates.append(os.environ["GMSH"])
    found = shutil.which("gmsh")
    if found:
        candidates.append(found)
    sibling = Path(sys.executable).resolve().parent / "gmsh"
    if sibling.exists():
        candidates.append(str(sibling))
    return candidates


def run_check(script: Path, gmsh: str, cwd: Path | None = None):
    return subprocess.run(
        [sys.executable, str(script), "--check", "--gmsh", gmsh],
        capture_output=True,
        text=True,
        cwd=str(cwd) if cwd else None,
        timeout=900,
    )


@pytest.fixture(scope="module")
def clean_check():
    """`makeGeometry.py --check` on the committed tree, or a clean skip."""
    candidates = gmsh_candidates()
    if not candidates:
        pytest.skip(GMSH_SKIP + " (no gmsh executable found)")
    completed = run_check(GENERATOR, candidates[0])
    if completed.returncode == 2:
        detail = (completed.stderr or completed.stdout).strip().splitlines()
        pytest.skip(GMSH_SKIP + (" — " + detail[0] if detail else ""))
    return SimpleNamespace(gmsh=candidates[0], completed=completed)


def test_check_regeneration_is_clean(clean_check):
    assert clean_check.completed.returncode == 0, (
        clean_check.completed.stdout + clean_check.completed.stderr
    )
    assert "up to date" in clean_check.completed.stdout


def test_check_detects_a_stale_stl_without_touching_the_tree(clean_check, tmp_path):
    copy = tmp_path / "geometry"
    shutil.copytree(GEOMETRY, copy, ignore=shutil.ignore_patterns("__pycache__"))

    stale_stl = copy / "stl" / "nacelle.stl"
    data = bytearray(stale_stl.read_bytes())
    data[84 + 12] ^= 0x01  # flip one byte inside the first vertex
    stale_stl.write_bytes(bytes(data))

    def snapshot():
        return {
            path: sha256(path)
            for path in sorted(copy.rglob("*"))
            if path.is_file()
        }

    before = snapshot()
    committed = sha256(STL)

    completed = run_check(copy / "src" / "makeGeometry.py", clean_check.gmsh, cwd=copy)
    assert completed.returncode == 1, completed.stdout + completed.stderr
    assert "stale STL" in completed.stderr
    # --check is non-destructive: neither the copied nor the committed tree changed
    assert snapshot() == before
    assert sha256(STL) == committed


def test_check_detects_a_stale_provenance(clean_check, tmp_path):
    copy = tmp_path / "geometry"
    shutil.copytree(GEOMETRY, copy, ignore=shutil.ignore_patterns("__pycache__"))

    provenance = copy / "PROVENANCE.md"
    text = provenance.read_text(encoding="utf-8")
    metadata = load_metadata()
    provenance.write_text(
        text.replace(metadata["input_sha256"], "0" * 64), encoding="utf-8"
    )
    before = provenance.read_bytes()

    completed = run_check(copy / "src" / "makeGeometry.py", clean_check.gmsh, cwd=copy)
    assert completed.returncode == 1, completed.stdout + completed.stderr
    assert "stale PROVENANCE.md" in completed.stderr
    assert provenance.read_bytes() == before
