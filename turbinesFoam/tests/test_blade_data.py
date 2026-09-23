# SPDX-License-Identifier: GPL-3.0-or-later
"""Pure-Python checks for the committed Phase VI blade geometry (W2).

Covers the committed S809 coordinate dataset and its provenance, the generated
`phaseVI_blade` artefact (binary STL shape/manifoldness, per-node
`radial_station`/`chord_fraction` metadata, input/artefact hashes) and the
builder registry / CLI contract of `geometry/src/makeGeometry.py`.

The blade route is pure Python (no gmsh, no subprocess), so none of these tests
needs gmsh or OpenFOAM; the module is never added to `conftest.SOLVER_DRIVEN`.
The gmsh-gated nacelle regeneration checks live in `test_nacelle_data.py`.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).resolve().parent
TURBINESFOAM = TESTS_DIR.parent
GEOMETRY = TURBINESFOAM / "geometry"
SRC = GEOMETRY / "src"
GENERATOR = SRC / "makeGeometry.py"
STL = GEOMETRY / "stl" / "phaseVI_blade.stl"
METADATA = GEOMETRY / "metadata" / "phaseVI_blade.json"
PACKAGE = TURBINESFOAM / "validation" / "phaseVI"
DATA = PACKAGE / "data"
S809_DIR = DATA / "s809"
S809_CSV = S809_DIR / "s809_somers_nlr.csv"
S809_PROVENANCE = S809_DIR / "PROVENANCE.md"
BLADE_CSV = DATA / "geometry" / "phaseVI_blade.csv"
CONFIG = PACKAGE / "config" / "case.yaml"

sys.path.insert(0, str(SRC))
sys.path.insert(0, str(PACKAGE / "tools"))

from blade_phasevi import (  # noqa: E402
    PITCH_DEG,
    input_sha256 as blade_input_sha256,
    load_s809,
)
from case_config import load_config  # noqa: E402
from makeGeometry import (  # noqa: E402
    BUILDERS,
    COMPONENTS,
    GMSH_BUILDER,
    PYTHON_BUILDER,
    RESERVED_COMPONENTS,
    RETIRED_COMPONENTS,
    component,
)

# Committed bytes are pinned: the S809 digest in its provenance, the blade
# hashes in `geometry/PROVENANCE.md` (and in the blade metadata itself).
S809_SHA256 = "861f3fe5bbda4d5e5d8985e47c5805bada76fc8e1b2157b944ff78d500f4af8b"
S809_PDF_SHA256 = "a7496aad5680d9d4001e36cf8b13d9603ea31956093cded60976dd8fcaffeeb3"
BLADE_INPUT_SHA256 = "7a43047dc0728a04fb3ea5ca9510119242e9aca390bde6895c023140500fe456"
BLADE_STL_SHA256 = "77def499047fe0633e57564247a0fdd1330afbf4148b707db08ec79af0a2ec27"
BLADE_METADATA_SHA256 = "01f75d17e2febb296a7b093f20ff011eb583de52e8f510815f876fcb8169e7dc"

#: A path that can never be executed: the Python route must ignore `--gmsh`
#: (and never consult `PATH`), so any attempt to use it fails loudly.
NO_GMSH = "/nonexistent/gmsh-not-installed"


# --------------------------------------------------------------------------
# independent readers: the artefacts must be readable without the generator
# --------------------------------------------------------------------------
def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
        triangles.append((values[0:3], values[3:6], values[6:9], values[9:12]))
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


def centroid_z(triangle):
    _, v0, v1, v2 = triangle
    return (v0[2] + v1[2] + v2[2]) / 3.0


def _sig6(value: float) -> float:
    """The metadata rounding convention (6 significant digits)."""
    return float(f"{value:.6g}")


def data_rows(path: Path) -> list[list[str]]:
    """CSV rows without the `#` comment/header lines."""
    with path.open(encoding="utf-8") as stream:
        return [row for row in csv.reader(stream) if row and not row[0].startswith("#")]


def _stations() -> list[float]:
    return [float(row[0]) for row in data_rows(BLADE_CSV)]


def tree_snapshot(root: Path) -> dict:
    """sha256 per file under `root` (bytecode caches excluded)."""
    return {
        path.relative_to(root): sha256(path)
        for path in sorted(root.rglob("*"))
        if path.is_file() and "__pycache__" not in path.parts
    }


def run_generator(script: Path, arguments: list[str], cwd: Path, tmp_dir: Path):
    """Run the generator with bytecode caching off and a local temp dir."""
    tmp_dir.mkdir(parents=True, exist_ok=True)
    return subprocess.run(
        [sys.executable, str(script), *arguments],
        capture_output=True,
        text=True,
        cwd=str(cwd),
        timeout=900,
        env={
            **os.environ,
            "PYTHONDONTWRITEBYTECODE": "1",
            "TMPDIR": str(tmp_dir),
        },
    )


# --------------------------------------------------------------------------
# S809 dataset and provenance
# --------------------------------------------------------------------------
def test_s809_provenance_complete():
    text = S809_PROVENANCE.read_text(encoding="utf-8")
    digest = sha256(S809_CSV)

    # the committed dataset is the pinned artefact
    assert digest == S809_SHA256

    # report, DOI and table identifier
    assert "Somers" in text
    assert "NREL/SR-440-6918" in text
    assert "10.2172/437668" in text
    assert "Table 2" in text
    # the verified source PDF, not only the local text dump
    assert S809_PDF_SHA256 in text
    assert "s809_ramsay.txt" in text
    assert "unverified text dump" in text
    # extraction method, tool, version and date
    assert "PyMuPDF" in text
    assert "1.27.2.3" in text
    assert "pdftotext" in text
    assert "2026-09-23" in text
    # per-directory provenance fields (modified requirement)
    for field in (
        "Source PDF sha256",
        "Table",
        "Row selection",
        "Units",
        "Extraction date",
    ):
        assert field in text, f"s809/PROVENANCE.md does not record '{field}'"
    # committed CSV sha256, recomputed from the committed bytes
    assert digest in text
    # independent cross-checks
    assert "TP-500-29955" in text
    assert "TP-442-7817" in text


def test_s809_parses_as_airfoil():
    rows = [
        (row[0], float(row[1]), float(row[2]))
        for row in data_rows(S809_CSV)
        if row[0] in ("upper", "lower")
    ]
    upper = [(x, z) for surface, x, z in rows if surface == "upper"]
    lower = [(x, z) for surface, x, z in rows if surface == "lower"]
    assert (len(upper), len(lower)) == (31, 30)

    # both surfaces run LE -> TE at strictly increasing x/c and share the TE
    for surface in (upper, lower):
        assert surface[-1] == (1.0, 0.0)
        xs = [x for x, _ in surface]
        assert 0.0 <= min(xs) and max(xs) <= 1.0
        assert all(b > a for a, b in zip(xs, xs[1:]))
        assert all(math.isfinite(z) for _, z in surface)
    assert upper[-1] == lower[-1] == (1.0, 0.0)

    # the design table carries no explicit leading-edge point: the profile
    # closes across the documented small blunt-LE gap
    assert upper[0] == pytest.approx((0.00037, 0.00275), abs=1e-12)
    assert lower[0] == pytest.approx((0.00140, -0.00498), abs=1e-12)
    assert abs(upper[0][1] - lower[0][1]) < 0.02

    # the committed CSV header repeats the source and the recorded recipe
    header = S809_CSV.read_text(encoding="utf-8")
    assert "NREL/SR-440-6918" in header
    assert S809_PDF_SHA256 in header

    # the generator's parser agrees with the independent parse
    parsed_upper, parsed_lower = load_s809()
    assert parsed_upper == tuple(upper)
    assert parsed_lower == tuple(lower)


# --------------------------------------------------------------------------
# committed blade artefact (always run, no gmsh, no OpenFOAM)
# --------------------------------------------------------------------------
def test_blade_metadata_complete():
    metadata = json.loads(METADATA.read_text(encoding="utf-8"))
    assert metadata["component"] == "phaseVI_blade"
    assert metadata["format_version"] == 1
    assert metadata["generator"] == {
        "script": "src/makeGeometry.py",
        "builder": "phaseVI_blade",
    }
    assert metadata["units"] == "m"
    # the blade mapping is populated: no reserved null placeholder remains
    assert "_reserved" not in metadata

    inputs = metadata["inputs"]
    assert inputs["blade_csv"] == "validation/phaseVI/data/geometry/phaseVI_blade.csv"
    assert inputs["s809_csv"] == "validation/phaseVI/data/s809/s809_somers_nlr.csv"
    assert inputs["blade_csv_sha256"] == sha256(BLADE_CSV)
    assert inputs["s809_csv_sha256"] == sha256(S809_CSV) == S809_SHA256
    assert inputs["pitch_deg"] == PITCH_DEG
    # the declared Phase VI pitch is the case configuration's pitch
    assert float(load_config(CONFIG)["turbine"]["pitch_deg"]) == PITCH_DEG

    # the declared digests match the committed inputs and artefacts
    assert metadata["input_sha256"] == blade_input_sha256(
        BLADE_CSV.read_bytes(), S809_CSV.read_bytes()
    )
    assert metadata["input_sha256"] == BLADE_INPUT_SHA256
    assert metadata["stl_sha256"] == sha256(STL) == BLADE_STL_SHA256
    assert sha256(METADATA) == BLADE_METADATA_SHA256

    triangles = read_binary_stl(STL)
    nodes = metadata["nodes"]
    assert metadata["n_triangles"] == len(nodes) == len(triangles) == 3000

    fractions = []
    stations = []
    for position, node in enumerate(nodes):
        assert node["index"] == position
        assert len(node["normal"]) == 3
        assert _norm(node["normal"]) == pytest.approx(1.0, abs=1e-5)
        assert node["area"] > 0.0
        # the metadata describes the same triangles as the committed STL
        assert _dot(node["normal"], triangles[position][0]) == pytest.approx(
            1.0, abs=5e-4
        )
        assert node["area"] == pytest.approx(
            triangle_area(triangles[position]), rel=1e-5
        )
        # radial_station is the triangle centroid's span coordinate (m),
        # chord_fraction the mean construction-time LE-based chord fraction
        assert node["radial_station"] == _sig6(centroid_z(triangles[position]))
        assert 0.0 <= node["chord_fraction"] <= 1.0
        stations.append(node["radial_station"])
        fractions.append(node["chord_fraction"])

    # every centroid lies inside the wetted span, and the fractions span the
    # chord (LE-side to TE-side means)
    vertices_z = [
        vertex[2] for triangle in triangles for vertex in triangle[1:]
    ]
    assert min(vertices_z) == pytest.approx(0.5083, abs=1e-6)
    assert max(vertices_z) == pytest.approx(5.029, abs=1e-6)
    assert min(stations) >= min(vertices_z)
    assert max(stations) <= max(vertices_z)
    assert min(fractions) < 0.01
    assert max(fractions) > 0.99


def test_stl_manifold_wetted_only():
    triangles = read_binary_stl(STL)
    assert len(triangles) == 3000

    seen = set()
    for facet_normal, v0, v1, v2 in triangles:
        # no degenerate or duplicate triangles, and the facet normal follows
        # the winding (outward)
        area_vector = _cross(_sub(v1, v0), _sub(v2, v0))
        assert _norm(area_vector) > 0.0
        assert _dot(area_vector, facet_normal) > 0.0
        assert len({v0, v1, v2}) == 3
        key = tuple(sorted((v0, v1, v2)))
        assert key not in seen
        seen.add(key)

    directed: dict[tuple, list] = {}
    for _, v0, v1, v2 in triangles:
        for first, second in ((v0, v1), (v1, v2), (v2, v0)):
            directed.setdefault(tuple(sorted((first, second))), []).append(
                (first, second)
            )

    boundary = []
    for edge, traversals in directed.items():
        assert len(traversals) in (1, 2), edge
        if len(traversals) == 1:
            boundary.append(edge)
        else:
            (a, b), (c, d) = traversals
            assert (a, b) == (d, c), edge  # interior edges are opposed

    # wetted only: the boundary is exactly two 60-edge rings, one per end
    assert len(boundary) == 2 * 60
    neighbours: dict[tuple, set] = {}
    for a, b in boundary:
        neighbours.setdefault(a, set()).add(b)
        neighbours.setdefault(b, set()).add(a)
    assert all(len(incident) == 2 for incident in neighbours.values())

    rings = []
    unseen = set(neighbours)
    while unseen:
        start = min(unseen)
        ring = [start]
        unseen.remove(start)
        previous, current = None, start
        while True:
            candidates = [point for point in neighbours[current] if point != previous]
            assert candidates
            step = candidates[0]
            if step == start:
                break
            assert step in unseen
            unseen.remove(step)
            ring.append(step)
            previous, current = current, step
        rings.append(ring)

    assert sorted(len(ring) for ring in rings) == [60, 60]
    root, tip = sorted(ring[0][2] for ring in rings)
    assert root == pytest.approx(0.5083, abs=1e-6)
    assert tip == pytest.approx(5.029, abs=1e-6)

    # the full wetted span: every committed station plane carries vertices
    stations = _stations()
    assert len(stations) == 26
    planes = {
        round(vertex[2], 4) for triangle in triangles for vertex in triangle[1:]
    }
    assert planes == {round(station, 4) for station in stations}


# --------------------------------------------------------------------------
# deterministic regeneration through the Python route (no gmsh)
# --------------------------------------------------------------------------
def _package_copy(destination: Path) -> Path:
    """The minimal package layout the blade route resolves through `PACKAGE`.

    The blade builder reads its inputs from
    `validation/phaseVI/data/{geometry,s809}/` next to the `geometry/` tree, so
    a `geometry/`-only copy (as the nacelle stale tests use) cannot run
    `--check`; the two input directories are copied alongside.
    """
    shutil.copytree(
        GEOMETRY,
        destination / "geometry",
        ignore=shutil.ignore_patterns("__pycache__"),
    )
    for name in ("geometry", "s809"):
        shutil.copytree(
            DATA / name,
            destination / "validation" / "phaseVI" / "data" / name,
            ignore=shutil.ignore_patterns("__pycache__"),
        )
    return destination


def _corrupt_stl(root: Path) -> None:
    stale = root / "geometry" / "stl" / "phaseVI_blade.stl"
    data = bytearray(stale.read_bytes())
    data[84 + 12] ^= 0x01  # flip one byte inside the first vertex
    stale.write_bytes(bytes(data))


def _corrupt_metadata(root: Path) -> None:
    stale = root / "geometry" / "metadata" / "phaseVI_blade.json"
    stale.write_text(
        stale.read_text(encoding="utf-8").replace(BLADE_STL_SHA256, "0" * 64),
        encoding="utf-8",
    )


def _corrupt_provenance(root: Path) -> None:
    stale = root / "geometry" / "PROVENANCE.md"
    stale.write_text(
        stale.read_text(encoding="utf-8").replace(BLADE_INPUT_SHA256, "0" * 64),
        encoding="utf-8",
    )


def test_check_byte_identical_without_gmsh(tmp_path):
    tmp_dir = tmp_path / "tmp"

    # clean tree, in place: byte-identical regeneration with no gmsh at all
    before = tree_snapshot(GEOMETRY)
    completed = run_generator(
        GENERATOR,
        ["--check", "--component", "phaseVI_blade", "--gmsh", NO_GMSH],
        GEOMETRY,
        tmp_dir,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "up to date" in completed.stdout
    assert BLADE_INPUT_SHA256 in completed.stdout
    assert tree_snapshot(GEOMETRY) == before

    # stale copies: the artifact is detected, the exit code is non-zero and
    # neither the copied tree nor the committed tree is modified
    committed = tree_snapshot(GEOMETRY)
    for label, corrupt in (
        ("stale STL", _corrupt_stl),
        ("stale metadata", _corrupt_metadata),
        ("stale PROVENANCE.md", _corrupt_provenance),
    ):
        copy = _package_copy(tmp_path / label.replace(" ", "_").replace(".", ""))
        corrupt(copy)
        before_copy = tree_snapshot(copy)

        completed = run_generator(
            copy / "geometry" / "src" / "makeGeometry.py",
            ["--check", "--component", "phaseVI_blade", "--gmsh", NO_GMSH],
            copy / "geometry",
            tmp_dir,
        )
        assert completed.returncode == 1, completed.stdout + completed.stderr
        assert label in completed.stderr, completed.stderr
        assert tree_snapshot(copy) == before_copy
        assert tree_snapshot(GEOMETRY) == committed


# --------------------------------------------------------------------------
# builder registry and CLI contract
# --------------------------------------------------------------------------
def test_builder_registry_and_retired_reservation(tmp_path):
    tmp_dir = tmp_path / "tmp"

    # the registry is the single source of routes and names
    assert BUILDERS == {"nacelle": GMSH_BUILDER, "phaseVI_blade": PYTHON_BUILDER}
    assert COMPONENTS == ("nacelle", "phaseVI_blade")
    assert RESERVED_COMPONENTS == ()  # the S2 blade reservation is retired
    assert RETIRED_COMPONENTS == ("blade0", "blade1", "blade2")

    blade = component("phaseVI_blade")
    assert blade.builder == PYTHON_BUILDER
    assert blade.geo is None
    assert blade.stl == STL and blade.metadata == METADATA
    nacelle = component("nacelle")
    assert nacelle.builder == GMSH_BUILDER and nacelle.geo == SRC / "nacelle.geo"

    # phaseVI_blade is accepted through the Python route without any gmsh
    completed = run_generator(
        GENERATOR,
        ["--check", "--component", "phaseVI_blade", "--gmsh", NO_GMSH],
        GEOMETRY,
        tmp_dir,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    assert "up to date" in completed.stdout

    # the help text advertises exactly the registry
    completed = run_generator(GENERATOR, ["--help"], GEOMETRY, tmp_dir)
    assert completed.returncode == 0
    assert "available: nacelle, phaseVI_blade" in " ".join(completed.stdout.split())

    # the retired names answer with the S2 retirement, never as deferred
    for name in RETIRED_COMPONENTS:
        completed = run_generator(GENERATOR, ["--component", name], GEOMETRY, tmp_dir)
        assert completed.returncode == 1
        assert "retired" in completed.stderr
        assert name in completed.stderr
        assert "phaseVI_blade" in completed.stderr
        assert "deferred" not in completed.stderr

    # unknown names are rejected listing the registry
    completed = run_generator(
        GENERATOR, ["--component", "bogus"], GEOMETRY, tmp_dir
    )
    assert completed.returncode == 2
    assert "unknown component 'bogus'" in completed.stderr
    assert "available: nacelle, phaseVI_blade" in completed.stderr

    # the nacelle route still requires the pinned gmsh
    completed = run_generator(
        GENERATOR,
        ["--check", "--component", "nacelle", "--gmsh", NO_GMSH],
        GEOMETRY,
        tmp_dir,
    )
    assert completed.returncode == 2
    assert "cannot execute gmsh" in completed.stderr
