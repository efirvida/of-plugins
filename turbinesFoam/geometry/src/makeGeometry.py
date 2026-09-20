#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Generate the nacelle/hub surface (binary STL + metadata JSON) for the ASM.

The geometry is described by the committed gmsh source `src/nacelle.geo`
(hemisphere + cylinder, Yang & Sotiropoulos, arXiv:1702.02108v4 Sec. 4.1 and
Fig. 3) and is meshed by the pinned gmsh version with fixed options, so that
regeneration from the committed inputs is byte-identical.

Usage:
    makeGeometry.py                        # generate the nacelle STL + metadata
    makeGeometry.py --check                # regenerate and compare; never writes
    makeGeometry.py --component nacelle    # explicit component (S1: nacelle only)
    makeGeometry.py --gmsh /path/to/gmsh   # explicit gmsh executable

`--check` regenerates the component into a temporary directory, byte-compares
the STL and the metadata JSON against the committed artefacts and verifies that
PROVENANCE.md records their sha256 values.  It never modifies the tree.

Exit codes:
    0   success (generated, or `--check` found no stale artefact)
    1   generation failure, or `--check` found a stale/missing artefact
    2   gmsh is missing, not runnable, or not the pinned version (`--check`)

The two sha256 values that PROVENANCE.md must record are printed by every run.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import struct
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # geometry/
SRC = ROOT / "src"
STL_DIR = ROOT / "stl"
METADATA_DIR = ROOT / "metadata"
PROVENANCE = ROOT / "PROVENANCE.md"

#: gmsh version the committed artefacts were generated with (PROVENANCE.md).
GMSH_VERSION = "4.15.2"
#: metadata schema version.
FORMAT_VERSION = 1
#: STL units, documented in the metadata.
UNITS = "m"
#: fixed 80-byte binary STL header (no gmsh version: the geometry is what is
#: committed, and the version is recorded in the metadata).
STL_HEADER = "turbinesFoam geometry pipeline - nacelle (hemisphere + 6R cylinder)"


class GenerationError(Exception):
    """A component could not be generated from the committed inputs."""


@dataclass(frozen=True)
class Component:
    """One generated surface: its gmsh source, STL and metadata paths."""

    name: str
    geo: Path
    stl: Path
    metadata: Path

    def __str__(self) -> str:  # pragma: no cover - trivial
        return self.name


def component(name: str) -> Component:
    return Component(
        name=name,
        geo=SRC / f"{name}.geo",
        stl=STL_DIR / f"{name}.stl",
        metadata=METADATA_DIR / f"{name}.json",
    )


#: Components generated in S1.
COMPONENTS = ("nacelle",)
#: Components the layout reserves for S2 (blade surface model); not generated.
RESERVED_COMPONENTS = ("blade0", "blade1", "blade2")

#: `name = <number>;` declarations in the .geo that carry the canonical
#: geometry parameters (derived expressions such as `L = cylinderRatio*R;` do
#: not match, and neither do `Mesh.*` options or commented-out lines).
PARAMETER_RE = re.compile(
    r"^[ \t]*([A-Za-z_]\w*)[ \t]*=[ \t]*"
    r"([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)[ \t]*;[ \t]*(?://.*)?$",
    re.MULTILINE,
)


# --------------------------------------------------------------------------
# small vector helpers (plain tuples, no numpy dependency)
# --------------------------------------------------------------------------
def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a) -> float:
    return (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]) ** 0.5


def _sig6(value: float) -> float:
    """Round to 6 significant digits (within the binary STL's float32 data)."""
    return float(f"{value:.6g}")


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------
def geo_parameters(geo_text: str) -> dict[str, float]:
    """The primary geometry parameters declared in the .geo source."""
    return {name: float(value) for name, value in PARAMETER_RE.findall(geo_text)}


def input_sha256(geo_bytes: bytes, parameters: dict[str, float]) -> str:
    """sha256 of the gmsh source plus its canonical parameter block."""
    payload = (
        geo_bytes
        + b"\0"
        + json.dumps(parameters, sort_keys=True, separators=(",", ":")).encode()
    )
    return hashlib.sha256(payload).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# --------------------------------------------------------------------------
# STL handling
# --------------------------------------------------------------------------
def read_binary_stl(path: Path) -> list[tuple]:
    """Read a binary STL into `[(v0, v1, v2, facet_normal), ...]`."""
    data = path.read_bytes()
    if len(data) < 84:
        raise ValueError(f"{path}: too short to be a binary STL")
    count = struct.unpack_from("<I", data, 80)[0]
    if len(data) != 84 + 50 * count:
        raise ValueError(
            f"{path}: not a binary STL ({len(data)} bytes for {count} triangles)"
        )
    triangles = []
    offset = 84
    for _ in range(count):
        values = struct.unpack_from("<12f", data, offset)
        offset += 50
        triangles.append(
            (
                (values[3], values[4], values[5]),
                (values[6], values[7], values[8]),
                (values[9], values[10], values[11]),
                (values[0], values[1], values[2]),
            )
        )
    return triangles


def node_entry(index: int, vertices: tuple) -> dict:
    """Per-node (per-triangle-centroid) normal and area, as the ASM consumes them."""
    normal, area = triangle_normal_area(*vertices)
    return {
        "index": index,
        "normal": [_sig6(component) for component in normal],
        "area": _sig6(area),
    }


def triangle_normal_area(v0: tuple, v1: tuple, v2: tuple) -> tuple[tuple, float]:
    """Outward unit normal (right-hand rule on the winding) and area."""
    area_vector = _cross(_sub(v1, v0), _sub(v2, v0))
    length = _norm(area_vector)
    if length <= 0.0:
        return (0.0, 0.0, 0.0), 0.0
    return tuple(component / length for component in area_vector), 0.5 * length


def canonical_triangles(triangles: list[tuple]) -> list[tuple]:
    """Order triangles and their vertices deterministically.

    gmsh's element ordering is not stable across environments (the same mesh is
    emitted in a different order), so the generator re-serializes it: each
    triangle is rotated to its lexicographically smallest cyclic form (which
    preserves the winding and therefore the outward normal) and the triangles
    are sorted.  The result depends only on the float32 vertex values, so the
    committed STL is byte-identical everywhere.
    """
    canonical = []
    for v0, v1, v2, _ in triangles:
        canonical.append(min((v0, v1, v2), (v1, v2, v0), (v2, v0, v1)))
    if len(set(canonical)) != len(canonical):
        raise GenerationError("the mesher produced duplicate triangles")
    return sorted(canonical)


def write_binary_stl(path: Path, triangles: list[tuple], header: str) -> None:
    """Write a binary STL with recomputed facet normals."""
    data = bytearray(header.encode("ascii")[:80].ljust(80, b"\0"))
    data += struct.pack("<I", len(triangles))
    for v0, v1, v2 in triangles:
        normal, _ = triangle_normal_area(v0, v1, v2)
        data += struct.pack("<12fH", *normal, *v0, *v1, *v2, 0)
    path.write_bytes(bytes(data))


def stl_statistics(triangles: list[tuple]) -> dict:
    """Bounds, total area and signed volume of a closed triangle surface."""
    xs = [vertex[0] for triangle in triangles for vertex in triangle]
    ys = [vertex[1] for triangle in triangles for vertex in triangle]
    zs = [vertex[2] for triangle in triangles for vertex in triangle]
    area = 0.0
    volume = 0.0
    for v0, v1, v2 in triangles:
        area += 0.5 * _norm(_cross(_sub(v1, v0), _sub(v2, v0)))
        volume += (
            v0[0] * (v1[1] * v2[2] - v1[2] * v2[1])
            - v0[1] * (v1[0] * v2[2] - v1[2] * v2[0])
            + v0[2] * (v1[0] * v2[1] - v1[1] * v2[0])
        ) / 6.0
    return {
        "bounds": (
            (min(xs), max(xs)),
            (min(ys), max(ys)),
            (min(zs), max(zs)),
        ),
        "area": area,
        "volume": volume,
    }


# --------------------------------------------------------------------------
# gmsh
# --------------------------------------------------------------------------
def gmsh_version(gmsh: str) -> tuple[str | None, str]:
    """Return `(version, reason)`; `version` is None when gmsh is unusable."""
    try:
        completed = subprocess.run(
            [gmsh, "--version"],
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
    except OSError as exc:
        return None, f"{gmsh}: cannot execute gmsh ({exc})"
    except subprocess.TimeoutExpired:
        return None, f"{gmsh}: `gmsh --version` timed out"
    text = (completed.stdout or completed.stderr or "").strip()
    if completed.returncode != 0:
        return None, f"{gmsh}: `gmsh --version` exited {completed.returncode}: {text}"
    if not text:
        return None, f"{gmsh}: `gmsh --version` produced no version"
    return text.splitlines()[-1].strip(), ""


def run_gmsh(gmsh: str, geo: Path, output: Path) -> subprocess.CompletedProcess:
    """Mesh `geo` and write the binary STL to `output`."""
    return subprocess.run(
        [gmsh, "-2", "-format", "stl", "-o", str(output), str(geo)],
        capture_output=True,
        text=True,
        timeout=600,
        check=False,
    )


# --------------------------------------------------------------------------
# metadata rendering
# --------------------------------------------------------------------------
def render_metadata(metadata: dict) -> str:
    """Render the metadata JSON with one node per line (stable, reviewable)."""
    keys = (
        "component",
        "format_version",
        "generator",
        "input_sha256",
        "units",
        "n_triangles",
        "parameters",
        "stl_sha256",
    )
    lines = ["{"]
    for key in keys:
        lines.append(f' "{key}": ' + json.dumps(metadata[key], separators=(", ", ": ")))
        lines[-1] += ","
    lines.append(' "nodes": [')
    nodes = [
        "  " + json.dumps(node, separators=(", ", ": ")) for node in metadata["nodes"]
    ]
    lines.append(",\n".join(nodes))
    lines.append(" ],")
    lines.append(
        ' "_reserved": '
        + json.dumps(metadata["_reserved"], separators=(", ", ": "))
    )
    lines.append("}")
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
# generation
# --------------------------------------------------------------------------
def build_component(target: Component, gmsh: str, version: str) -> tuple[str, str]:
    """Generate `target` into the tree; return `(metadata text, summary text)`."""
    geo_bytes = target.geo.read_bytes()
    parameters = geo_parameters(geo_bytes.decode("utf-8"))
    if not parameters:
        raise GenerationError(f"{target.geo}: no primary parameters declared")

    target.stl.parent.mkdir(parents=True, exist_ok=True)
    target.metadata.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="makeGeometry-") as temp_dir:
        generated_stl = Path(temp_dir) / target.stl.name
        completed = run_gmsh(gmsh, target.geo, generated_stl)
        if completed.returncode != 0 or not generated_stl.exists():
            sys.stderr.write((completed.stdout or "") + (completed.stderr or ""))
            raise GenerationError(
                f"gmsh failed to mesh {target.geo} (exit {completed.returncode})"
            )

        triangles = read_binary_stl(generated_stl)
        for index, (v0, v1, v2, facet_normal) in enumerate(triangles):
            computed = _cross(_sub(v1, v0), _sub(v2, v0))
            if (
                _norm(facet_normal) > 0.0
                and computed[0] * facet_normal[0]
                + computed[1] * facet_normal[1]
                + computed[2] * facet_normal[2]
                <= 0.0
            ):
                raise GenerationError(
                    f"{target.geo}: triangle {index} has an inconsistent normal"
                )
        canonical = canonical_triangles(triangles)
        write_binary_stl(generated_stl, canonical, STL_HEADER)

        stl_bytes = generated_stl.read_bytes()
        metadata = {
            "component": target.name,
            "format_version": FORMAT_VERSION,
            "generator": {
                "script": str(target.geo.relative_to(ROOT)),
                "gmsh_version": version,
            },
            "input_sha256": input_sha256(geo_bytes, parameters),
            "units": UNITS,
            "n_triangles": len(canonical),
            "parameters": parameters,
            "stl_sha256": hashlib.sha256(stl_bytes).hexdigest(),
            "nodes": [
                node_entry(index, triangle)
                for index, triangle in enumerate(canonical)
            ],
            "_reserved": {
                "blade": {"radial_station": None, "chord_fraction": None},
            },
        }
        target.stl.write_bytes(stl_bytes)

    text = render_metadata(metadata)
    target.metadata.write_text(text, encoding="utf-8")

    stats = stl_statistics(canonical)
    summary = (
        f"{target.name}: {len(canonical)} triangles, "
        f"x [{stats['bounds'][0][0]:.4g}, {stats['bounds'][0][1]:.4g}] "
        f"y [{stats['bounds'][1][0]:.4g}, {stats['bounds'][1][1]:.4g}] "
        f"z [{stats['bounds'][2][0]:.4g}, {stats['bounds'][2][1]:.4g}], "
        f"area {stats['area']:.6g}, volume {stats['volume']:.6g}"
    )
    return text, summary


def _gmsh_hint() -> None:
    print(
        "hint: the pinned gmsh (pip package, PROVENANCE.md) needs a runnable "
        "system library environment; on SDumont load the `glu` module.",
        file=sys.stderr,
    )


def generate(component_name: str, gmsh: str) -> int:
    version, reason = gmsh_version(gmsh)
    if version is None:
        print(f"error: {reason}", file=sys.stderr)
        _gmsh_hint()
        return 2
    if version != GMSH_VERSION:
        print(
            f"warning: gmsh {version} differs from the pinned {GMSH_VERSION} "
            f"(PROVENANCE.md)",
            file=sys.stderr,
        )

    target = component(component_name)
    try:
        _, summary = build_component(target, gmsh, version)
    except GenerationError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    parameters = geo_parameters(target.geo.read_text(encoding="utf-8"))
    input_digest = input_sha256(target.geo.read_bytes(), parameters)
    print(summary)
    print(f"wrote {target.stl.relative_to(ROOT)}")
    print(f"wrote {target.metadata.relative_to(ROOT)}")
    print(f"input sha256: {input_digest}")
    print(f"stl sha256:   {file_sha256(target.stl)}")
    print("PROVENANCE.md must record both sha256 values above.")
    return 0


def check(component_name: str, gmsh: str) -> int:
    version, reason = gmsh_version(gmsh)
    if version is None:
        print(f"error: {reason}", file=sys.stderr)
        _gmsh_hint()
        return 2
    if version != GMSH_VERSION:
        print(
            f"error: gmsh {version} is not the pinned {GMSH_VERSION} "
            "(PROVENANCE.md); cannot verify deterministic regeneration",
            file=sys.stderr,
        )
        return 2

    target = component(component_name)
    problems: list[str] = []

    geo_text = target.geo.read_text(encoding="utf-8")
    geo_bytes = target.geo.read_bytes()
    parameters = geo_parameters(geo_text)
    expected_input_sha256 = input_sha256(geo_bytes, parameters)

    with tempfile.TemporaryDirectory(prefix="makeGeometry-check-") as temp_dir:
        candidate = Component(
            name=target.name,
            geo=target.geo,
            stl=Path(temp_dir) / target.stl.name,
            metadata=Path(temp_dir) / target.metadata.name,
        )
        try:
            rendered, _ = build_component(candidate, gmsh, version)
        except GenerationError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 1

        if not target.stl.exists():
            problems.append(f"missing STL: {target.stl.relative_to(ROOT)}")
        elif candidate.stl.read_bytes() != target.stl.read_bytes():
            problems.append(
                f"stale STL: {target.stl.relative_to(ROOT)} "
                f"(committed {file_sha256(target.stl)}, regenerated "
                f"{file_sha256(candidate.stl)})"
            )

        if not target.metadata.exists():
            problems.append(f"missing metadata: {target.metadata.relative_to(ROOT)}")
        elif target.metadata.read_text(encoding="utf-8") != rendered:
            problems.append(
                f"stale metadata: {target.metadata.relative_to(ROOT)} "
                "(regenerated content differs)"
            )

        if not PROVENANCE.exists():
            problems.append("missing PROVENANCE.md")
        else:
            provenance = PROVENANCE.read_text(encoding="utf-8")
            for label, value in (
                ("input sha256", expected_input_sha256),
                ("STL sha256", file_sha256(candidate.stl)),
            ):
                if value not in provenance:
                    problems.append(
                        f"stale PROVENANCE.md: it does not record the current "
                        f"{label} {value}"
                    )

    if problems:
        for problem in problems:
            print(f"stale: {problem}", file=sys.stderr)
        return 1
    print(f"{target.name}: committed STL, metadata and PROVENANCE.md are up to date")
    print(f"input sha256: {expected_input_sha256}")
    return 0


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "command",
        nargs="?",
        default="generate",
        choices=("generate",),
        help="generate (default); use --check to verify without writing",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="regenerate to a temporary directory and compare; never writes",
    )
    parser.add_argument(
        "--component",
        default="nacelle",
        help="component to build; S1 provides 'nacelle' "
        f"(S2 reserves {', '.join(RESERVED_COMPONENTS)})",
    )
    parser.add_argument(
        "--gmsh",
        default=None,
        help="gmsh executable (default: the first `gmsh` found on PATH)",
    )
    args = parser.parse_args(argv)

    if args.component in RESERVED_COMPONENTS:
        print(
            f"error: component '{args.component}' is deferred to S2 (blade surface "
            "model); S1 generates only: " + ", ".join(COMPONENTS),
            file=sys.stderr,
        )
        return 1
    if args.component not in COMPONENTS:
        parser.error(
            f"unknown component '{args.component}' "
            f"(available: {', '.join(COMPONENTS)}; "
            f"deferred: {', '.join(RESERVED_COMPONENTS)})"
        )

    gmsh = args.gmsh or shutil.which("gmsh")
    if not gmsh:
        print(
            "error: no gmsh executable found; pass --gmsh PATH (pinned "
            f"{GMSH_VERSION}, see PROVENANCE.md)",
            file=sys.stderr,
        )
        return 2

    if args.check:
        return check(args.component, gmsh)
    return generate(args.component, gmsh)


if __name__ == "__main__":
    raise SystemExit(main())
