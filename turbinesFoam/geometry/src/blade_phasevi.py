#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Pure-Python structured loft for the `phaseVI_blade` geometry component.

The blade surface is built from the committed inputs only — the 26-station
Phase VI blade geometry (`validation/phaseVI/data/geometry/phaseVI_blade.csv`,
NREL/TP-500-29955 Table A-1) and the committed S809 design coordinates
(`validation/phaseVI/data/s809/s809_somers_nlr.csv`, NREL/SR-440-6918
Table 2) — by a structured loft with a fixed profile point count, a fixed
triangulation and no gmsh involvement (design D10/§7.3).

Sections
--------
Exactly the 26 committed stations, in order.  Every section is a closed
60-point profile in section-local coordinates `(u, w)`: `u` is the LE-based
chord fraction (0 at the leading edge, 1 at the trailing edge) and `w` the
S809 ordinate `z/c`.

* `r <= 0.8835 m` (cylinder stations): a circle of radius `chord/2` centred on
  the chord line at mid-chord, points placed at the same fixed perimeter
  fractions and in the same order as the S809 polygon, so the rings correspond
  1:1 (design §7.3).
* `r >= 1.0085 m` (airfoil stations): the S809 polygon at `(u, w) = (x/c, z/c)`.
* the single transition segment `0.8835 -> 1.0085 m` lofts circle to S809.

The closed profile is the 31 upper-surface points (LE -> TE, ending at the
shared trailing edge) followed by the 29 lower-surface points back towards the
leading edge (the trailing-edge duplicate dropped).  The data contain no
explicit leading-edge point, so the polygon closes across the small blunt
leading-edge gap; no coordinate is invented (see the dataset PROVENANCE.md).

Generation frame (design §4.4)
------------------------------
`+z` is the outward blade span (root -> tip), `+y` the reference chord
direction (trailing edge -> leading edge before twist/pitch) and
`+x = y x z` the section normal (the suction side for the Phase VI/twin
orientation).  A section point is placed as

    y = (0.25 - u) * chord,   x = w * chord,   z = r

and rotated about `+z` by `sigma * (-(twist_report + pitch_deg))` with
`sigma = -1` for Phase VI (the element pitch axis is radially inward), i.e. by
`+(twist_report + PITCH_DEG)`.  The quarter-chord point (`u = 0.25`) then lies
on the radial reference line, exactly as the element placement does.

`PITCH_DEG = 3.0` is the declared Phase VI pitch; `test_blade_data.py` asserts
it against `config/case.yaml` (`turbine.pitch_deg`).

Output
------
One wetted-surface binary STL (25 segments x 60 points = 3000 triangles, fixed
outward winding, no root or tip caps) written through the pipeline's
`canonical_triangles` / `write_binary_stl` helpers with the per-component
header, plus the §4.6 metadata JSON with a populated `radial_station` (the
triangle centroid's span coordinate, m) and `chord_fraction` (the mean of the
three vertices' construction-time LE-based chord fractions) for every node.
Regeneration from the committed inputs is byte-identical: all arithmetic is
deterministic Python over committed inputs and the triangle order is
canonicalized.

Winding note: the design's loft pattern is annotated "(outward)".  The closed
profile is clockwise in the generation frame (measured signed area < 0), so
the outward triangles are `(a, c, b)` and `(b, c, d)` for the ring quad
`a = ring[i][k], b = ring[i][k+1], c = ring[i+1][k], d = ring[i+1][k+1]`;
the literal `(a, b, c)` / `(b, d, c)` pattern yields inward normals (negative
closed-shell volume) for this orientation.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import struct
from dataclasses import dataclass
from pathlib import Path

# geometry/src -> geometry -> turbinesFoam
HERE = Path(__file__).resolve().parent
GEOMETRY = HERE.parent
PACKAGE = GEOMETRY.parent

#: Committed inputs (read-only); the generator consumes these and nothing else.
BLADE_CSV = PACKAGE / "validation/phaseVI/data/geometry/phaseVI_blade.csv"
S809_CSV = PACKAGE / "validation/phaseVI/data/s809/s809_somers_nlr.csv"

#: Registry name, metadata schema version and units (design §4.6).
COMPONENT = "phaseVI_blade"
FORMAT_VERSION = 1
UNITS = "m"
#: Per-component 80-byte binary STL header (the nacelle header bytes are
#: unchanged by the registry change).
STL_HEADER = "turbinesFoam geometry pipeline - phaseVI_blade (S809, wetted)"

#: Declared Phase VI pitch angle, degrees; asserted against
#: `config/case.yaml` (`turbine.pitch_deg`) by `tests/test_blade_data.py`.
PITCH_DEG = 3.0
#: Phase VI/twin orientation: the element pitch axis is radially inward, so the
#: generation-frame rotation about `+z` is `sigma * (-(twist + pitch))`.
SIGMA = -1.0

#: Section rule cut radii: circle at or below, S809 at or above (design §7.3).
CYLINDER_RADIUS = 0.8835
#: Airfoil sections start at this station; the gap is the transition segment.
AIRFOIL_RADIUS = 1.0085
#: Fixed profile point count (31 upper + 29 lower) and station count.
N_PROFILE_POINTS = 60
N_STATIONS = 26

#: Canonical parameter block folded into `input_sha256` together with the
#: committed input bytes (design §4.6: "input bytes + canonical parameter
#: block").
PARAMETERS = {
    "airfoil_radius_m": AIRFOIL_RADIUS,
    "cylinder_radius_m": CYLINDER_RADIUS,
    "n_profile_points": N_PROFILE_POINTS,
    "n_stations": N_STATIONS,
    "pitch_deg": PITCH_DEG,
    "sigma": SIGMA,
}


@dataclass(frozen=True)
class Station:
    """One committed radial station (`phaseVI_blade.csv` row)."""

    radius: float
    chord: float
    twist_deg: float
    chord_mount: float


@dataclass(frozen=True)
class Triangle:
    """One lofted triangle with its construction-time node metadata."""

    v0: tuple[float, float, float]
    v1: tuple[float, float, float]
    v2: tuple[float, float, float]
    radial_station: float
    chord_fraction: float

    @property
    def vertices(self) -> tuple:
        return (self.v0, self.v1, self.v2)


# --------------------------------------------------------------------------
# committed inputs
# --------------------------------------------------------------------------
def _data_rows(path: Path) -> list[str]:
    """CSV lines without the `#` comment/header lines."""
    return [
        line
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]


def load_stations(path: Path = BLADE_CSV) -> tuple[Station, ...]:
    """Parse the committed 26-station blade geometry (strictly increasing r)."""
    stations = []
    for row in csv.reader(_data_rows(path)):
        if len(row) != 4:
            raise ValueError(f"{path}: expected 4 columns, got {len(row)}")
        stations.append(Station(*(float(value) for value in row)))
    if len(stations) != N_STATIONS:
        raise ValueError(
            f"{path}: expected {N_STATIONS} radial stations, got {len(stations)}"
        )
    if any(b.radius <= a.radius for a, b in zip(stations, stations[1:])):
        raise ValueError(f"{path}: radial stations must be strictly increasing")
    for station in stations:
        if station.chord <= 0.0:
            raise ValueError(f"{path}: chord must be positive")
    return tuple(stations)


def load_s809(
    path: Path = S809_CSV,
) -> tuple[tuple[tuple[float, float], ...], tuple[tuple[float, float], ...]]:
    """Parse the committed S809 coordinates into `(upper, lower)`, LE -> TE."""
    rows = list(csv.DictReader(_data_rows(path)))
    upper = tuple(
        (float(row["x_over_c"]), float(row["z_over_c"]))
        for row in rows
        if row["surface"] == "upper"
    )
    lower = tuple(
        (float(row["x_over_c"]), float(row["z_over_c"]))
        for row in rows
        if row["surface"] == "lower"
    )
    if not upper or not lower:
        raise ValueError(f"{path}: both surfaces are required")
    if upper[-1] != (1.0, 0.0) or lower[-1] != (1.0, 0.0):
        raise ValueError(f"{path}: the surfaces must share the trailing edge (1, 0)")
    for name, surface in (("upper", upper), ("lower", lower)):
        if any(b[0] <= a[0] for a, b in zip(surface, surface[1:])):
            raise ValueError(f"{path}: {name} x/c must be strictly increasing")
    return upper, lower


# --------------------------------------------------------------------------
# closed profile (shared by every section)
# --------------------------------------------------------------------------
def closed_profile(
    upper: tuple[tuple[float, float], ...],
    lower: tuple[tuple[float, float], ...],
) -> tuple[tuple[float, float], ...]:
    """The closed 60-point profile: upper LE->TE, then lower back towards the LE.

    The trailing edge appears once; the profile closes across the small blunt
    leading-edge gap between the first upper and the first lower point.
    """
    profile = tuple(upper) + tuple(reversed(lower[:-1]))
    if len(profile) != N_PROFILE_POINTS:
        raise ValueError(
            f"expected {N_PROFILE_POINTS} profile points, got {len(profile)}"
        )
    return profile


def perimeter_fractions(
    profile: tuple[tuple[float, float], ...],
) -> tuple[float, ...]:
    """Cumulative perimeter fraction of every profile point, 0 at the upper LE.

    The closing segment (lower leading edge -> upper leading edge) counts
    towards the perimeter, so the last fraction lies just below 1.0.
    """
    count = len(profile)
    lengths = [
        math.hypot(
            profile[(index + 1) % count][0] - profile[index][0],
            profile[(index + 1) % count][1] - profile[index][1],
        )
        for index in range(count)
    ]
    total = sum(lengths)
    if total <= 0.0:
        raise ValueError("the profile has zero perimeter")
    fractions = []
    walked = 0.0
    for length in lengths:
        fractions.append(walked / total)
        walked += length
    return tuple(fractions)


# --------------------------------------------------------------------------
# sections
# --------------------------------------------------------------------------
def section_uv(
    station: Station,
    profile: tuple[tuple[float, float], ...],
    fractions: tuple[float, ...],
) -> tuple[tuple[float, float], ...]:
    """Section-local `(u, w)` points: circle (cylinder) or S809 (airfoil)."""
    if station.radius <= CYLINDER_RADIUS:
        points = []
        for fraction in fractions:
            phi = math.pi - 2.0 * math.pi * fraction
            points.append((0.5 + 0.5 * math.cos(phi), 0.5 * math.sin(phi)))
        return tuple(points)
    if station.radius >= AIRFOIL_RADIUS:
        return profile
    raise ValueError(
        f"radius {station.radius} m lies between the cylinder "
        f"({CYLINDER_RADIUS} m) and airfoil ({AIRFOIL_RADIUS} m) cut radii"
    )


def section_points(
    station: Station,
    uv: tuple[tuple[float, float], ...],
) -> tuple[tuple[float, float, float], ...]:
    """Map `(u, w)` into the generation frame at the station (design §4.4)."""
    angle = math.radians(SIGMA * -(station.twist_deg + PITCH_DEG))
    cos_angle, sin_angle = math.cos(angle), math.sin(angle)
    points = []
    for u, w in uv:
        x = w * station.chord
        y = (0.25 - u) * station.chord
        points.append(
            (
                x * cos_angle - y * sin_angle,
                x * sin_angle + y * cos_angle,
                station.radius,
            )
        )
    return tuple(points)


def blade_triangles(
    stations: tuple[Station, ...],
    profile: tuple[tuple[float, float], ...],
    fractions: tuple[float, ...],
) -> tuple[Triangle, ...]:
    """Loft the 26 rings into `2 x 60 x 25 = 3000` outward triangles."""
    rings = []
    for station in stations:
        uv = section_uv(station, profile, fractions)
        rings.append((uv, section_points(station, uv)))

    triangles = []
    for index in range(len(rings) - 1):
        (uv_i, ring_i), (uv_j, ring_j) = rings[index], rings[index + 1]
        for k in range(N_PROFILE_POINTS):
            nk = (k + 1) % N_PROFILE_POINTS
            a, b = ring_i[k], ring_i[nk]
            c, d = ring_j[k], ring_j[nk]
            ua, ub, uc, ud = uv_i[k][0], uv_i[nk][0], uv_j[k][0], uv_j[nk][0]
            # outward winding for the clockwise profile (module docstring)
            triangles.append(
                Triangle(a, c, b, (a[2] + c[2] + b[2]) / 3.0, (ua + uc + ub) / 3.0)
            )
            triangles.append(
                Triangle(b, c, d, (b[2] + c[2] + d[2]) / 3.0, (ub + uc + ud) / 3.0)
            )
    return tuple(triangles)


# --------------------------------------------------------------------------
# canonical output
# --------------------------------------------------------------------------
def _float32(value: float) -> float:
    """Round to the binary STL's float32 storage (metadata matches the bytes)."""
    return struct.unpack("<f", struct.pack("<f", value))[0]


def _sig6(value: float) -> float:
    """Round to 6 significant digits (the pipeline's metadata convention)."""
    return float(f"{value:.6g}")


def _rounded(triangle: Triangle) -> Triangle:
    """The triangle as the binary STL stores it (float32 vertices)."""
    rounded = tuple(
        tuple(_float32(value) for value in vertex) for vertex in triangle.vertices
    )
    return Triangle(
        rounded[0],
        rounded[1],
        rounded[2],
        (rounded[0][2] + rounded[1][2] + rounded[2][2]) / 3.0,
        triangle.chord_fraction,
    )


def canonical_records(triangles: tuple[Triangle, ...]) -> tuple[Triangle, ...]:
    """Order the records exactly as `makeGeometry.canonical_triangles` does.

    The reused canonicalizer rotates each triangle to its lexicographically
    smallest cyclic form (preserving the winding) and sorts; every rotation of
    a stored triangle maps back to its record, so the metadata order is the
    STL order.
    """
    from makeGeometry import canonical_triangles

    rounded = tuple(_rounded(triangle) for triangle in triangles)
    canonical = canonical_triangles(
        [(tri.v0, tri.v1, tri.v2, None) for tri in rounded]
    )
    by_rotation = {}
    for triangle in rounded:
        for rotation in (
            triangle.vertices,
            (triangle.v1, triangle.v2, triangle.v0),
            (triangle.v2, triangle.v0, triangle.v1),
        ):
            by_rotation[rotation] = triangle
    records = tuple(by_rotation[triple] for triple in canonical)
    if len(records) != len(rounded):
        raise ValueError("the loft produced duplicate triangles")
    return records


def node_metadata(index: int, triangle: Triangle) -> dict:
    """One metadata node: normal/area plus `radial_station`/`chord_fraction`."""
    from makeGeometry import node_entry

    entry = node_entry(index, triangle.vertices)
    entry["radial_station"] = _sig6(triangle.radial_station)
    entry["chord_fraction"] = _sig6(triangle.chord_fraction)
    return entry


def input_sha256(blade_bytes: bytes, s809_bytes: bytes) -> str:
    """sha256 over the committed input bytes and the canonical parameters."""
    payload = (
        blade_bytes
        + b"\0"
        + s809_bytes
        + b"\0"
        + json.dumps(PARAMETERS, sort_keys=True, separators=(",", ":")).encode()
    )
    return hashlib.sha256(payload).hexdigest()


def render_metadata(metadata: dict) -> str:
    """Render the metadata JSON with one node per line (stable, reviewable)."""
    keys = (
        "component",
        "format_version",
        "generator",
        "inputs",
        "input_sha256",
        "units",
        "n_triangles",
        "stl_sha256",
    )
    lines = ["{"]
    for key in keys:
        line = f' "{key}": ' + json.dumps(metadata[key], separators=(", ", ": "))
        lines.append(line + ",")
    lines.append(' "nodes": [')
    nodes = [
        "  " + json.dumps(node, separators=(", ", ": ")) for node in metadata["nodes"]
    ]
    lines.append(",\n".join(nodes))
    lines.append(" ]")
    lines.append("}")
    return "\n".join(lines) + "\n"


def _relative(path: Path) -> str:
    """The path relative to the `turbinesFoam` package root when it is inside."""
    try:
        return str(path.resolve().relative_to(PACKAGE))
    except ValueError:
        return str(path)


def build(
    stl_path: Path,
    metadata_path: Path,
    *,
    blade_csv: Path = BLADE_CSV,
    s809_csv: Path = S809_CSV,
) -> tuple[str, str]:
    """Generate the STL and metadata; return `(metadata text, summary)`."""
    from makeGeometry import stl_statistics, write_binary_stl

    stations = load_stations(blade_csv)
    upper, lower = load_s809(s809_csv)
    profile = closed_profile(upper, lower)
    fractions = perimeter_fractions(profile)
    records = canonical_records(blade_triangles(stations, profile, fractions))

    stl_path.parent.mkdir(parents=True, exist_ok=True)
    metadata_path.parent.mkdir(parents=True, exist_ok=True)
    triangles = [record.vertices for record in records]
    write_binary_stl(stl_path, triangles, STL_HEADER)
    stl_bytes = stl_path.read_bytes()

    metadata = {
        "component": COMPONENT,
        "format_version": FORMAT_VERSION,
        "generator": {"script": "src/makeGeometry.py", "builder": COMPONENT},
        "inputs": {
            "blade_csv": _relative(blade_csv),
            "blade_csv_sha256": hashlib.sha256(blade_csv.read_bytes()).hexdigest(),
            "s809_csv": _relative(s809_csv),
            "s809_csv_sha256": hashlib.sha256(s809_csv.read_bytes()).hexdigest(),
            "pitch_deg": PITCH_DEG,
        },
        "input_sha256": input_sha256(
            blade_csv.read_bytes(), s809_csv.read_bytes()
        ),
        "units": UNITS,
        "n_triangles": len(records),
        "stl_sha256": hashlib.sha256(stl_bytes).hexdigest(),
        "nodes": [
            node_metadata(index, triangle) for index, triangle in enumerate(records)
        ],
    }
    text = render_metadata(metadata)
    metadata_path.write_text(text, encoding="utf-8")

    stats = stl_statistics(triangles)
    summary = (
        f"{COMPONENT}: {len(records)} triangles, "
        f"x [{stats['bounds'][0][0]:.4g}, {stats['bounds'][0][1]:.4g}] "
        f"y [{stats['bounds'][1][0]:.4g}, {stats['bounds'][1][1]:.4g}] "
        f"z [{stats['bounds'][2][0]:.4g}, {stats['bounds'][2][1]:.4g}], "
        f"area {stats['area']:.6g}, wetted surface (no caps)"
    )
    return text, summary
