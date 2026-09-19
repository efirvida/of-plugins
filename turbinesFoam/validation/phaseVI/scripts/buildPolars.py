#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Build the committed S809 polar files for the NREL Phase VI case.

Reads the committed raw table transcriptions under `data/polars/raw/` (from
NREL/TP-500-29955 Tables A-3..A-8 and NREL/TP-442-7817 Tables B1..B4) and
emits:

- `S809_OSU_Re1M_total.dat`  single-Re baseline (OSU Re 1e6, total drag),
- `S809_multiRe.dat`         multi-Re sensitivity set (OSU clean series).

The baseline uses the total drag coefficient `Cdw` from the wake traverse
where the table reports it and falls back to the pressure drag `Cdp` where the
wake column is blank (the report leaves it blank in the separated region).
Beyond the measured +-20.1/26.1 deg range the file carries a clearly marked
static extension taken from the CSU Table A-3 measurements (mirrored for
negative angles) plus a flat-plate closure at +-120/150/180 deg, so that the
`profileData` lookup never silently returns zero outside the table.

Usage:
    buildPolars.py [--check]

`--check` re-renders in memory and exits 1 when a committed file is missing or
stale, without writing anything.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "polars" / "raw"
POLARS_DIR = ROOT / "data" / "polars"

BASELINE = POLARS_DIR / "S809_OSU_Re1M_total.dat"
MULTI_RE = POLARS_DIR / "S809_multiRe.dat"

# Multi-Re set: nominal Reynolds levels of TP-442-7817 Tables B1..B4.
MULTI_RE_LEVELS = [0.75e6, 1.0e6, 1.25e6, 1.5e6]

# Static extension closure (flat plate); documented as non-measured.
FLAT_PLATE = {
    120.0: (-0.75, 1.5, 0.0),
    150.0: (-0.75, 0.7, 0.0),
    180.0: (0.0, 0.02, 0.0),
}


def read_raw(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"{path}: no data rows")
    rows.sort(key=lambda row: float(row["alpha_deg"]))
    return rows


def drag(row: dict[str, str]) -> float:
    """Total drag where reported, else pressure drag (fallback documented)."""
    cdw = row.get("cdw")
    if cdw:
        return float(cdw)
    cdp = row.get("cdp")
    if cdp:
        return float(cdp)
    raise ValueError("raw table row has neither cdw nor cdp")


def baseline_rows() -> list[tuple[float, float, float, float]]:
    table = read_raw(RAW_DIR / "table_A7.csv")
    rows = [
        (float(row["alpha_deg"]), float(row["cl"]), drag(row), float(row["cm"]))
        for row in table
    ]
    return rows


def extension_rows(core: list[tuple[float, float, float, float]]):
    """CSU Table A-3 measured extension, mirrored for negative alpha."""
    csu = read_raw(RAW_DIR / "table_A3.csv")
    high = [
        (float(row["alpha_deg"]), float(row["cl"]), drag(row), 0.0)
        for row in csu
    ]
    core_high = max(alpha for alpha, _, _, _ in core)
    core_low = min(alpha for alpha, _, _, _ in core)
    rows = [(alpha, cl, cd, cm) for alpha, cl, cd, cm in high if alpha > core_high]
    mirrored = [
        (-alpha, -cl, cd, cm) for alpha, cl, cd, cm in high if alpha > -core_low
    ]
    mirrored.sort()
    rows = mirrored + rows
    for alpha, (cl, cd, cm) in FLAT_PLATE.items():
        rows.append((alpha, cl, cd, cm))
        rows.append((-alpha, -cl, cd, cm))
    rows.sort()
    return rows


def render_single_re(rows) -> str:
    lines = [
        "// S809 polar, OSU Re = 1e6, TOTAL drag (NREL/TP-500-29955 Table A-7;",
        "// cross-checked against NREL/TP-442-7817 Table B2, which is the same",
        "// OSU clean series). `Cdw` (wake traverse) is used where reported; rows",
        "// where the report leaves the wake column blank fall back to `Cdp`",
        "// (pressure drag) -- verified by scripts/buildPolars.py from the raw",
        "// table transcription in data/polars/raw/table_A7.csv.",
        "// Rows with |alpha| beyond the measured OSU range are a STATIC EXTENSION",
        "// from CSU Table A-3 (Re = 0.3e6, mirrored for negative alpha) plus a",
        "// flat-plate closure, added so profileData never reads outside the table.",
        "// (alpha_deg Cl Cd Cm)",
    ]
    for alpha, cl, cd, cm in rows:
        lines.append(f"({alpha:.4g} {cl:.6g} {cd:.6g} {cm:.6g})")
    return "\n".join(lines) + "\n"


def multi_re_tables():
    """Pivot the TP-442-7817 clean series onto a shared angle grid.

    The shared grid is the Re = 1e6 angle list restricted to the range every
    level measures, so no value is extrapolated. Off-grid values are linearly
    interpolated within each level's own measured samples.
    """
    raw = read_raw(RAW_DIR / "table_TP442-7817.csv")
    levels: dict[float, dict[float, tuple[float, float, float, float]]] = {
        level: {} for level in MULTI_RE_LEVELS
    }
    for row in raw:
        level = float(row["re_millions"]) * 1e6
        if level not in levels:
            continue
        alpha = float(row["alpha_deg"])
        levels[level][alpha] = (
            float(row["cl"]),
            drag(row),
            float(row["cm"]),
        )
    reference = sorted(levels[1.0e6])
    low = max(min(samples) for samples in levels.values())
    high = min(max(samples) for samples in levels.values())
    grid = [alpha for alpha in reference if low <= alpha <= high]
    if len(grid) < 10:
        raise ValueError("multi-Re angle grid is too small")
    tables = {level: [] for level in MULTI_RE_LEVELS}
    for alpha in grid:
        for level in MULTI_RE_LEVELS:
            cl, cd, cm = interpolate(levels[level], alpha)
            tables[level].append((cl, cd, cm))
    return grid, tables


def interpolate(samples, x):
    keys = sorted(samples)
    if x <= keys[0]:
        return samples[keys[0]]
    if x >= keys[-1]:
        return samples[keys[-1]]
    for lower, upper in zip(keys, keys[1:]):
        if lower <= x <= upper:
            ratio = (x - lower) / (upper - lower)
            a = samples[lower]
            b = samples[upper]
            return tuple(a[i] + ratio * (b[i] - a[i]) for i in range(3))
    raise AssertionError("unreachable")


def render_multi_re() -> str:
    grid, tables = multi_re_tables()
    lines = [
        "// S809 multi-Re sensitivity set, OSU clean series from NREL/TP-442-7817",
        "// Tables B1..B4 (Re = 0.75/1/1.25/1.5e6). `Cdw` where reported, else",
        "// `Cdp`; values linearly interpolated onto the shared Re = 1e6 grid",
        "// restricted to the range every level measures (no extrapolation).",
        "// Not blended with CSU (Tables A-3..A-5) or DUT (Table A-8).",
        "tableType multiRe;",
        "ReList",
        "(",
        "    " + " ".join(f"{level:.4g}" for level in MULTI_RE_LEVELS),
        ");",
        "clData",
        "(",
    ]
    for index, alpha in enumerate(grid):
        values = " ".join(
            f"{tables[level][index][0]:.6g}" for level in MULTI_RE_LEVELS
        )
        lines.append(f"    ({alpha:.4g} {values})")
    lines.append(");")
    lines.append("cdData")
    lines.append("(")
    for index, alpha in enumerate(grid):
        values = " ".join(
            f"{tables[level][index][1]:.6g}" for level in MULTI_RE_LEVELS
        )
        lines.append(f"    ({alpha:.4g} {values})")
    lines.append(");")
    lines.append("cmData")
    lines.append("(")
    for index, alpha in enumerate(grid):
        values = " ".join(
            f"{tables[level][index][2]:.6g}" for level in MULTI_RE_LEVELS
        )
        lines.append(f"    ({alpha:.4g} {values})")
    lines.append(");")
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail (exit 1) if the committed polar files are missing or stale",
    )
    args = parser.parse_args(argv)

    core = baseline_rows()
    rows = core + extension_rows(core)
    rows.sort(key=lambda row: row[0])
    rendered = {
        BASELINE: render_single_re(rows),
        MULTI_RE: render_multi_re(),
    }

    stale = False
    for path, content in rendered.items():
        if args.check:
            if not path.exists():
                print(f"missing polar file: {path}", file=sys.stderr)
                stale = True
            elif path.read_text(encoding="utf-8") != content:
                print(f"stale polar file: {path}", file=sys.stderr)
                stale = True
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content, encoding="utf-8")
            print(f"wrote {path.relative_to(ROOT)}")
    if args.check:
        return 1 if stale else 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
