#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Extract the committed Phase VI experimental CSVs from the WDH workbook.

Reads the NWTC `wt_loads_statistics.xls` workbook (public, DOI
10.21947/WDH-DAP/1910052), sheet `ldsmean`, keeps only the canonical 0 deg-yaw
Sequence H and Sequence S rows, and writes the derived performance and spanwise
CSVs under `data/experiment/`.

The workbook is NOT committed (17 MB); this script documents the extraction so
the committed numbers can be traced to the source sheet, rows, and columns.

Usage:
    extractExperimentalData.py --workbook /path/to/wt_loads_statistics.xls [--check]
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT_DIR = ROOT / "data" / "experiment"

SHEET = "ldsmean"

# Frozen column manifest (verified against the sheet header row, 2026-09-19).
COLUMNS = {
    "file_name": 0,
    "windspeed": 2,
    "yaw_token": 3,
    "repetition": 4,
    "wtbaro": 30,
    "wtatemp": 34,
    "yaw_deg": 47,
    "rpm": 52,
    "eaeroth": 54,
    "lsstqcor": 58,
    "rotpow": 61,
    "cn30": 63,
    "ct30": 64,
    "cm30": 65,
    "cn47": 67,
    "ct47": 68,
    "cm47": 69,
    "cn63": 71,
    "ct63": 72,
    "cm63": 73,
    "cn80": 75,
    "ct80": 76,
    "cm80": 77,
    "cn95": 79,
    "ct95": 80,
    "cm95": 81,
}

SPANWISE_COLUMNS = {
    0.30: ("cn30", "ct30", "cm30"),
    0.47: ("cn47", "ct47", "cm47"),
    0.63: ("cn63", "ct63", "cm63"),
    0.80: ("cn80", "ct80", "cm80"),
    0.95: ("cn95", "ct95", "cm95"),
}

# Canonical rows: Sequence H one row per speed, Sequence S 7 m/s repeat.
CANONICAL_ROWS = {
    "H": {
        7.0: "h0700000",
        10.0: "h1000000",
        13.0: "h1300000",
        15.0: "h1500000",
        20.0: "h20m0000",
        25.0: "h2500000",
    },
    "S": {
        7.0: "s0700000",
    },
}

ROTOR_RADIUS = 5.029
MAX_ABS_YAW_DEG = 0.5
RHO_GAS_CONSTANT = 287.058  # J/(kg K), dry air
KELVIN_OFFSET = 273.15

PERFORMANCE_HEADER = [
    "file_name",
    "sequence",
    "wind_speed_m_s",
    "yaw_token",
    "yaw_deg",
    "repetition",
    "rpm",
    "tsr",
    "rotpow_kw",
    "lsstqcor_nm",
    "eaeroth_n",
    "wtbaro_pa",
    "wtatemp_degc",
    "rho_kg_m3",
]
SPANWISE_HEADER = [
    "file_name",
    "sequence",
    "wind_speed_m_s",
    "r_over_R",
    "cn",
    "ct",
    "cm",
]


def load_rows(workbook: Path):
    try:
        import xlrd
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise SystemExit(
            "xlrd is required to read the .xls workbook (pip install xlrd)"
        ) from exc
    book = xlrd.open_workbook(str(workbook), on_demand=True)
    if SHEET not in book.sheet_names():
        raise SystemExit(f"{workbook}: sheet {SHEET!r} not found")
    sheet = book.sheet_by_name(SHEET)
    rows = {}
    for index in range(1, sheet.nrows):
        name = sheet.cell_value(index, COLUMNS["file_name"])
        if isinstance(name, str) and name:
            rows.setdefault(name, []).append(index)
    return sheet, rows


def cell(sheet, row, key):
    return sheet.cell_value(row, COLUMNS[key])


def select_rows(sheet, rows):
    """Resolve the canonical rows and enforce the 0-deg-yaw rule."""
    selected = []
    for sequence, speeds in CANONICAL_ROWS.items():
        for speed, name in sorted(speeds.items()):
            candidates = rows.get(name, [])
            if not candidates:
                raise SystemExit(f"canonical row {name!r} not found in {SHEET}")
            viable = []
            for index in candidates:
                yaw = float(cell(sheet, index, "yaw_deg"))
                if abs(yaw) > MAX_ABS_YAW_DEG:
                    continue
                repetition = float(cell(sheet, index, "repetition"))
                viable.append((repetition, index, yaw))
            if not viable:
                raise SystemExit(
                    f"no |YAW| <= {MAX_ABS_YAW_DEG} row for {name!r} in {SHEET}"
                )
            repetition, index, yaw = min(viable)
            selected.append(
                {
                    "file_name": name,
                    "sequence": sequence,
                    "wind_speed_m_s": speed,
                    "repetition": int(repetition),
                    "yaw_deg": yaw,
                    "row": index,
                }
            )
    return selected


def performance_record(sheet, entry):
    rpm = float(cell(sheet, entry["row"], "rpm"))
    tsr = rpm * 2.0 * math.pi * ROTOR_RADIUS / (60.0 * entry["wind_speed_m_s"])
    pressure = float(cell(sheet, entry["row"], "wtbaro"))
    temperature = float(cell(sheet, entry["row"], "wtatemp"))
    rho = pressure / (RHO_GAS_CONSTANT * (temperature + KELVIN_OFFSET))
    return {
        "file_name": entry["file_name"],
        "sequence": entry["sequence"],
        "wind_speed_m_s": entry["wind_speed_m_s"],
        "yaw_token": cell(sheet, entry["row"], "yaw_token"),
        "yaw_deg": entry["yaw_deg"],
        "repetition": entry["repetition"],
        "rpm": rpm,
        "tsr": tsr,
        "rotpow_kw": float(cell(sheet, entry["row"], "rotpow")),
        "lsstqcor_nm": float(cell(sheet, entry["row"], "lsstqcor")),
        "eaeroth_n": float(cell(sheet, entry["row"], "eaeroth")),
        "wtbaro_pa": pressure,
        "wtatemp_degc": temperature,
        "rho_kg_m3": rho,
    }


def spanwise_records(sheet, entry):
    records = []
    for r_over_r, (cn_key, ct_key, cm_key) in sorted(SPANWISE_COLUMNS.items()):
        records.append(
            {
                "file_name": entry["file_name"],
                "sequence": entry["sequence"],
                "wind_speed_m_s": entry["wind_speed_m_s"],
                "r_over_R": r_over_r,
                "cn": float(cell(sheet, entry["row"], cn_key)),
                "ct": float(cell(sheet, entry["row"], ct_key)),
                "cm": float(cell(sheet, entry["row"], cm_key)),
            }
        )
    return records


def write_csv(path: Path, header, records) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=header)
        writer.writeheader()
        writer.writerows(records)
    print(f"wrote {path.relative_to(ROOT)} ({len(records)} rows)")


def render(workbook: Path):
    sheet, rows = load_rows(workbook)
    selected = select_rows(sheet, rows)
    files = {
        "sequence_H_performance.csv": (PERFORMANCE_HEADER, []),
        "sequence_H_spanwise.csv": (SPANWISE_HEADER, []),
        "sequence_S_performance.csv": (PERFORMANCE_HEADER, []),
        "sequence_S_spanwise.csv": (SPANWISE_HEADER, []),
    }
    for entry in selected:
        sequence = entry["sequence"]
        files[f"sequence_{sequence}_performance.csv"][1].append(
            performance_record(sheet, entry)
        )
        files[f"sequence_{sequence}_spanwise.csv"][1].extend(
            spanwise_records(sheet, entry)
        )
        print(
            f"{entry['file_name']}: {entry['wind_speed_m_s']:g} m/s, "
            f"yaw {entry['yaw_deg']:+.3f} deg, repetition "
            f"{entry['repetition']}, sheet row {entry['row'] + 1}"
        )
    return files


def render_text(header, records) -> str:
    import io

    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=header, lineterminator="\n")
    writer.writeheader()
    writer.writerows(records)
    return stream.getvalue()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workbook",
        type=Path,
        required=False,
        help="path to wt_loads_statistics.xls (17 MB, not committed)",
    )
    parser.add_argument(
        "--out-dir", type=Path, default=EXPERIMENT_DIR,
        help="output directory for the derived CSVs",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="compare against the committed CSVs instead of writing",
    )
    args = parser.parse_args(argv)
    if args.workbook is None:
        parser.error("--workbook is required")

    files = render(args.workbook)
    exit_code = 0
    for name, (header, records) in files.items():
        path = args.out_dir / name
        rendered = render_text(header, records)
        if args.check:
            if not path.exists():
                print(f"missing {path}", file=sys.stderr)
                exit_code = 1
            elif path.read_text(encoding="utf-8") != rendered:
                print(f"stale {path}", file=sys.stderr)
                exit_code = 1
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(rendered, encoding="utf-8")
            print(f"wrote {path.relative_to(ROOT)} ({len(records)} rows)")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
