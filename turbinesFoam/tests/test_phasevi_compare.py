# SPDX-License-Identifier: GPL-3.0-or-later
"""Pure-Python contracts for the Phase VI comparison tool.

Covers the Tooling layer of the design testing strategy: the F1 metric
definitions, the `root_dist` -> r/R mapping and station interpolation, the
drift flag, the fail-loud exit codes and the 7 m/s sign gate. The run data
below is synthetic placeholders generated under pytest's `tmp_path`; no
fixture file is committed and no result is a validation result. These tests
never require OpenFOAM.
"""

from __future__ import annotations

import csv
import json
import math
import re
import sys
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).resolve().parent
PACKAGE = TESTS_DIR.parent / "validation" / "phaseVI"
sys.path.insert(0, str(PACKAGE / "scripts"))

import comparePhaseVI as compare  # noqa: E402

SPEED = 7.0
TSR = 5.408
TURBINE_COLUMNS = [
    "time", "angle_deg", "tsr", "cp", "cd", "ct", "cd_blade1", "ct_blade1",
]
ELEMENT_COLUMNS = [
    "time", "root_dist", "x", "y", "z", "rel_vel_mag", "Re", "alpha_deg",
    "alpha_geom_deg", "cl", "cd", "fx", "fy", "fz", "end_effect_factor",
    "c_ref_t", "c_ref_n", "f_ref_t", "f_ref_n",
]
#: The real `bladeSurface` CSV schemas written by
#: `src/fvOptions/bladeSurface/bladeSurfaceSource.C` (station table at
#: `:61-63`, opt-in node table at `:70-72`). The three-way fixture below must
#: use these exact columns; `test_surface_csv_schema_matches_writer` pins the
#: constants against the writer source so a schema drift cannot silently
#: invalidate the comparison contract.
SURFACE_STATION_COLUMNS = [
    "time", "station", "root_dist", "area", "force_x", "force_y", "force_z",
    "c_ref_n", "c_ref_t", "f_ref_n", "f_ref_t",
]
SURFACE_NODE_COLUMNS = [
    "time", "node", "x", "y", "z", "nx", "ny", "nz", "fx", "fy", "fz",
    "area", "station", "chord_fraction",
]
SURFACE_STATION_NAME = "turbine.blade1.surface.csv"
SURFACE_NODES_NAME = "turbine.blade1.surface_nodes.csv"


def write_run(root: Path, sign: float = 1.0, tsr: float = TSR,
              steps: int = 4, elements: int = 3,
              surface_stations: int = 0) -> Path:
    """Write a synthetic run directory (placeholder values, not results).

    `surface_stations > 0` additionally writes the ASM-mesh surface output
    (`postProcessing/bladeSurface/<name>.csv` with the real station schema,
    plus the opt-in `*_nodes.csv` table the comparison must ignore).
    """
    turbine = root / compare.TURBINE_CSV
    turbine.parent.mkdir(parents=True, exist_ok=True)
    with turbine.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(TURBINE_COLUMNS)
        cp = 0.30
        ct = cp / tsr
        for step in range(steps):
            writer.writerow([
                f"{step * 0.008:.6f}", f"{step * 3.0:.4f}", f"{tsr:.6f}",
                f"{cp:.6f}", f"{sign * 0.55:.6f}", f"{sign * ct:.6f}",
                f"{sign * 0.27:.6f}", f"{sign * ct:.6f}",
            ])
    element_dir = root / compare.ELEMENT_DIR
    element_dir.mkdir(parents=True, exist_ok=True)
    for index in range(elements):
        root_dist = index / (elements - 1)
        with (element_dir / f"blade1Element{index}.csv").open(
            "w", newline="", encoding="utf-8"
        ) as stream:
            writer = csv.writer(stream)
            writer.writerow(ELEMENT_COLUMNS)
            for step in range(steps):
                writer.writerow([
                    f"{step * 0.008:.6f}", f"{root_dist:.6f}", "0", "0",
                    "12.192", "70", "1e6", "6", "6", "0.9", "0.01",
                    f"{sign * 200.0:.4f}", "0", "0", "1",
                    f"{sign * 0.08:.6f}", f"{sign * 0.60:.6f}", "1.5", "12",
                ])
    if surface_stations:
        surface_dir = root / compare.SURFACE_DIR
        surface_dir.mkdir(parents=True, exist_ok=True)
        with (surface_dir / SURFACE_STATION_NAME).open(
            "w", newline="", encoding="utf-8"
        ) as stream:
            writer = csv.writer(stream)
            writer.writerow(SURFACE_STATION_COLUMNS)
            for step in range(steps):
                for index in range(surface_stations):
                    root_dist = index / (surface_stations - 1)
                    writer.writerow([
                        f"{step * 0.008:.6f}", index, f"{root_dist:.6f}",
                        "0.050000",
                        f"{sign * 80.0:.4f}", "0", "0",
                        f"{sign * (0.60 - 0.40 * root_dist):.6f}",
                        f"{sign * (0.08 - 0.06 * root_dist):.6f}",
                        "1.5", "12",
                    ])
        # The opt-in per-node table is not a station table: its basename is the
        # station name plus `_nodes.csv`, and `read_surface_stations` must
        # exclude it instead of trying to read `root_dist`/`c_ref_n` from it.
        with (surface_dir / SURFACE_NODES_NAME).open(
            "w", newline="", encoding="utf-8"
        ) as stream:
            writer = csv.writer(stream)
            writer.writerow(SURFACE_NODE_COLUMNS)
            writer.writerow([
                "0.000000", "0", "0.1", "0", "12.192", "0", "0", "1",
                f"{sign * 8.0:.4f}", "0", "0", "0.050000", "0.5", "0.25",
            ])
    return root


def linear_at(points, x):
    """Independent linear interpolation oracle for the surface conversion."""
    points = sorted(points)
    if x <= points[0][0]:
        return points[0][1]
    if x >= points[-1][0]:
        return points[-1][1]
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 <= x <= x1:
            ratio = 0.0 if x1 == x0 else (x - x0) / (x1 - x0)
            return y0 + ratio * (y1 - y0)
    raise AssertionError(f"{x} outside the profile")


def fake_model(tsr: float = TSR, cp: float = 0.30, cd: float = 0.55,
               c_ref_n: float = 0.60, c_ref_t: float = 0.08):
    return {
        "stats": {"tsr": tsr},
        "metrics": {"tsr_config": TSR, "cp": cp, "cd": cd},
        "spanwise": {
            f"{station:.2f}": (c_ref_n, c_ref_t) for station in compare.STATIONS
        },
    }


def test_metric_definitions_match_f1():
    cfg = compare.load_config(compare.DEFAULT_CONFIG)
    radius = float(cfg["turbine"]["radius"])
    area = math.pi * radius**2
    rho = 1.2
    stats = {
        "cp": 0.30, "cd": 0.55, "ct": 0.06,
        "blade_cd": {"blade1": 0.27, "blade2": 0.27},
    }
    metrics = compare.turbine_metrics(stats, rho, SPEED, cfg, "blade")
    q_dyn = 0.5 * rho * area * SPEED**2
    # F1: Q = ct * q_dyn * R with R = rotorRadius
    assert metrics["torque_nm"] == pytest.approx(stats["ct"] * q_dyn * radius)
    # F1: P = cp * q_dyn * Uinf with cp = ct * TSR
    assert metrics["power_w"] == pytest.approx(stats["cp"] * q_dyn * SPEED)
    # blade-only thrust headline; rotor (hub included) reported secondary
    assert metrics["thrust_n"] == pytest.approx(0.54 * q_dyn)
    assert metrics["thrust_scope"] == "blade"
    assert metrics["thrust_rotor_n"] == pytest.approx(stats["cd"] * q_dyn)
    rotor = compare.turbine_metrics(stats, rho, SPEED, cfg, "rotor")
    assert rotor["thrust_scope"] == "rotor"
    assert rotor["thrust_n"] == pytest.approx(stats["cd"] * q_dyn)


def test_spanwise_mapping_and_interpolation():
    cfg = compare.load_config(compare.DEFAULT_CONFIG)
    radius = float(cfg["turbine"]["radius"])
    rows_root = [
        {"time": "0.5", "root_dist": "0.0", "c_ref_n": "0.6", "c_ref_t": "0.08"},
    ]
    rows_tip = [
        {"time": "0.5", "root_dist": "1.0", "c_ref_n": "0.2", "c_ref_t": "0.02"},
    ]
    profile = compare.spanwise_profile(
        [(Path("root.csv"), rows_root), (Path("tip.csv"), rows_tip)],
        start=0.0, end=1.0, allow_short=True, radius=radius,
    )
    # root_dist is blade-normalized: 0 maps to the root cutout, 1 to the tip
    assert profile[0][0] == pytest.approx(compare.ROOT_CUTOUT_RADIUS / radius)
    assert profile[-1][0] == pytest.approx(1.0)

    simple = [(0.0, 0.6, 0.08), (1.0, 0.2, 0.02)]
    assert compare.interpolate(simple, 0.5) == pytest.approx((0.4, 0.05))
    assert compare.interpolate(simple, -1.0) == pytest.approx((0.6, 0.08))
    assert compare.interpolate(simple, 2.0) == pytest.approx((0.2, 0.02))


def test_drift_flag_threshold():
    period = 1.0
    stable = [
        {"time": f"{4.0 + index * 0.1:.2f}", "cp": "0.30", "cd": "0.55"}
        for index in range(30)
    ]
    assert compare.drift_flags(stable, start=4.0, period=period)["flag"] is False
    drifting = [
        {"time": f"{4.0 + index * 0.1:.2f}",
         "cp": "0.30" if index < 10 else ("0.32" if index < 20 else "0.30"),
         "cd": "0.55"}
        for index in range(30)
    ]
    flags = compare.drift_flags(drifting, start=4.0, period=period)
    assert flags["cp"]["variation"] > compare.DRIFT_TOLERANCE
    assert flags["flag"] is True


def test_sign_gate_pass_and_fail():
    passing = compare.sign_gate({"alm": fake_model(), "asm": fake_model()})
    assert passing["pass"] is True
    assert all(passing["alm"]["checks"].values())
    mirrored = compare.sign_gate({"alm": fake_model(c_ref_n=-0.6)})
    assert mirrored["pass"] is False
    assert mirrored["alm"]["checks"]["c_ref_n_positive"] is False
    no_power = compare.sign_gate({"alm": fake_model(cp=-0.3)})
    assert no_power["pass"] is False


def test_main_dry_merge_and_definitions(tmp_path):
    run = write_run(tmp_path / "run-alm")
    out = tmp_path / "out"
    rc = compare.main([
        "--run-dir", str(run), "--out", str(out),
        "--allow-short-window", "--revolutions", "0", "1",
    ])
    assert rc == compare.EXIT_OK
    metrics = json.loads((out / "metrics.json").read_text(encoding="utf-8"))
    assert "R = rotorRadius" in metrics["definitions"]["torque"]
    assert metrics["definitions"]["power"] == (
        "P = cp * q_dyn * Uinf with cp = ct * TSR, A = pi R^2"
    )
    assert metrics["case"]["rho_source"].startswith("WTBARO/WTATEMP")
    assert metrics["cm_excluded"]["follow_up"].startswith("add momentCoefficient_")
    assert any("sub-cell" in item for item in metrics["limitations"])
    rows = list(csv.DictReader((out / "spanwise_comparison.csv").open()))
    assert [row["r_over_R"] for row in rows] == [
        f"{station:.2f}" for station in compare.STATIONS
    ]
    assert (out / "turbine_comparison.csv").is_file()
    assert (out / "report.txt").is_file()


def test_main_fail_loud_exit_codes(tmp_path):
    missing = compare.main([
        "--alm-dir", str(tmp_path / "missing"), "--out", str(tmp_path / "out1"),
    ])
    assert missing == compare.EXIT_MISSING_INPUT

    run = write_run(tmp_path / "run-tsr", tsr=4.0)
    mismatch = compare.main([
        "--run-dir", str(run), "--out", str(tmp_path / "out2"),
        "--allow-short-window",
    ])
    assert mismatch == compare.EXIT_TSR_MISMATCH

    short = compare.main([
        "--run-dir", str(run), "--out", str(tmp_path / "out3"),
    ])
    assert short == compare.EXIT_SHORT_WINDOW

    empty = tmp_path / "empty-experiment"
    empty.mkdir()
    absent = compare.main([
        "--run-dir", str(run), "--experiment", str(empty),
        "--out", str(tmp_path / "out4"), "--allow-short-window",
    ])
    assert absent == compare.EXIT_MISSING_INPUT


def quoted_header(source: str, marker: str) -> list[str]:
    """Rebuild a C++ `<< "a,b" << "c,d" << endl` header as a column list."""
    match = re.search(rf"\*{marker}_\s*<<(.*?)<<\s*endl", source, re.DOTALL)
    assert match is not None, f"no {marker} header in the writer"
    chunks = re.findall(r'"([^"]*)"', match.group(1))
    return [
        column for chunk in chunks for column in chunk.split(",") if column
    ]


def test_surface_csv_schema_matches_writer():
    """The fixture schema is the one `bladeSurfaceSource.C` actually writes.

    The comparison consumes a C++-produced CSV; if the writer schema drifts,
    the synthetic fixture would keep passing while real runs break. Pin both
    headers against the writer source.
    """
    source = (
        TESTS_DIR.parent
        / "src" / "fvOptions" / "bladeSurface" / "bladeSurfaceSource.C"
    ).read_text(encoding="utf-8")
    assert quoted_header(source, "stationFile") == SURFACE_STATION_COLUMNS
    assert quoted_header(source, "nodeFile") == SURFACE_NODE_COLUMNS


def test_main_three_way_surface_merge(tmp_path):
    """Three models in both tables; surface rows use the element r/R mapping."""
    alm = write_run(tmp_path / "alm")
    asm = write_run(tmp_path / "asm")
    asm_mesh = write_run(tmp_path / "asm-mesh", surface_stations=5)
    out = tmp_path / "out"
    rc = compare.main([
        "--alm-dir", str(alm),
        "--asm-dir", str(asm),
        "--asm-mesh-dir", str(asm_mesh),
        "--out", str(out),
        "--allow-short-window", "--revolutions", "0", "1",
        "--sign-gate",
    ])
    assert rc == compare.EXIT_OK

    turbine_rows = list(
        csv.DictReader((out / "turbine_comparison.csv").open(encoding="utf-8"))
    )
    spanwise_rows = list(
        csv.DictReader((out / "spanwise_comparison.csv").open(encoding="utf-8"))
    )
    models = ["alm", "asm", "asm-mesh"]
    assert sorted({row["model"] for row in turbine_rows}) == models
    assert sorted({row["model"] for row in spanwise_rows}) == models
    assert len(turbine_rows) == 3 * len(compare.TURBINE_BANDS)
    for model in models:
        stations = [
            row["r_over_R"] for row in spanwise_rows if row["model"] == model
        ]
        assert stations == [f"{station:.2f}" for station in compare.STATIONS]
        assert all(
            row["within_band"] in ("yes", "no")
            for row in spanwise_rows
            if row["model"] == model
        )

    metrics = json.loads((out / "metrics.json").read_text(encoding="utf-8"))
    assert "surface_conversion" in metrics["definitions"]
    surface = metrics["models"]["asm-mesh"]["metrics"]
    assert surface["surface_stations"] == 5
    # The synthetic directory has no run.json/twin, so the audit fields are
    # best-effort nulls; a real run carries the staged hash (test_blade_stage).
    assert surface["staged_stl_sha256"] is None
    assert surface["surface_kernel"] is None

    # Independent conversion: root_dist -> r/R with the element formula, then
    # interpolation to the five measured stations.
    radius = float(compare.load_config(compare.DEFAULT_CONFIG)["turbine"]["radius"])
    span = 1.0 - compare.ROOT_CUTOUT_RADIUS / radius
    stations = compare.read_surface_stations(asm_mesh, 0.0, 1.0, True)
    assert [round(root_dist, 6) for root_dist, _, _ in stations] == [
        0.0, 0.25, 0.5, 0.75, 1.0,
    ]
    mapped_n = [
        (compare.ROOT_CUTOUT_RADIUS / radius + root_dist * span, c_ref_n)
        for root_dist, c_ref_n, _ in stations
    ]
    mapped_t = [
        (compare.ROOT_CUTOUT_RADIUS / radius + root_dist * span, c_ref_t)
        for root_dist, _, c_ref_t in stations
    ]
    for station in compare.STATIONS:
        values = metrics["models"]["asm-mesh"]["spanwise"][f"{station:.2f}"]
        assert values[0] == pytest.approx(linear_at(mapped_n, station), rel=1e-12)
        assert values[1] == pytest.approx(linear_at(mapped_t, station), rel=1e-12)

    gate = json.loads((out / "sign_gate.json").read_text(encoding="utf-8"))
    assert gate["pass"] is True
    assert sorted(key for key in gate if key != "pass") == models
    assert all(gate[model]["pass"] for model in models)


def test_main_three_way_missing_surface_input(tmp_path):
    """A missing or incomplete ASM-mesh directory exits non-zero."""
    alm = write_run(tmp_path / "alm")
    flags = ["--allow-short-window", "--revolutions", "0", "1"]

    missing = compare.main([
        "--alm-dir", str(alm),
        "--asm-mesh-dir", str(tmp_path / "missing-surface"),
        "--out", str(tmp_path / "out-missing"),
    ] + flags)
    assert missing == compare.EXIT_MISSING_INPUT

    # A directory with a turbine CSV but no `postProcessing/bladeSurface/`
    # station table is incomplete, not a silent ALM-style fallback.
    incomplete = write_run(tmp_path / "asm-mesh-no-surface")
    absent = compare.main([
        "--alm-dir", str(alm),
        "--asm-mesh-dir", str(incomplete),
        "--out", str(tmp_path / "out-incomplete"),
    ] + flags)
    assert absent == compare.EXIT_MISSING_INPUT
