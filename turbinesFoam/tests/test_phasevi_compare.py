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


def write_run(root: Path, sign: float = 1.0, tsr: float = TSR,
              steps: int = 4, elements: int = 3) -> Path:
    """Write a synthetic run directory (placeholder values, not results)."""
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
    return root


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
