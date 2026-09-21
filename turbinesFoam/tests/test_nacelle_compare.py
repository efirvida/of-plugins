# SPDX-License-Identifier: GPL-3.0-or-later
"""Pure-Python contracts for the periodic-nacelle comparison tool.

Covers the metric definitions (time-averaged <u>(z), resolved k(z), CD), the
station reconciliation (exact panels for 1R/3R/5R/7R, the bracketed 10R
stretch), the per-grid acceptance scope, the fail-loud exit codes and the
restart-file merge. The run data below is synthetic: the profiles are built
from the committed reference so the expected metrics are analytic, and no
OpenFOAM environment is needed.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).resolve().parent
TURBINESFOAM = TESTS_DIR.parent
PACKAGE = TURBINESFOAM / "validation" / "nacelle-asn"
sys.path.insert(0, str(PACKAGE / "scripts"))

import compareNacelle as compare  # noqa: E402

REFERENCE = PACKAGE / "data" / "reference"
RADIUS = 1.0
PERIOD_S = 10.0
# Window [60, 210] s from config/case.yaml: 150 one-second samples cover it.
DEFAULT_TIMES = [60.0 + index for index in range(150)]
# Explicit panel mapping of the fixture (independent of the module under test).
PANEL_MAP = {1.0: (1.0,), 3.0: (3.0,), 5.0: (5.0,), 7.0: (7.0,), 10.0: (9.0, 11.0)}


def read_panel(station: float, quantity: str) -> list[tuple[float, float]]:
    name = f"{quantity}_{station:g}R.csv"
    column = "u_over_U" if quantity == "u" else "k_over_U2"
    points = []
    with (REFERENCE / name).open(encoding="utf-8") as stream:
        reader = csv.reader(stream)
        assert next(reader) == ["z_over_D", column]
        for row in reader:
            if row:
                points.append((float(row[0]), float(row[1])))
    points.sort()
    return points


def interp(points: list[tuple[float, float]], x: float) -> float:
    if x <= points[0][0]:
        return points[0][1]
    if x >= points[-1][0]:
        return points[-1][1]
    for (x0, v0), (x1, v1) in zip(points, points[1:]):
        if x0 <= x <= x1:
            return v0 if x1 == x0 else v0 + (x - x0) / (x1 - x0) * (v1 - v0)
    return points[-1][1]


def write_run(
    root: Path,
    *,
    times: list[float] | None = None,
    fx: float = 0.754,
    perturb: dict[float, float] | None = None,
    split_at: int | None = None,
) -> Path:
    """Write a synthetic run directory: probe series + total-force CSV.

    The u profile is the mean of the mapped reference panels (plus the
    per-station `perturb` offset); the v component oscillates with amplitude
    `2 sqrt(k)` at a 10 s period, so the resolved TKE k = 1/2 <v'^2> recovers
    the reference k exactly when the window holds whole periods.
    """
    times = list(DEFAULT_TIMES if times is None else times)
    perturb = perturb or {}

    probes = []
    for station in sorted(PANEL_MAP):
        panels = PANEL_MAP[station]
        # The fixture samples the union of both figures' grids, so the metric
        # test isolates the averaging/reconstruction from interpolation error.
        z_values = sorted(
            {
                z
                for panel in panels
                for quantity in ("u", "k")
                for z, _ in read_panel(panel, quantity)
            }
        )
        for z_over_d in z_values:
            u = sum(interp(read_panel(panel, "u"), z_over_d) for panel in panels) / len(panels)
            k = sum(interp(read_panel(panel, "k"), z_over_d) for panel in panels) / len(panels)
            probes.append(
                {
                    # D = 2R: the figure axis is z/D, the probe coordinates are m.
                    "z": 2.0 * RADIUS * z_over_d,
                    "u": u + perturb.get(station, 0.0),
                    "k": k,
                    "station": station,
                }
            )

    header = [
        f"# Probe {index} ({probe['station'] + 0.05:g} 0.125 {probe['z']:.8g})"
        for index, probe in enumerate(probes)
    ]
    header.append("# Time" + "".join(f" {index}" for index in range(len(probes))))
    rows = []
    for time in times:
        fields = [f"{time:.8g}"]
        for probe in probes:
            amplitude = 2.0 * math.sqrt(probe["k"])
            v = amplitude * math.sin(2.0 * math.pi * time / PERIOD_S)
            fields.append(f"({probe['u']:.10g} {v:.10g} 0)")
        rows.append(" ".join(fields))

    if split_at is None:
        chunks = [(0, "0", rows)]
    else:
        chunks = [
            (0, "0", rows[:split_at]),
            (split_at, f"{times[split_at]:g}", rows[split_at:]),
        ]
    for _, directory, chunk in chunks:
        if not chunk:
            continue
        path = root / "postProcessing" / "profiles" / directory / "U"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(header + chunk) + "\n", encoding="utf-8")

    force = root / "postProcessing" / "nacelle" / "nacelle.csv"
    force.parent.mkdir(parents=True, exist_ok=True)
    with force.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(["time", "fx", "fy", "fz", "f_n_mag", "f_tau_mag", "cd"])
        for time in times:
            writer.writerow([
                f"{time:.8g}", f"{fx:.8g}", "0", "0", f"{fx:.8g}", "0",
                f"{abs(fx) / (0.5 * math.pi):.8g}",
            ])
    return root


def read_metrics(out: Path) -> dict:
    return json.loads((out / "metrics.json").read_text(encoding="utf-8"))


def test_reference_panels_and_station_mapping():
    cfg = compare.load_config(compare.DEFAULT_CONFIG)
    panels = compare.load_reference(compare.DEFAULT_REFERENCE)
    # Ten stations, each with the two figures (u and k): the paper has no 10R
    # panel, so no 10R reference file exists.
    assert len(panels) == 10
    assert 10.0 not in panels
    for station in (1.0, 9.0, 11.0, 19.0):
        assert panels[station]["u"].is_file()
        assert panels[station]["k"].is_file()

    assert compare.metrics_stations_R(cfg) == [1.0, 3.0, 5.0, 7.0, 10.0]
    mapping = compare.station_mapping(cfg, panels)
    for station in (1.0, 3.0, 5.0, 7.0):
        assert mapping[station]["mapping"] == "exact"
        assert mapping[station]["panels_R"] == [station]
    # The 10R stretch is bracketed by 9R and 11R: no 10R profile is invented.
    assert mapping[10.0]["mapping"] == "bracket"
    assert mapping[10.0]["panels_R"] == [9.0, 11.0]
    assert mapping[10.0]["offsets_R"] == [-1.0, 1.0]


def test_profile_interpolation_clamps_at_the_ends():
    xs = [-1.0, 0.0, 1.0]
    values = [0.0, 0.5, 1.0]
    assert compare.interpolate(xs, values, 0.5) == pytest.approx(0.75)
    # A run profile may not reach the outer reference points; the ends clamp
    # (documented) instead of extrapolating.
    assert compare.interpolate(xs, values, -2.0) == 0.0
    assert compare.interpolate(xs, values, 2.0) == 1.0


def test_clean_run_reproduces_reference_and_cd(tmp_path):
    run = write_run(tmp_path / "run-coarse")
    out = tmp_path / "out"
    assert compare.main(["--run-dir", str(run), "--out", str(out)]) == compare.EXIT_OK
    metrics = read_metrics(out)

    for key, entry in metrics["stations"].items():
        if entry["mapping"] == "exact":
            # A run built from a panel reproduces it to float noise.
            assert entry["worst"]["rms_u"] < 1.0e-6, key
            assert entry["worst"]["rms_k"] < 1.0e-6, key
        else:
            # 10R is bracketed by 9R and 11R: the best achievable agreement is
            # half their difference (~0.008 u/U, ~6e-4 k/U^2).
            assert entry["worst"]["rms_u"] < 0.02, key
            assert entry["worst"]["rms_k"] < 0.002, key
        assert entry["worst"]["within_band"] is True, key
    assert metrics["acceptance"]["claimed_R"] == [1.0, 10.0]
    assert metrics["acceptance"]["verdict"] == "PASS"
    assert metrics["window"]["mode"] == "full"
    assert "resolved" in metrics["definitions"]["k"]

    written_fx = float(f"{0.754:.8g}")  # the fixture writes 8 significant digits
    drag = metrics["drag"]
    assert drag["cd"] == pytest.approx(written_fx / (0.5 * math.pi), rel=1.0e-12)
    assert drag["cd_from_csv"] == pytest.approx(drag["cd"], rel=1.0e-5)
    assert metrics["datum"] == {"cd": pytest.approx(0.48), "is_target": False}

    rows = list(csv.DictReader((out / "profiles.csv").open(encoding="utf-8")))
    assert rows
    assert {row["quantity"] for row in rows} == {"u", "k"}


def test_cd_is_a_datum_not_an_acceptance_target(tmp_path):
    fx = 0.30 * 0.5 * math.pi  # CD = 0.30, far from the paper's 0.48 datum
    run = write_run(tmp_path / "run-coarse", fx=fx)
    out = tmp_path / "out"
    assert compare.main(["--run-dir", str(run), "--out", str(out)]) == compare.EXIT_OK
    metrics = read_metrics(out)
    written_fx = float(f"{fx:.8g}")
    assert metrics["drag"]["cd"] == pytest.approx(
        written_fx / (0.5 * math.pi), rel=1.0e-12
    )
    assert metrics["datum"]["is_target"] is False
    assert metrics["acceptance"]["verdict"] == "PASS"


def test_coarse_scope_ignores_the_mid_wake(tmp_path):
    # 5R is not claimed by the coarse grid, so a perturbed 5R must not fail.
    run = write_run(tmp_path / "run-coarse", perturb={5.0: 0.3})
    out = tmp_path / "out"
    assert compare.main(["--run-dir", str(run), "--out", str(out)]) == compare.EXIT_OK
    metrics = read_metrics(out)
    assert metrics["stations"]["5R"]["claimed"] is False
    assert metrics["stations"]["5R"]["worst"]["within_band"] is False

    # 1R is claimed: the same perturbation must fail the comparison.
    run_b = write_run(tmp_path / "run-coarse-b", perturb={1.0: 0.3})
    assert compare.main([
        "--run-dir", str(run_b), "--out", str(tmp_path / "out-b"),
    ]) == compare.EXIT_ACCEPTANCE


def test_medium_scope_claims_every_station(tmp_path):
    run = write_run(tmp_path / "run-medium")
    out = tmp_path / "out"
    assert compare.main(["--run-dir", str(run), "--out", str(out)]) == compare.EXIT_OK
    metrics = read_metrics(out)
    assert metrics["acceptance"]["claimed_R"] == [1.0, 3.0, 5.0, 7.0, 10.0]
    assert metrics["acceptance"]["source"].startswith("acceptance.medium")
    assert metrics["case"]["mesh"] == "medium"

    run_b = write_run(tmp_path / "run-medium-b", perturb={5.0: 0.3})
    assert compare.main([
        "--run-dir", str(run_b), "--out", str(tmp_path / "out-b"),
    ]) == compare.EXIT_ACCEPTANCE


def test_far_wake_stretch_uses_both_bracketing_panels(tmp_path):
    run = write_run(tmp_path / "run-coarse", perturb={10.0: 0.3})
    out = tmp_path / "out"
    assert compare.main(["--run-dir", str(run), "--out", str(out)]) == compare.EXIT_ACCEPTANCE
    metrics = read_metrics(out)

    reconciliation = metrics["station_reconciliation"]
    assert 10.0 not in reconciliation["paper_panels_R"]
    assert "no 10R profile is invented" in reconciliation["note"]
    assert reconciliation["mapping"]["10R"]["panels_R"] == [9.0, 11.0]
    panels = metrics["stations"]["10R"]["panels"]
    assert {entry["station_panel_R"] for entry in panels} == {9.0, 11.0}
    assert all(entry["within_band"] is False for entry in panels)

    rows = [
        row for row in csv.DictReader((out / "profiles.csv").open(encoding="utf-8"))
        if row["station_R"] == "10"
    ]
    assert rows
    assert {row["station_panel_R"] for row in rows} == {"9", "11"}


def test_restart_files_are_merged_by_time(tmp_path):
    run = write_run(tmp_path / "run-coarse", split_at=60)
    out = tmp_path / "out"
    assert compare.main(["--run-dir", str(run), "--out", str(out)]) == compare.EXIT_OK
    metrics = read_metrics(out)
    assert metrics["case"]["samples"] == len(DEFAULT_TIMES)
    assert len(metrics["case"]["probe_files"]) == 2
    assert metrics["stations"]["1R"]["worst"]["rms_u"] < 1.0e-6


def test_short_window_exit_codes(tmp_path):
    times = [100.0 + index for index in range(11)]
    run = write_run(tmp_path / "run-coarse", times=times)
    out = tmp_path / "out"
    assert compare.main([
        "--run-dir", str(run), "--out", str(out),
    ]) == compare.EXIT_SHORT_WINDOW

    out_allowed = tmp_path / "out-allowed"
    assert compare.main([
        "--run-dir", str(run), "--out", str(out_allowed), "--allow-short-window",
    ]) == compare.EXIT_OK
    assert read_metrics(out_allowed)["window"]["mode"] == "short"


def test_fail_loud_exit_codes(tmp_path):
    missing = compare.main([
        "--run-dir", str(tmp_path / "missing"), "--out", str(tmp_path / "o1"),
    ])
    assert missing == compare.EXIT_MISSING_INPUT

    run = write_run(tmp_path / "run-coarse")
    (run / "postProcessing" / "profiles" / "0" / "U").unlink()
    assert compare.main([
        "--run-dir", str(run), "--out", str(tmp_path / "o2"),
    ]) == compare.EXIT_MISSING_INPUT

    run_b = write_run(tmp_path / "run-coarse-b")
    (run_b / "postProcessing" / "nacelle" / "nacelle.csv").unlink()
    assert compare.main([
        "--run-dir", str(run_b), "--out", str(tmp_path / "o3"),
    ]) == compare.EXIT_MISSING_INPUT


def test_band_override_can_fail_a_claimed_station(tmp_path):
    run = write_run(tmp_path / "run-coarse", perturb={1.0: 0.05})
    assert compare.main([
        "--run-dir", str(run), "--out", str(tmp_path / "out"),
    ]) == compare.EXIT_OK
    assert compare.main([
        "--run-dir", str(run), "--out", str(tmp_path / "out-tight"),
        "--band-u", "0.01",
    ]) == compare.EXIT_ACCEPTANCE
