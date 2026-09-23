#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Compare NREL Phase VI simulation runs with the measured WDH data.

Reads the turbine-level CSV and the element-level CSVs for the ALM and ASM
models and, for the mesh-backed ASM variant (`asm-mesh`, `--asm-mesh-dir`), the
turbine-level CSV plus the per-station surface CSVs
(`postProcessing/bladeSurface/*.csv`), merges them with the committed
experimental CSVs under `data/experiment/`, and writes turbine-level and
spanwise simulation-versus-experiment tables plus `metrics.json` and
`report.txt`.

The surface station rows already carry the element public `c_ref_n`/`c_ref_t`
definitions, so the surface conversion only applies the element `root_dist` ->
r/R mapping; no coefficient definition is redefined.

Metric definitions (fixed against the NREL report before any result is
claimed; F1):

    q_dyn   = 1/2 * rho * A * Uinf^2              [N]
    Q       = ct * q_dyn * R                      [Nm]  shaft torque
    P       = cp * q_dyn * Uinf                   [W]   rotor power
    cp      = ct * TSR                            [-]
    T_blade = sum_b cd_blade * q_dyn              [N]   blade-only thrust
    T_rotor = cd * q_dyn                          [N]   rotor thrust (hub
                                                        included; labelled
                                                        secondary)
    rho     = WTBARO / (287.058 * (WTATEMP + 273.15))  [kg/m^3]

with `R = rotorRadius` and `A = pi * R^2` (the turbine.csv `ct` is the torque
coefficient and `cd` the axial force coefficient). Spanwise `c_ref_n`/`c_ref_t`
are compared against the measured CN/CT at 30/47/63/80/95 % span. Spanwise CM
is excluded because the element CSV has no `cm` column (documented follow-up).

Exit codes:
    0  comparison produced (band/drift flags are reported, not hidden)
    1  missing input or failed sign gate
    2  averaging window shorter than four revolutions without
       --allow-short-window
    3  configured TSR does not match the run
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
ROOT = SCRIPTS.parent
sys.path.insert(0, str(ROOT / "tools"))

from case_config import load_config, select_speed  # noqa: E402

DEFAULT_EXPERIMENT = ROOT / "data" / "experiment"
DEFAULT_CONFIG = ROOT / "config" / "case.yaml"

TURBINE_CSV = Path("postProcessing") / "turbines" / "0" / "turbine.csv"
ELEMENT_DIR = Path("postProcessing") / "actuatorLineElements" / "0"
SURFACE_DIR = Path("postProcessing") / "bladeSurface"

STATIONS = (0.30, 0.47, 0.63, 0.80, 0.95)
TURBINE_BANDS = {"power": 0.15, "torque": 0.15, "thrust": 0.15}
SPANWISE_BAND_MIN = 0.15
SPANWISE_BAND_FRACTION = 0.20
TSR_TOLERANCE = 1.0e-3
DRIFT_TOLERANCE = 0.01
ROOT_CUTOUT_RADIUS = 0.5083
MAX_ABS_YAW_DEG = 0.5
RHO_GAS_CONSTANT = 287.058
KELVIN_OFFSET = 273.15

EXIT_OK = 0
EXIT_MISSING_INPUT = 1
EXIT_SHORT_WINDOW = 2
EXIT_TSR_MISMATCH = 3

LIMITATIONS = [
    "At every affordable mesh (D/32..D/64) the ASM chordwise strips are "
    "sub-cell, so the strip-wise force distribution is sub-grid: the "
    "comparison tests chord-averaged inflow and the load shift versus ALM, "
    "not a chord-resolved surface.",
    "The imported blade surface (asm-mesh) is sub-grid at the affordable "
    "meshes: background cells are 0.314/0.210/0.157 m (D/32..D/64) against a "
    "0.218-0.744 m chord, so the three-way comparison tests the model form "
    "(how the load is distributed), not resolved chordwise physics.",
    "The surface kernel/width is a confound with the imported geometry: the "
    "paper-cosine surface and the no-mesh Gaussian differ in kernel and width "
    "as well as in geometry, so the kernel-matched ASM-mesh configuration "
    "(`generate_case.py --surface-kernel gaussian`) is an ablation, not a "
    "model.",
    "The ASM-mesh variant is distribution-only: the element BEM chain, "
    "chord-averaged inflow, dynamic stall and end effects are unchanged, and "
    "the surface does not sample per-node inflow or recompute loads.",
    "URANS k-omega SST cannot capture deep-stall unsteadiness or hysteresis; "
    "the separated high-speed points (13-25 m/s) are trend and stall-onset "
    "evidence only.",
    "The tower is not modelled (rotor + hub only); EAEROTH is a blade "
    "pressure integration and the measured loads include tower shadow.",
    "Spanwise CM is excluded: the element CSV has no `cm` column (follow-up "
    "change adds `momentCoefficient_` to the element writer).",
    "Blockage is not simulated (uniform free box versus the NASA-Ames test "
    "section); documented, not corrected.",
]


def _exit(message: str, code: int) -> int:
    print(f"error: {message}", file=sys.stderr)
    return code


def read_rows(path: Path) -> list[dict[str, float]]:
    with path.open(encoding="utf-8") as stream:
        rows = []
        for row in csv.DictReader(stream):
            parsed = {}
            for key, value in row.items():
                if value is None or value == "":
                    continue
                try:
                    parsed[key] = float(value)
                except ValueError:
                    parsed[key] = value
            rows.append(parsed)
    return rows


def load_experiment(experiment_dir: Path, sequence: str, speed: float):
    performance = experiment_dir / f"sequence_{sequence}_performance.csv"
    spanwise = experiment_dir / f"sequence_{sequence}_spanwise.csv"
    if not performance.exists():
        raise FileNotFoundError(str(performance))
    if not spanwise.exists():
        raise FileNotFoundError(str(spanwise))
    row = None
    for candidate in read_rows(performance):
        if math.isclose(float(candidate["wind_speed_m_s"]), float(speed), abs_tol=1e-9):
            row = candidate
            break
    if row is None:
        raise KeyError(f"{performance}: no row for {speed:g} m/s")
    if abs(float(row["yaw_deg"])) > MAX_ABS_YAW_DEG:
        raise ValueError(f"{performance}: row yaw exceeds {MAX_ABS_YAW_DEG} deg")
    stations = {}
    for candidate in read_rows(spanwise):
        if math.isclose(float(candidate["wind_speed_m_s"]), float(speed), abs_tol=1e-9):
            stations[float(candidate["r_over_R"])] = {
                "cn": float(candidate["cn"]),
                "ct": float(candidate["ct"]),
                "cm": float(candidate["cm"]),
            }
    if not stations:
        raise KeyError(f"{spanwise}: no rows for {speed:g} m/s")
    return row, stations


def window_bounds(cfg, speed: float, sequence: str = "H") -> tuple[float, float, float]:
    """(discard time, end time, revolution period) in seconds."""
    entry = select_speed(cfg, speed, sequence)
    omega = float(entry["tsr"]) * float(speed) / float(cfg["turbine"]["radius"])
    period = 2.0 * math.pi / omega
    solver = cfg["solver"]
    return (
        float(solver["discard_revolutions"]) * period,
        float(solver["end_revolutions"]) * period,
        period,
    )


def mean(values):
    values = [value for value in values if value is not None]
    if not values:
        return None
    return sum(values) / len(values)


def analyse_turbine(rows, start, end, period, allow_short, tsr_config):
    if not rows:
        raise ValueError("turbine CSV has no rows")
    times = [float(row["time"]) for row in rows]
    window = [row for row in rows if start <= float(row["time"]) <= end]
    available = (max(times) - start) / period
    mode = "full"
    if not window:
        if not allow_short:
            raise ShortWindow(
                f"averaging window starts at {start:.4g} s (rev "
                f"{start / period:g}) but the run ends at {max(times):.4g} s"
            )
        window = rows
        mode = "short"
    elif available < 4.0:
        if not allow_short:
            raise ShortWindow(
                f"only {available:.2f} revolutions after the discarded "
                "revolutions (four required)"
            )
        mode = "short"
    revolutions = (max(float(row["time"]) for row in window) - start) / period
    tsr = mean([float(row["tsr"]) for row in window])
    if tsr is None or abs(tsr - tsr_config) > TSR_TOLERANCE:
        raise TsrMismatch(
            f"run TSR {tsr if tsr is not None else float('nan'):.6f} differs "
            f"from the configured {tsr_config:.6f} by more than "
            f"{TSR_TOLERANCE:g}"
        )
    cp = mean([float(row["cp"]) for row in window])
    cd = mean([float(row["cd"]) for row in window])
    ct = mean([float(row["ct"]) for row in window])
    blade_names = sorted(
        key[3:] for key in rows[0] if key.startswith("cd_") and key[3:]
    )
    blade_cd = {
        name: mean([float(row[f"cd_{name}"]) for row in window])
        for name in blade_names
    }
    return {
        "rows": len(window),
        "tsr": tsr,
        "cp": cp,
        "cd": cd,
        "ct": ct,
        "cp_from_ct": ct * tsr,
        "blade_cd": blade_cd,
        "revolutions_averaged": revolutions,
        "window_mode": mode,
        "window_start_s": start,
        "window_end_s": end,
        "time_max_s": max(times),
    }


def drift_flags(rows, start, period):
    """Per-revolution means of cp/cd across the window and their variation."""
    result = {}
    for key in ("cp", "cd"):
        per_rev = {}
        for row in rows:
            time = float(row["time"])
            if time < start:
                continue
            index = int((time - start) / period)
            per_rev.setdefault(index, []).append(float(row[key]))
        means = [mean(values) for index, values in sorted(per_rev.items())]
        means = [value for value in means if value is not None]
        if len(means) < 2:
            result[key] = {"variation": None, "flag": True,
                           "reason": "fewer than two revolutions in the window"}
            continue
        baseline = abs(mean(means))
        variation = (max(means) - min(means)) / baseline if baseline else None
        result[key] = {
            "per_revolution": means,
            "variation": variation,
            "flag": variation is not None and variation > DRIFT_TOLERANCE,
        }
    result["flag"] = any(result[key]["flag"] for key in ("cp", "cd"))
    return result


def read_elements(element_dir: Path):
    files = sorted(
        path for path in element_dir.glob("*.csv") if path.is_file()
    )
    if not files:
        raise FileNotFoundError(str(element_dir / "*.csv"))
    elements = []
    for path in files:
        rows = read_rows(path)
        if not rows:
            raise MissingInput(f"{path}: no data rows")
        elements.append((path, rows))
    return elements


def r_over_r_from_root_dist(root_dist, radius):
    """Map the element blade-normalized root distance to r/R.

    `root_dist` is 0 at the root cutout and 1 at the tip; the mapping is the
    one the element profile has always used (`spanwise_profile`).
    """
    span = 1.0 - ROOT_CUTOUT_RADIUS / radius
    return ROOT_CUTOUT_RADIUS / radius + root_dist * span


def spanwise_profile(elements, start, end, allow_short, radius):
    """Time-mean c_ref_n/c_ref_t per element mapped to r/R."""
    profile = []
    for path, rows in elements:
        window = [row for row in rows if start <= float(row["time"]) <= end]
        if not window:
            if not allow_short:
                raise ShortWindow(f"{path}: no samples in the averaging window")
            window = rows
        root_dist = mean([float(row["root_dist"]) for row in window])
        c_ref_n = mean([float(row["c_ref_n"]) for row in window])
        c_ref_t = mean([float(row["c_ref_t"]) for row in window])
        r_over_r = r_over_r_from_root_dist(root_dist, radius)
        profile.append((r_over_r, c_ref_n, c_ref_t))
    profile.sort()
    return profile


def read_surface_stations(run_dir: Path, start, end, allow_short):
    """Time-mean per-station c_ref_n/c_ref_t of an ASM-mesh run.

    Reads the surface station CSVs (`postProcessing/bladeSurface/*.csv`, the
    per-node `*_nodes.csv` output is not a station table), filters on the
    averaging window like `spanwise_profile` and returns
    `(root_dist, c_ref_n, c_ref_t)` tuples sorted by root distance. The rows
    already carry the element public reference coefficients, so the caller only
    applies the element r/R mapping. A missing surface directory or an empty
    station table fails loudly.
    """
    surface_dir = run_dir / SURFACE_DIR
    files = sorted(
        path for path in surface_dir.glob("*.csv")
        if path.is_file() and not path.name.endswith("_nodes.csv")
    )
    if not files:
        raise FileNotFoundError(str(surface_dir / "*.csv"))
    stations: dict[int, list[dict]] = {}
    for path in files:
        rows = read_rows(path)
        if not rows:
            raise MissingInput(f"{path}: no data rows")
        window = [row for row in rows if start <= float(row["time"]) <= end]
        if not window:
            if not allow_short:
                raise ShortWindow(f"{path}: no samples in the averaging window")
            window = rows
        for row in window:
            stations.setdefault(int(float(row["station"])), []).append(row)
    # The Phase VI twin writes station output for blade1 only (`writePerf false`
    # on blade2), but merging any station files by station id keeps exactly one
    # profile per station (identical values for symmetric blades) instead of
    # double-counting the load.
    profile = []
    for station in sorted(stations):
        rows = stations[station]
        root_dist = mean([float(row["root_dist"]) for row in rows])
        c_ref_n = mean([float(row["c_ref_n"]) for row in rows])
        c_ref_t = mean([float(row["c_ref_t"]) for row in rows])
        profile.append((root_dist, c_ref_n, c_ref_t))
    profile.sort()
    return profile


def read_surface_audit(run_dir: Path) -> dict:
    """Fair-comparison audit fields of an ASM-mesh run directory.

    `run.json` (written by `runPhaseVI.sh`) records the staged surface STL
    hash; `system/fvOptions` is the installed twin whose `kernel` key selects
    the surface kernel (`cosine` is the default and is not rendered). Both are
    best-effort: a synthetic directory without them reports null.
    """
    audit = {"staged_stl_sha256": None, "surface_kernel": None}
    run_json = run_dir / "run.json"
    if run_json.is_file():
        try:
            payload = json.loads(run_json.read_text(encoding="utf-8"))
        except ValueError:
            payload = {}
        if isinstance(payload, dict):
            audit["staged_stl_sha256"] = payload.get("staged_stl_sha256")
    fv_options = run_dir / "system" / "fvOptions"
    if fv_options.is_file():
        match = re.search(
            r"^\s*kernel\s+(\w+)\s*;",
            fv_options.read_text(encoding="utf-8"),
            re.MULTILINE,
        )
        audit["surface_kernel"] = match.group(1) if match else "cosine"
    return audit


def interpolate(profile, x):
    if not profile:
        return None
    if x <= profile[0][0]:
        return profile[0][1], profile[0][2]
    if x >= profile[-1][0]:
        return profile[-1][1], profile[-1][2]
    for (x0, n0, t0), (x1, n1, t1) in zip(profile, profile[1:]):
        if x0 <= x <= x1:
            ratio = 0.0 if x1 == x0 else (x - x0) / (x1 - x0)
            return n0 + ratio * (n1 - n0), t0 + ratio * (t1 - t0)
    return None


def turbine_metrics(stats, rho, speed, cfg, thrust_scope):
    radius = float(cfg["turbine"]["radius"])
    area = math.pi * radius**2
    q_dyn = 0.5 * rho * area * speed**2
    power = stats["cp"] * q_dyn * speed
    torque = stats["ct"] * q_dyn * radius
    thrust_rotor = stats["cd"] * q_dyn
    blade_sum = None
    if stats["blade_cd"]:
        blade_sum = sum(stats["blade_cd"].values()) * q_dyn
    scope = "blade" if (thrust_scope == "blade" and blade_sum is not None) \
        else "rotor"
    thrust = blade_sum if scope == "blade" else thrust_rotor
    return {
        "q_dyn_pa": q_dyn,
        "power_w": power,
        "power_kw": power / 1.0e3,
        "torque_nm": torque,
        "thrust_n": thrust,
        "thrust_scope": scope,
        "thrust_rotor_n": thrust_rotor,
        "thrust_blade_n": blade_sum,
        "cp": stats["cp"],
        "ct": stats["ct"],
        "cd": stats["cd"],
    }


class ShortWindow(Exception):
    pass


class TsrMismatch(Exception):
    pass


class MissingInput(Exception):
    pass


def span_limited_blade_thrust(elements, start, end, allow_short, radius):
    """Blade-only thrust from element axial forces restricted to r/R >= 0.25.

    EAEROTH integrates pressure over the panels from 25 % span to the tip and
    excludes the hub. `fx` is the dimensional axial element force; the two
    blades are symmetric, so the blade-1 sum is doubled.
    """
    span = 1.0 - ROOT_CUTOUT_RADIUS / radius
    total = 0.0
    included = 0
    for path, rows in elements:
        window = [row for row in rows if start <= float(row["time"]) <= end]
        if not window:
            if not allow_short:
                raise ShortWindow(f"{path}: no samples in the averaging window")
            window = rows
        root_dist = mean([float(row["root_dist"]) for row in window])
        r_over_r = ROOT_CUTOUT_RADIUS / radius + root_dist * span
        if r_over_r >= 0.25:
            total += mean([float(row["fx"]) for row in window])
            included += 1
    return 2.0 * total, included, len(elements)


def analyse_model(run_dir: Path, cfg, speed, sequence, rho, allow_short,
                  thrust_scope, match_span=False, model=None):
    surface = model == "asm-mesh"
    turbine_path = run_dir / TURBINE_CSV
    element_dir = run_dir / ELEMENT_DIR
    if not turbine_path.is_file():
        raise FileNotFoundError(str(turbine_path))
    if not surface and not element_dir.is_dir():
        raise FileNotFoundError(str(element_dir))
    start, end, period = window_bounds(cfg, speed, sequence)
    rows = read_rows(turbine_path)
    if not rows:
        raise MissingInput(f"{turbine_path}: no data rows")
    entry = select_speed(cfg, speed, sequence)
    stats = analyse_turbine(
        rows, start, end, period, allow_short, float(entry["tsr"])
    )
    drift = drift_flags(rows, start, period)
    radius = float(cfg["turbine"]["radius"])
    if surface:
        # The station rows carry the element c_ref_n/c_ref_t definitions; the
        # only conversion is the element root_dist -> r/R mapping.
        profile = [
            (r_over_r_from_root_dist(root_dist, radius), c_ref_n, c_ref_t)
            for root_dist, c_ref_n, c_ref_t in read_surface_stations(
                run_dir, start, end, allow_short
            )
        ]
        profile.sort()
        # The element CSVs still exist (suppression stops only the strip
        # projection) and back element_profiles/--match-eaeroth-span.
        elements = read_elements(element_dir) if element_dir.is_dir() else []
    else:
        elements = read_elements(element_dir)
        profile = spanwise_profile(elements, start, end, allow_short, radius)
    metrics = turbine_metrics(stats, rho, speed, cfg, thrust_scope)
    metrics["run_dir"] = str(run_dir)
    metrics["element_profiles"] = len(elements)
    if surface:
        metrics["surface_stations"] = len(profile)
        metrics.update(read_surface_audit(run_dir))
    if match_span:
        if not elements:
            raise MissingInput(
                f"{element_dir}: element CSVs are required for "
                "--match-eaeroth-span"
            )
        span_thrust, included, total = span_limited_blade_thrust(
            elements, start, end, allow_short, radius
        )
        metrics["thrust_blade_span_n"] = span_thrust
        metrics["span_elements_included"] = included
        metrics["span_elements_total"] = total
        if metrics["thrust_scope"] == "blade":
            metrics["thrust_n"] = span_thrust
            metrics["thrust_scope"] = "blade (r/R >= 0.25)"
    return {
        "stats": stats,
        "drift": drift,
        "profile": profile,
        "spanwise": {
            f"{station:.2f}": interpolate(profile, station)
            for station in STATIONS
        },
        "metrics": metrics,
    }


def band_for(value, measured, relative):
    if measured == 0:
        return None
    return abs(value - measured) / abs(measured) / relative


def write_turbine_comparison(path, results, experiment, thrust_scope):
    measured = {
        "power": float(experiment["rotpow_kw"]),
        "torque": float(experiment["lsstqcor_nm"]),
        "thrust": float(experiment["eaeroth_n"]),
    }
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ["model", "quantity", "sim", "experiment", "delta_frac",
             "band_frac", "within_band", "units"]
        )
        units = {"power": "kW", "torque": "Nm", "thrust": "N"}
        for model in sorted(results):
            metrics = results[model]["metrics"]
            sim = {
                "power": metrics["power_kw"],
                "torque": metrics["torque_nm"],
                "thrust": metrics["thrust_n"],
            }
            for quantity in ("power", "torque", "thrust"):
                band = TURBINE_BANDS[quantity]
                fraction = band_for(sim[quantity], measured[quantity], band)
                writer.writerow([
                    model,
                    quantity,
                    f"{sim[quantity]:.8g}",
                    f"{measured[quantity]:.8g}",
                    f"{(sim[quantity] - measured[quantity]) / measured[quantity]:.6g}",
                    f"{fraction:.6g}" if fraction is not None else "",
                    "yes" if fraction is not None and fraction <= 1.0 else "no",
                    units[quantity],
                ])
    return measured


def write_spanwise_comparison(path, results, experiment_spanwise):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "model", "r_over_R", "c_ref_n", "cn_exp", "dcn_frac",
            "c_ref_t", "ct_exp", "dct_frac", "cm_exp_excluded",
            "within_band",
        ])
        for model in sorted(results):
            spanwise = results[model]["spanwise"]
            for station in STATIONS:
                values = spanwise.get(f"{station:.2f}")
                measured = experiment_spanwise[station]
                if values is None:
                    writer.writerow([model, f"{station:.2f}", "", "", "", "",
                                     "", "", f"{measured['cm']:.6g}", "no"])
                    continue
                c_ref_n, c_ref_t = values
                band_n = max(
                    SPANWISE_BAND_MIN,
                    SPANWISE_BAND_FRACTION * abs(measured["cn"]),
                )
                band_t = max(
                    SPANWISE_BAND_MIN,
                    SPANWISE_BAND_FRACTION * abs(measured["ct"]),
                )
                within = (
                    abs(c_ref_n - measured["cn"]) <= band_n
                    and abs(c_ref_t - measured["ct"]) <= band_t
                )
                writer.writerow([
                    model,
                    f"{station:.2f}",
                    f"{c_ref_n:.6g}",
                    f"{measured['cn']:.6g}",
                    f"{(c_ref_n - measured['cn']) / measured['cn']:.6g}",
                    f"{c_ref_t:.6g}",
                    f"{measured['ct']:.6g}",
                    f"{(c_ref_t - measured['ct']) / measured['ct']:.6g}",
                    f"{measured['cm']:.6g}",
                    "yes" if within else "no",
                ])


def sign_gate(results):
    """7 m/s sign verification: TSR, cp/cd positivity, spanwise c_ref signs."""
    report = {}
    failed = False
    for model in sorted(results):
        entry = results[model]
        metrics = entry["metrics"]
        spanwise = entry["spanwise"]
        checks = {
            "tsr_matches_config": abs(
                entry["stats"]["tsr"] - metrics["tsr_config"]
            ) <= TSR_TOLERANCE,
            "cp_positive": metrics["cp"] > 0.0,
            "cd_positive": metrics["cd"] > 0.0,
            "c_ref_n_positive": all(
                values is not None and values[0] > 0.0
                for values in spanwise.values()
            ),
            "c_ref_t_positive": all(
                values is not None and values[1] > 0.0
                for values in spanwise.values()
            ),
        }
        passed = all(checks.values())
        failed = failed or not passed
        report[model] = {"checks": checks, "pass": passed}
    report["pass"] = not failed
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alm-dir", type=Path, default=None)
    parser.add_argument("--asm-dir", type=Path, default=None)
    parser.add_argument("--asm-mesh-dir", type=Path, default=None,
                        help="ASM-mesh run directory (imported surface variant)")
    parser.add_argument("--run-dir", type=Path, default=None,
                        help="single run directory labelled as the model name")
    parser.add_argument("--experiment", type=Path, default=DEFAULT_EXPERIMENT)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--speed", type=float, default=7.0)
    parser.add_argument("--sequence", choices=("H", "S"), default="H")
    parser.add_argument(
        "--revolutions", nargs=2, type=float,
        metavar=("DISCARD", "END"), default=None,
    )
    parser.add_argument("--rho", default="auto",
                        help="'auto' or a value in kg/m^3")
    parser.add_argument("--thrust-scope", choices=("blade", "rotor"),
                        default="blade")
    parser.add_argument("--match-eaeroth-span", action="store_true",
                        help="restrict the blade thrust sum to r/R >= 0.25")
    parser.add_argument("--sign-gate", action="store_true")
    parser.add_argument("--allow-short-window", action="store_true")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args(argv)

    directories = {}
    if args.alm_dir:
        directories["alm"] = args.alm_dir
    if args.asm_dir:
        directories["asm"] = args.asm_dir
    if args.asm_mesh_dir:
        directories["asm-mesh"] = args.asm_mesh_dir
    if args.run_dir:
        directories[args.run_dir.name] = args.run_dir
    if not directories:
        return _exit("no run directory given (--alm-dir/--asm-dir/"
                     "--asm-mesh-dir/--run-dir)", EXIT_MISSING_INPUT)

    out_dir = args.out if args.out is not None else (
        ROOT / "results" / f"U{args.speed:g}-{args.sequence}"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        cfg = load_config(args.config)
        if args.revolutions is not None:
            cfg["solver"]["discard_revolutions"] = int(args.revolutions[0])
            cfg["solver"]["end_revolutions"] = int(args.revolutions[1])
        experiment, experiment_spanwise = load_experiment(
            args.experiment, args.sequence, args.speed
        )
    except (FileNotFoundError, KeyError, ValueError) as exc:
        return _exit(f"experiment/config input missing: {exc}", EXIT_MISSING_INPUT)

    rho_value = (
        float(experiment["rho_kg_m3"]) if args.rho == "auto" else float(args.rho)
    )
    rho_source = (
        f"WTBARO/WTATEMP ideal gas (experiment row {experiment['file_name']})"
        if args.rho == "auto" else "user override (--rho)"
    )

    results = {}
    exit_code = EXIT_OK
    for model, directory in sorted(directories.items()):
        try:
            results[model] = analyse_model(
                directory.resolve(), cfg, args.speed, args.sequence, rho_value,
                args.allow_short_window, args.thrust_scope,
                match_span=args.match_eaeroth_span, model=model,
            )
        except (FileNotFoundError, MissingInput) as exc:
            return _exit(f"missing input for {model}: {exc}", EXIT_MISSING_INPUT)
        except ShortWindow as exc:
            return _exit(f"{model}: short averaging window: {exc}",
                         EXIT_SHORT_WINDOW)
        except (TsrMismatch, ValueError) as exc:
            return _exit(f"{model}: {exc}", EXIT_TSR_MISMATCH)
        results[model]["metrics"]["tsr_config"] = float(
            select_speed(cfg, args.speed, args.sequence)["tsr"]
        )
        results[model]["metrics"]["rho_kg_m3"] = rho_value
        results[model]["metrics"]["rho_source"] = rho_source

    measured = write_turbine_comparison(
        out_dir / "turbine_comparison.csv", results, experiment,
        args.thrust_scope,
    )
    write_spanwise_comparison(
        out_dir / "spanwise_comparison.csv", results, experiment_spanwise
    )

    definitions = {
        "torque": "Q = ct * q_dyn * R with R = rotorRadius, "
                  "q_dyn = 1/2 rho A Uinf^2",
        "power": "P = cp * q_dyn * Uinf with cp = ct * TSR, "
                 "A = pi R^2",
        "thrust_blade": "T_blade = sum_blades cd_blade * q_dyn (blade-only)",
        "thrust_rotor": "T_rotor = cd * q_dyn (hub included; secondary)",
        "spanwise": "time-mean c_ref_n/c_ref_t interpolated to "
                    "30/47/63/80/95 % span",
        "density": "rho = WTBARO / (287.058 * (WTATEMP + 273.15))",
    }
    if "asm-mesh" in results:
        # Recorded only when a surface model is actually compared, so an
        # ALM/ASM-only metrics.json keeps its delivered definition block.
        definitions["surface_conversion"] = (
            "asm-mesh spanwise rows come from "
            "postProcessing/bladeSurface/*.csv and keep the element "
            "c_ref_n/c_ref_t definitions with the element root_dist -> r/R "
            "mapping"
        )

    metrics = {
        "case": {
            "speed_m_s": args.speed,
            "sequence": args.sequence,
            "tsr_configured": float(
                select_speed(cfg, args.speed, args.sequence)["tsr"]
            ),
            "rotor_radius_m": float(cfg["turbine"]["radius"]),
            "frontal_area_m2": math.pi * float(cfg["turbine"]["radius"]) ** 2,
            "experiment_row": experiment["file_name"],
            "rho_kg_m3": rho_value,
            "rho_source": rho_source,
            "thrust_scope": args.thrust_scope,
            "match_eaeroth_span": args.match_eaeroth_span,
        },
        "definitions": definitions,
        "bands": {
            "turbine_relative": TURBINE_BANDS,
            "spanwise": "max(0.15, 20 % of measured)",
            "drift_flag": DRIFT_TOLERANCE,
        },
        "experiment": {
            "file_name": experiment["file_name"],
            "rotpow_kw": measured["power"],
            "lsstqcor_nm": measured["torque"],
            "eaeroth_n": measured["thrust"],
            "wtbaro_pa": float(experiment["wtbaro_pa"]),
            "wtatemp_degc": float(experiment["wtatemp_degc"]),
            "spanwise": {
                f"{station:.2f}": experiment_spanwise[station]
                for station in STATIONS
            },
        },
        "window": {
            "discard_revolutions": int(cfg["solver"]["discard_revolutions"]),
            "end_revolutions": int(cfg["solver"]["end_revolutions"]),
            "allow_short_window": args.allow_short_window,
        },
        "models": {
            model: {
                "run_dir": results[model]["metrics"]["run_dir"],
                "metrics": results[model]["metrics"],
                "drift": results[model]["drift"],
                "revolutions_averaged": results[model]["stats"][
                    "revolutions_averaged"
                ],
                "window_mode": results[model]["stats"]["window_mode"],
                "spanwise": {
                    key: values for key, values in results[model]["spanwise"].items()
                },
            }
            for model in results
        },
        "cm_excluded": {
            "reason": "the element CSV has no `cm` column",
            "follow_up": "add momentCoefficient_ to the element writer in a "
                         "separate opt-in change",
            "measured_cm_stations": {
                f"{station:.2f}": experiment_spanwise[station]["cm"]
                for station in STATIONS
            },
        },
        "limitations": LIMITATIONS,
    }

    (out_dir / "metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    lines = [
        f"NREL Phase VI comparison - {args.sequence} sequence, "
        f"{args.speed:g} m/s",
        f"experiment row: {experiment['file_name']}",
        f"rho: {rho_value:.6f} kg/m^3 ({rho_source})",
        f"  source row: WTBARO {float(experiment['wtbaro_pa']):.4f} Pa, "
        f"WTATEMP {float(experiment['wtatemp_degc']):.4f} degC",
        f"averaging window: revolutions "
        f"{cfg['solver']['discard_revolutions']}.."
        f"{cfg['solver']['end_revolutions']} after discarding revolutions 0.."
        f"{cfg['solver']['discard_revolutions']}",
        "",
        "Metric definitions (fixed against the NREL report):",
        "  q_dyn = 1/2 rho A Uinf^2",
        "  Q = ct * q_dyn * R            (shaft torque; R = rotorRadius)",
        "  P = cp * q_dyn * Uinf         (cp = ct * TSR; A = pi R^2)",
        "  T_blade = sum_blades cd_blade * q_dyn (blade-only)",
        "  T_rotor = cd * q_dyn          (hub included; secondary)",
        "",
        "Turbine level (measured: ROTPOW / LSSTQCOR / EAEROTH):",
    ]
    for model in sorted(results):
        model_metrics = results[model]["metrics"]
        lines.append(
            f"  {model.upper()}: P {model_metrics['power_kw']:.4g} kW vs "
            f"{measured['power']:.4g} kW; Q {model_metrics['torque_nm']:.4g} Nm "
            f"vs {measured['torque']:.4g} Nm; T {model_metrics['thrust_n']:.4g} N "
            f"vs {measured['thrust']:.4g} N "
            f"(scope {model_metrics['thrust_scope']})"
        )
        drift = results[model]["drift"]
        lines.append(
            f"    drift flag: {drift['flag']} "
            f"(cp variation {drift['cp']['variation']}, "
            f"cd variation {drift['cd']['variation']})"
        )
    lines += [
        "",
        "Bands: turbine power/torque/thrust +-15 %; spanwise max(0.15, 20 %).",
        "Spanwise CM excluded (element CSV has no cm column).",
        "",
        "Limitations:",
    ]
    lines += [f"  - {item}" for item in LIMITATIONS]
    report = "\n".join(lines) + "\n"
    (out_dir / "report.txt").write_text(report, encoding="utf-8")

    gate = None
    if args.sign_gate:
        gate = sign_gate(results)
        (out_dir / "sign_gate.json").write_text(
            json.dumps(gate, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        if not gate["pass"]:
            print("sign gate FAIL", file=sys.stderr)
            return EXIT_MISSING_INPUT
        print("sign gate PASS")

    print(report)
    print(f"outputs written to {out_dir}")
    if gate is not None:
        print(f"sign gate: {'PASS' if gate['pass'] else 'FAIL'}")
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
