#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Compare a periodic-nacelle actuator-surface run with the paper's reference.

Reads the run's probe time series and total-force CSV, forms the metrics the
`nacelle-validation-case` spec defines, and enforces the per-station/per-grid
acceptance recorded in `config/case.yaml`. No solver is required: the reference
profiles are the committed digitized CSVs under `data/reference/` and the run
output is read as plain text.

Metric definitions (fixed against the paper before any result is claimed):

    <u>(z)  = time mean of the streamwise velocity on the station line,
              normalised by U_inf                                [u/U_inf]
    k(z)    = 1/2 <u'_i u'_i> from the same series, the RESOLVED
              turbulent kinetic energy (the URANS fallback's modelled k
              field is not used), normalised by U_inf^2          [k/U_inf^2]
    CD      = |F_drag| / (0.5 rho U_inf^2 pi R^2)                [-]

with F_drag the time-mean streamwise force on the body (the source's `fx`),
rho the configured `inflow.density` and R the nacelle radius. The paper's
permeable-disk CD = 0.48 is recorded as a comparison DATUM, never as an
acceptance target.

Station reconciliation (documented in the package README): `config/case.yaml`
declares stations 1R/3R/5R/7R plus a 10R far-wake stretch, but the paper's
panels are at ODD multiples of R (1R..19R) and there is no 10R panel. A station
with an exact panel uses it; the 10R stretch is BRACKETED by the 9R and 11R
panels and its acceptance uses the worse of the two (the offsets are recorded
in `metrics.json`). No 10R reference profile is invented.

Averaging window: [wash_out_T * T_ft, (wash_out_T + average_T) * T_ft] from
`config/case.yaml`. The probes write every time step, so k uses the full 0.1 s
series rather than the 1 s field-write cadence.

Exit codes:
    0  comparison produced (per-station flags are reported, not hidden)
    1  missing input (run directory, probe file, force CSV, reference panel)
    2  averaging window too short (unless --allow-short-window)
    3  a claimed station failed its acceptance band
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from pathlib import Path

try:
    import yaml
except ImportError:  # pragma: no cover - the package tooling requires PyYAML
    yaml = None

SCRIPTS = Path(__file__).resolve().parent
ROOT = SCRIPTS.parent
DEFAULT_CONFIG = ROOT / "config" / "case.yaml"
DEFAULT_REFERENCE = ROOT / "data" / "reference"
DEFAULT_RESULTS = ROOT / "results"

FORCE_CSV = Path("postProcessing") / "nacelle" / "nacelle.csv"
PROBE_GLOB = "postProcessing/profiles/*/U"
PROBE_HEADER = re.compile(r"^#\s*Probe\s+(\d+)\s+\(([^)]*)\)")

# Acceptance bands: RMS deviation over the profile, in normalised units. They
# are package choices (the paper compares graphically), exposed as CLI
# overrides for calibration studies.
BAND_U_RMS = 0.10
BAND_K_RMS = 0.02

# A station line is matched to the nearest declared station within this
# distance (R): the rendered probes sit at the cell centre nearest the declared
# line, i.e. at most half the coarsest cell (0.5 * 30R/115 ~= 0.13R) away.
STATION_MATCH_TOL_R = 0.15

MIN_SAMPLES = 20
MIN_WINDOW_COVERAGE = 0.9

EXIT_OK = 0
EXIT_MISSING_INPUT = 1
EXIT_SHORT_WINDOW = 2
EXIT_ACCEPTANCE = 3

LIMITATIONS = [
    "The reference is a digitized figure, not raw data: the per-pixel "
    "resolution is ~0.003-0.007 u/U, ~1e-4 k/U^2 and ~0.004 z/D, and the "
    "committed points inherit the marker size and line thickness.",
    "The reference simulation is the paper's 502x348x348 wall-resolved LES "
    "with a dynamic SGS model that standard OpenFOAM does not ship; the case "
    "runs the documented WALE adaptation.",
    "The nacelle is represented only by the actuator surface: no body-fitted "
    "mesh and no resolved boundary layer, so only the wake profiles and the "
    "integral drag are claimed.",
    "The probes sit at the cell centres nearest the declared station lines "
    "(at most half a cell away, recorded per station) so that no probe lands "
    "on a cell face or a processor boundary.",
    "The paper's panels are at odd multiples of R: the 10R far-wake stretch "
    "has no panel and is bracketed by the 9R and 11R profiles (no 10R "
    "reference is invented).",
    "k is the RESOLVED turbulent kinetic energy from the 0.1 s probe series; "
    "it excludes the subgrid contribution and is not the URANS k field.",
    "The paper's permeable-disk CD = 0.48 is a comparison datum, not an "
    "acceptance target for the actuator-surface case.",
]


class MissingInput(Exception):
    pass


class ShortWindow(Exception):
    pass


def _exit(message: str, code: int) -> int:
    print(f"error: {message}", file=sys.stderr)
    return code


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


def load_config(path: Path) -> dict:
    if yaml is None:
        raise MissingInput("PyYAML is required to read config/case.yaml")
    with Path(path).open(encoding="utf-8") as stream:
        data = yaml.safe_load(stream)
    if not isinstance(data, dict):
        raise MissingInput(f"{path} does not contain a YAML mapping")
    return data


def metrics_stations_R(cfg: dict) -> list[float]:
    """Declared metric stations: stations_R plus the far-wake stretch_R."""
    return [float(value) for value in cfg["metrics"]["stations_R"]] + [
        float(cfg["metrics"]["stretch_R"])
    ]


def load_reference(reference_dir: Path) -> dict[float, dict[str, Path]]:
    """Digitized panel files by station: {station_R: {"u": path, "k": path}}."""
    panels: dict[float, dict[str, Path]] = {}
    for path in sorted(reference_dir.glob("u_*R.csv")):
        match = re.fullmatch(r"u_([0-9]+)R\.csv", path.name)
        if not match:
            continue
        station = float(match.group(1))
        panels[station] = {"u": path, "k": reference_dir / f"k_{match.group(1)}R.csv"}
    if not panels:
        raise MissingInput(f"no reference panels under {reference_dir}")
    for station, files in sorted(panels.items()):
        for kind, path in files.items():
            if not path.is_file():
                raise MissingInput(f"missing {kind} profile for {station:g}R: {path}")
    return panels


def read_reference_profile(path: Path, value_column: str) -> list[tuple[float, float]]:
    """(z_over_D, value) points of one committed reference CSV."""
    with path.open(encoding="utf-8") as stream:
        reader = csv.reader(stream)
        try:
            header = next(reader)
        except StopIteration:
            raise MissingInput(f"{path}: empty reference file") from None
        if [column.strip() for column in header] != ["z_over_D", value_column]:
            raise MissingInput(
                f"{path}: expected header ['z_over_D', '{value_column}'], "
                f"found {header}"
            )
        points = []
        for row in reader:
            if not row:
                continue
            points.append((float(row[0]), float(row[1])))
    if len(points) < 2:
        raise MissingInput(f"{path}: fewer than two reference points")
    points.sort()
    return points


def station_mapping(cfg: dict, panels: dict[float, dict[str, Path]]) -> dict:
    """Declared station -> exact panel or bracketing panel pair."""
    available = sorted(panels)
    mapping: dict[float, dict] = {}
    for station in metrics_stations_R(cfg):
        if station in panels:
            mapping[station] = {
                "mapping": "exact",
                "panels_R": [station],
                "offsets_R": [0.0],
            }
            continue
        lower = [panel for panel in available if panel < station]
        upper = [panel for panel in available if panel > station]
        if not lower or not upper:
            raise MissingInput(
                f"station {station:g}R has no exact or bracketing reference "
                f"panel (available: {available})"
            )
        bracket = [max(lower), min(upper)]
        mapping[station] = {
            "mapping": "bracket",
            "panels_R": bracket,
            "offsets_R": [panel - station for panel in bracket],
        }
    return mapping


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


def parse_probe_file(path: Path) -> tuple[list, list, list]:
    """(locations, times, rows) of one OpenFOAM `probes` output file."""
    locations: list[tuple[float, float, float]] = []
    times: list[float] = []
    rows: list[list[float]] = []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#"):
                match = PROBE_HEADER.match(line)
                if match:
                    coordinates = tuple(float(v) for v in match.group(2).split())
                    if len(coordinates) != 3:
                        raise MissingInput(f"{path}: malformed probe header {line!r}")
                    locations.append(coordinates)
                continue
            fields = line.replace("(", " ").replace(")", " ").split()
            if not fields:
                continue
            numbers = [float(field) for field in fields]
            times.append(numbers[0])
            rows.append(numbers[1:])
    return locations, times, rows


def read_probe_series(paths: list[Path], unset_limit: float = 1.0e30):
    """Merged (times, locations, series) of every probe file.

    `series[probe][index]` is the (ux, uy, uz) tuple at `times[index]`.
    A restarted run writes a new file under the new start-time directory, so
    all matching files are merged by time (first occurrence wins).
    """
    locations: list | None = None
    merged: dict[float, list[float]] = {}
    for path in paths:
        file_locations, times, rows = parse_probe_file(path)
        if not file_locations:
            raise MissingInput(f"{path}: no probe locations")
        if locations is None:
            locations = file_locations
        elif len(locations) != len(file_locations) or any(
            abs(a - b) > 1.0e-9
            for first, second in zip(locations, file_locations)
            for a, b in zip(first, second)
        ):
            raise MissingInput(f"{path}: probe locations differ from {paths[0]}")
        expected = 3 * len(file_locations)
        for time, row in zip(times, rows):
            if len(row) != expected:
                raise MissingInput(
                    f"{path}: row at t = {time:g} has {len(row)} values, "
                    f"expected {expected}"
                )
            if any(abs(value) > unset_limit for value in row):
                raise MissingInput(
                    f"{path}: a probe was not found in any cell at t = {time:g} "
                    "(unset value)"
                )
            if time not in merged:
                merged[time] = row
    if locations is None:
        raise MissingInput("no probe files found")
    times_sorted = sorted(merged)
    if len(times_sorted) < 2:
        raise MissingInput("fewer than two probe samples")
    columns = list(zip(*(merged[time] for time in times_sorted)))
    series = [
        list(zip(columns[3 * index], columns[3 * index + 1], columns[3 * index + 2]))
        for index in range(len(locations))
    ]
    return times_sorted, locations, series


def group_probes(
    locations: list, stations: list[float], tol: float
) -> tuple[dict[float, list[int]], dict[float, float]]:
    """Probe indices per declared station, and the mean x offset (R)."""
    groups: dict[float, list[int]] = {station: [] for station in stations}
    for index, location in enumerate(locations):
        nearest = min(stations, key=lambda station: abs(station - location[0]))
        if abs(nearest - location[0]) > tol:
            raise MissingInput(
                f"probe {index} at x = {location[0]:g}R is not within {tol:g}R "
                f"of any declared station"
            )
        groups[nearest].append(index)
    offsets = {}
    for station, indices in groups.items():
        if not indices:
            raise MissingInput(f"no probes were sampled at station {station:g}R")
        offsets[station] = (
            sum(locations[index][0] for index in indices) / len(indices)
        ) - station
    return groups, offsets


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


def window_bounds(cfg: dict) -> tuple[float, float]:
    """(start, end) of the averaging window, in seconds."""
    run = cfg["run"]
    t_ft = float(run["flow_through_time_s"])
    start = float(run["wash_out_T"]) * t_ft
    end = start + float(run["average_T"]) * t_ft
    return start, end


def mean_variance(values: list[float]) -> tuple[float, float]:
    count = len(values)
    mean = sum(values) / count
    variance = sum((value - mean) ** 2 for value in values) / count
    return mean, variance


def station_profile(
    series: list,
    indices: list[int],
    locations: list,
    velocity: float,
) -> list[tuple[float, float, float, int]]:
    """(z, <u>/U_inf, k/U_inf^2, samples) per probe, sorted by z.

    `series` is already restricted to the averaging window.
    """
    profile = []
    for index in indices:
        samples = series[index]
        count = len(samples)
        if count < 2:
            raise ShortWindow(
                f"probe {index} has {count} samples in the averaging window"
            )
        means = [0.0, 0.0, 0.0]
        variances = [0.0, 0.0, 0.0]
        for component in range(3):
            means[component], variances[component] = mean_variance(
                [value[component] for value in samples]
            )
        profile.append(
            (
                locations[index][2],
                means[0] / velocity,
                0.5 * sum(variances) / velocity**2,
                count,
            )
        )
    profile.sort(key=lambda point: point[0])
    return profile


def interpolate(xs: list[float], values: list[float], x: float) -> float:
    if x <= xs[0]:
        return values[0]
    if x >= xs[-1]:
        return values[-1]
    for (x0, v0), (x1, v1) in zip(zip(xs, values), zip(xs[1:], values[1:])):
        if x0 <= x <= x1:
            ratio = 0.0 if x1 == x0 else (x - x0) / (x1 - x0)
            return v0 + ratio * (v1 - v0)
    return values[-1]


def compare_station(
    station: float,
    profile: list[tuple[float, float, float, int]],
    mapping: dict,
    panels: dict[float, dict[str, Path]],
    radius: float,
    band_u: float,
    band_k: float,
) -> dict:
    """Pointwise deviations and RMS metrics against every mapped panel.

    A bracketed station (no exact panel) uses the WORSE of its two panels for
    the acceptance decision; both are reported.
    """
    zs = [point[0] for point in profile]
    us = [point[1] for point in profile]
    ks = [point[2] for point in profile]
    per_panel = []
    points = []
    for panel in mapping["panels_R"]:
        # The u and k figures are separate marker chains, so each quantity is
        # compared on its own reference grid.
        reference_u = read_reference_profile(panels[panel]["u"], "u_over_U")
        reference_k = read_reference_profile(panels[panel]["k"], "k_over_U2")
        deviations_u = []
        for z_over_d, u_ref in reference_u:
            z_r = z_over_d * 2.0 * radius
            u_run = interpolate(zs, us, z_r)
            deviations_u.append(u_run - u_ref)
            points.append(
                {
                    "station_panel_R": panel,
                    "quantity": "u",
                    "z_over_D": z_over_d,
                    "z_R": z_r,
                    "run": u_run,
                    "ref": u_ref,
                    "delta": u_run - u_ref,
                }
            )
        deviations_k = []
        for z_over_d, k_ref in reference_k:
            z_r = z_over_d * 2.0 * radius
            k_run = interpolate(zs, ks, z_r)
            deviations_k.append(k_run - k_ref)
            points.append(
                {
                    "station_panel_R": panel,
                    "quantity": "k",
                    "z_over_D": z_over_d,
                    "z_R": z_r,
                    "run": k_run,
                    "ref": k_ref,
                    "delta": k_run - k_ref,
                }
            )
        rms_u = math.sqrt(sum(value**2 for value in deviations_u) / len(deviations_u))
        rms_k = math.sqrt(sum(value**2 for value in deviations_k) / len(deviations_k))
        per_panel.append(
            {
                "station_panel_R": panel,
                "offset_R": panel - station,
                "reference_points_u": len(deviations_u),
                "reference_points_k": len(deviations_k),
                "rms_u": rms_u,
                "max_abs_u": max(abs(value) for value in deviations_u),
                "rms_k": rms_k,
                "max_abs_k": max(abs(value) for value in deviations_k),
                "within_band": rms_u <= band_u and rms_k <= band_k,
            }
        )
    return {
        "mapping": mapping["mapping"],
        "offsets_R": mapping["offsets_R"],
        "panels": per_panel,
        "points": points,
        "worst": {
            "rms_u": max(entry["rms_u"] for entry in per_panel),
            "rms_k": max(entry["rms_k"] for entry in per_panel),
            "within_band": all(entry["within_band"] for entry in per_panel),
        },
    }


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


def read_force_csv(path: Path) -> list[dict[str, float]]:
    with path.open(encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        rows = []
        for row in reader:
            parsed = {}
            for key, value in row.items():
                if key is None or value is None or value == "":
                    continue
                parsed[key.strip()] = float(value)
            if "time" not in parsed:
                raise MissingInput(f"{path}: no time column")
            rows.append(parsed)
    if len(rows) < 2:
        raise MissingInput(f"{path}: fewer than two force samples")
    return rows


def drag_coefficient(rows, cfg, start: float, end: float) -> dict:
    """CD from the mean streamwise force on the body (config definition)."""
    window = [row for row in rows if start <= row["time"] <= end]
    if not window:
        raise ShortWindow(
            f"no force samples in the averaging window [{start:g}, {end:g}] s"
        )
    radius = float(cfg["nacelle"]["radius"])
    velocity = float(cfg["inflow"]["velocity"])
    density = float(cfg["inflow"]["density"])
    mean_fx = sum(row["fx"] for row in window) / len(window)
    # The source's own cd column uses |F| and its referenceArea; kept as a
    # cross-check of the config definition above.
    if all("cd" in row for row in window):
        cd_from_csv = sum(row["cd"] for row in window) / len(window)
    else:
        cd_from_csv = None
    reference = 0.5 * density * velocity**2 * math.pi * radius**2
    return {
        "samples": len(window),
        "mean_fx_n": mean_fx,
        "mean_fy_n": sum(row["fy"] for row in window) / len(window),
        "mean_fz_n": sum(row["fz"] for row in window) / len(window),
        "cd": abs(mean_fx) / reference,
        "cd_from_csv": cd_from_csv,
        "reference_force_n": reference,
        "definition": "|F_drag| / (0.5 rho U_inf^2 pi R^2)",
    }


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


def acceptance_scope(cfg: dict, mesh: str) -> tuple[list[float], str]:
    if mesh == "medium":
        return metrics_stations_R(cfg), "acceptance.medium = all_stations"
    return [float(value) for value in cfg["acceptance"]["coarse"]], "acceptance.coarse"


def infer_mesh(run_dir: Path) -> str:
    name = run_dir.name.lower()
    for mesh in ("medium", "coarse"):
        if mesh in name:
            return mesh
    raise MissingInput(
        f"cannot infer the grid from '{run_dir.name}'; pass --mesh coarse|medium"
    )


def write_profiles_csv(path: Path, stations: dict, scope: list[float]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "station_R", "station_panel_R", "mapping", "claimed", "quantity",
            "z_over_D", "z_R", "run", "ref", "delta",
        ])
        for station in sorted(stations):
            entry = stations[station]
            claimed = station in scope
            for point in entry["points"]:
                writer.writerow([
                    f"{station:g}",
                    f"{point['station_panel_R']:g}",
                    entry["mapping"],
                    "yes" if claimed else "no",
                    point["quantity"],
                    f"{point['z_over_D']:.8g}",
                    f"{point['z_R']:.8g}",
                    f"{point['run']:.8g}",
                    f"{point['ref']:.8g}",
                    f"{point['delta']:.8g}",
                ])


def build_report(summary: dict) -> str:
    lines = [
        "Periodic-nacelle actuator-surface comparison",
        f"run directory: {summary['case']['run_dir']}",
        f"grid: {summary['case']['mesh']} ({summary['case']['cells']} cells)",
        f"reference: {summary['case']['reference']}",
        f"probes: {summary['case']['probe_count']} points, "
        f"{summary['case']['samples']} samples",
        f"averaging window: [{summary['window']['start_s']:g}, "
        f"{summary['window']['end_s']:g}] s "
        f"({summary['window']['t_ft_s']:g} s flow-through time)",
        "",
        "Metric definitions:",
        "  <u>(z)/U_inf  time-mean streamwise velocity on the station line",
        "  k(z)/U_inf^2  1/2 <u'_i u'_i> from the 0.1 s probe series (resolved)",
        "  CD            |F_drag| / (0.5 rho U_inf^2 pi R^2)",
        "",
        f"Acceptance bands: rms(u) <= {summary['bands']['u_rms']:g} u/U_inf, "
        f"rms(k) <= {summary['bands']['k_rms']:g} k/U_inf^2",
        f"Acceptance scope: {summary['acceptance']['source']}",
        "Station reconciliation: "
        + summary["station_reconciliation"]["note"],
        "",
        f"{'station':>8} {'mapping':>8} {'claimed':>8} {'rms u':>8} {'rms k':>8} "
        f"{'within band':>12}",
    ]
    claimed = summary["acceptance"]["claimed_R"]
    for key in sorted(summary["stations"], key=lambda name: float(name[:-1])):
        entry = summary["stations"][key]
        lines.append(
            f"{key:>8} {entry['mapping']:>8} "
            f"{('yes' if entry['claimed'] else 'no'):>8} "
            f"{entry['worst']['rms_u']:>8.4g} {entry['worst']['rms_k']:>8.4g} "
            f"{('yes' if entry['worst']['within_band'] else 'NO'):>12}"
        )
    drag = summary["drag"]
    lines += [
        "",
        f"CD = {drag['cd']:.6g} (from mean fx = {drag['mean_fx_n']:.6g} N over "
        f"{drag['samples']} samples)",
        f"  source CSV cd (cross-check): {drag['cd_from_csv']:.6g}"
        if drag["cd_from_csv"] is not None else "  source CSV cd: n/a",
        f"  paper permeable-disk datum: {summary['datum']['cd']:g} "
        "(not an acceptance target)",
        "",
        f"verdict: {summary['acceptance']['verdict']}",
        "",
        "Limitations:",
    ]
    lines += [f"  - {item}" for item in LIMITATIONS]
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--mesh", choices=("coarse", "medium"), default=None)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--band-u", type=float, default=BAND_U_RMS)
    parser.add_argument("--band-k", type=float, default=BAND_K_RMS)
    parser.add_argument("--allow-short-window", action="store_true")
    args = parser.parse_args(argv)

    run_dir = args.run_dir.resolve()
    if not run_dir.is_dir():
        return _exit(f"run directory not found: {run_dir}", EXIT_MISSING_INPUT)

    try:
        cfg = load_config(args.config)
        mesh = args.mesh or infer_mesh(run_dir)
        panels = load_reference(args.reference)
        mapping = station_mapping(cfg, panels)
        stations_declared = metrics_stations_R(cfg)
        start, end = window_bounds(cfg)
        velocity = float(cfg["inflow"]["velocity"])
        radius = float(cfg["nacelle"]["radius"])

        probe_paths = sorted(run_dir.glob(PROBE_GLOB))
        if not probe_paths:
            raise MissingInput(str(run_dir / PROBE_GLOB))
        all_times, locations, series = read_probe_series(probe_paths)
        groups, offsets = group_probes(locations, stations_declared, STATION_MATCH_TOL_R)

        force_path = run_dir / FORCE_CSV
        if not force_path.is_file():
            raise MissingInput(str(force_path))
        force_rows = read_force_csv(force_path)

        window_indices = [
            index for index, time in enumerate(all_times) if start <= time <= end
        ]
        coverage = (
            (all_times[window_indices[-1]] - all_times[window_indices[0]])
            / (end - start)
            if len(window_indices) >= 2
            else 0.0
        )
        if (
            len(window_indices) < MIN_SAMPLES
            or coverage < MIN_WINDOW_COVERAGE
        ):
            if not args.allow_short_window:
                raise ShortWindow(
                    f"{len(window_indices)} samples cover {coverage:.0%} of the "
                    f"configured window [{start:g}, {end:g}] s (need >= "
                    f"{MIN_SAMPLES} samples and >= {MIN_WINDOW_COVERAGE:.0%} "
                    "coverage)"
                )
            mode = "short"
            window_indices = window_indices or list(range(len(all_times)))
        else:
            mode = "full"

        window_times = [all_times[index] for index in window_indices]
        # Per probe, the samples inside the window (series is probe-major).
        window_series = [
            [series[probe][index] for index in window_indices]
            for probe in range(len(series))
        ]

        stations = {}
        for station in stations_declared:
            profile = station_profile(
                window_series, groups[station], locations, velocity
            )
            entry = compare_station(
                station, profile, mapping[station], panels, radius,
                args.band_u, args.band_k,
            )
            entry["x_offset_R"] = offsets[station]
            entry["probes"] = len(groups[station])
            entry["samples_in_window"] = len(window_times)
            stations[station] = entry

        scope, scope_source = acceptance_scope(cfg, mesh)
        claimed = sorted(station for station in stations if station in scope)
        failures = [
            f"{station:g}R"
            for station in claimed
            if not stations[station]["worst"]["within_band"]
        ]
        verdict = "PASS" if not failures else "FAIL: " + ", ".join(failures)

        drag = drag_coefficient(force_rows, cfg, start, end)

        out_dir = args.out if args.out is not None else DEFAULT_RESULTS / run_dir.name
        out_dir.mkdir(parents=True, exist_ok=True)

        cells = [int(value) for value in cfg["mesh"]["resolutions"][mesh]["cells"]]
        summary = {
            "case": {
                "run_dir": str(run_dir),
                "mesh": mesh,
                "cells": cells[0] * cells[1] * cells[2],
                "solver": str(cfg["solver"]["headline"]),
                "reference": "digitized wall-resolved LES "
                             "(Yang & Sotiropoulos, arXiv:1702.02108v4, Sec. 4.1)",
                "probe_files": [str(path.relative_to(run_dir)) for path in probe_paths],
                "probe_count": len(locations),
                "samples": len(all_times),
            },
            "window": {
                "start_s": start,
                "end_s": end,
                "t_ft_s": float(cfg["run"]["flow_through_time_s"]),
                "mode": mode,
                "samples_in_window": len(window_indices),
            },
            "definitions": {
                "u": "time-mean streamwise velocity on the station line / U_inf",
                "k": "1/2 <u'_i u'_i> from the probe series (resolved) / U_inf^2",
                "cd": "|F_drag| / (0.5 rho U_inf^2 pi R^2)",
                "rho_kg_m3": float(cfg["inflow"]["density"]),
                "u_inf_m_s": velocity,
                "radius_m": radius,
            },
            "bands": {"u_rms": args.band_u, "k_rms": args.band_k},
            "acceptance": {
                "scope": scope,
                "source": scope_source,
                "claimed_R": claimed,
                "verdict": verdict,
            },
            "datum": {
                "cd": float(cfg["metrics"]["permeable_disk_datum"]),
                "is_target": False,
            },
            "station_reconciliation": {
                "declared_stations_R": stations_declared,
                "paper_panels_R": sorted(panels),
                "note": "the paper has no 10R panel; the 10R stretch is "
                        "bracketed by the 9R and 11R panels and no 10R profile "
                        "is invented",
                "mapping": {
                    f"{station:g}R": {
                        "mapping": mapping[station]["mapping"],
                        "panels_R": mapping[station]["panels_R"],
                        "offsets_R": mapping[station]["offsets_R"],
                    }
                    for station in stations_declared
                },
            },
            "stations": {
                f"{station:g}R": {
                    "mapping": stations[station]["mapping"],
                    "offsets_R": stations[station]["offsets_R"],
                    "x_offset_R": stations[station]["x_offset_R"],
                    "probes": stations[station]["probes"],
                    "samples_in_window": stations[station]["samples_in_window"],
                    "panels": stations[station]["panels"],
                    "worst": stations[station]["worst"],
                    "claimed": station in scope,
                }
                for station in sorted(stations)
            },
            "drag": drag,
            "limitations": LIMITATIONS,
        }

        (out_dir / "metrics.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        write_profiles_csv(out_dir / "profiles.csv", stations, scope)
        report = build_report(summary)
        (out_dir / "report.txt").write_text(report, encoding="utf-8")

        print(report)
        print(f"outputs written to {out_dir}")
        if failures:
            print(f"acceptance FAIL at {', '.join(failures)}", file=sys.stderr)
            return EXIT_ACCEPTANCE
        return EXIT_OK
    except MissingInput as exc:
        return _exit(str(exc), EXIT_MISSING_INPUT)
    except ShortWindow as exc:
        return _exit(f"short averaging window: {exc}", EXIT_SHORT_WINDOW)
    except (KeyError, ValueError, IndexError) as exc:
        return _exit(f"malformed input: {exc!r}", EXIT_MISSING_INPUT)


if __name__ == "__main__":
    raise SystemExit(main())
