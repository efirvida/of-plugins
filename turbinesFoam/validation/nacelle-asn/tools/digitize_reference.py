#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Digitize the wall-resolved-LES reference profiles from the paper figures.

Source: Yang & Sotiropoulos, "A new class of actuator surface models for wind
turbines", arXiv:1702.02108v4 (2018), Fig. 5 (`u_nacelle`, vertical profiles of
time-averaged streamwise velocity) and Fig. 6 (`k_nacelle`, vertical profiles
of turbulence kinetic energy), Sec. 4.1. Both figures carry ten panels at
x = 1R, 3R, ..., 19R; the reference series in each panel is the **open-circle
marker** series (the paper's wall-resolved LES).

The figures ship with the arXiv LaTeX source as raster EPS converted to PDF
(`u_nacelle-eps-converted-to.pdf`, `k_nacelle-eps-converted-to.pdf`); the
embedded image is the highest-resolution copy available and is what this tool
digizes. `--pdf` extracts it (requires PyMuPDF); `--figure` uses an image that
was extracted earlier.

Method (recorded in `data/reference/PROVENANCE.md`):

* Panels are located from the frame lines. Each panel is calibrated twice and
  independently:
  - `tick`: the labelled axis values are matched to the major tick marks via
    the x-axis label centroids, and a least-squares line is fitted through
    the tick positions;
  - `label`: a least-squares line through the label centroids themselves.
* The circle markers are detected twice and independently:
  - `hole`: the marker interior is an enclosed white region, so the holes of
    the filled mask are labelled and near-circular holes of marker size give
    the centres;
  - `filter`: a ring template (annulus) is correlated with the dark mask and
    local maxima above a score threshold are taken.
* The committed profile is the union of both detections (deduplicated), and
  the JSON report records the cross-check between the two detectors and the
  two calibrations, so every committed value is traceable to the figure.

Requires numpy, Pillow and scipy (and PyMuPDF for `--pdf`) — digitization only;
the comparison tool and its tests stay stdlib-only and never re-digitize.

Usage:
    digitize_reference.py --pdf u=u_nacelle-eps-converted-to.pdf \
                          --pdf k=k_nacelle-eps-converted-to.pdf \
                          --out data/reference
    digitize_reference.py --figure u=u_nacelle.png --figure k=k_nacelle.png \
                          --out data/reference --overlay /tmp/overlay
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.signal import fftconvolve

TOOL_DIR = Path(__file__).resolve().parent
DEFAULT_OUT = TOOL_DIR.parent / "data" / "reference"

DARK = 128
R_OUT = 14.5          # marker ring outer radius [px]
R_IN = 9.5            # marker ring inner radius [px]
FILTER_THRESHOLD = 1.05
MATCH_TOL_PX = 4.0    # same-marker agreement distance between detectors
LABEL_MATCH_PX = 12.0  # label centroid -> major tick search radius

# Per figure: station labels (row-major over the 2x5 panel grid), the labelled
# x values and the labelled y values.
FIGURES = {
    "u": {
        "field": "u_over_U",
        "stations": ["1R", "3R", "5R", "7R", "9R",
                     "11R", "13R", "15R", "17R", "19R"],
        "x_values": [0.0, 0.5, 1.0],
        "y_values": [1.0, 0.0, -1.0],
    },
    "k": {
        "field": "k_over_U2",
        "stations": ["1R", "3R", "5R", "7R", "9R",
                     "11R", "13R", "15R", "17R", "19R"],
        "x_values": [0.0, 0.02, 0.04],
        "y_values": [1.0, 0.0, -1.0],
    },
}


# ---------------------------------------------------------------------------
# image utilities
# ---------------------------------------------------------------------------
def load_mask(path: Path) -> np.ndarray:
    return np.asarray(Image.open(path).convert("L")) < DARK


def group_runs(indices):
    out = []
    for value in indices:
        if out and value == out[-1][-1] + 1:
            out[-1].append(value)
        else:
            out.append([value])
    return [(int(run[0]), int(run[-1])) for run in out]


def long_lines(mask: np.ndarray, axis: int, frac: float = 0.5):
    counts = mask.sum(axis=axis)
    return group_runs(np.where(counts > frac * mask.shape[1 - axis])[0])


def find_panels(mask: np.ndarray):
    """The 2 x 5 panel grid from the frame lines."""
    rows = long_lines(mask, axis=1)
    cols = long_lines(mask, axis=0)
    if len(rows) != 4 or len(cols) != 10:
        raise SystemExit(f"unexpected frame layout: rows={rows} cols={cols}")
    ybands = [(rows[0], rows[1]), (rows[2], rows[3])]
    xbands = [(cols[2 * i], cols[2 * i + 1]) for i in range(5)]
    panels = []
    for row, (top_frame, bottom_frame) in enumerate(ybands):
        for col, (left_frame, right_frame) in enumerate(xbands):
            panels.append({
                "row": row, "col": col,
                # inner plot area (frame lines excluded)
                "x0": left_frame[1] + 1, "x1": right_frame[0] - 1,
                "y0": top_frame[1] + 1, "y1": bottom_frame[0] - 1,
                "x_frame": right_frame[0] - 1,   # last non-frame column
                "frame_left": left_frame[1],     # last column of the left frame
                "y_frame": bottom_frame[0] - 1,  # last non-frame row
                "frame_bottom": bottom_frame[0],  # first row of the bottom frame
            })
    return panels


def tick_runs(mask: np.ndarray, panel, axis: str):
    """Tick runs inward from the bottom (x) or left (y) frame.

    Returns `(centre, width, depth)` per run. A tick is narrow (<= 7 px); a
    data marker that touches the frame produces a much wider run and is
    rejected by the caller.
    """
    runs = []
    if axis == "x":
        probe = [(x, _depth_up(mask, panel["y_frame"], x))
                 for x in range(panel["x0"], panel["x_frame"] + 1)]
    else:
        probe = [(y, _depth_right(mask, panel["frame_left"] + 1, y))
                 for y in range(panel["y0"], panel["y_frame"] + 1)]
    for (a, b) in group_runs([pos for pos, depth in probe if depth >= 5]):
        depths = [depth for pos, depth in probe if a <= pos <= b]
        runs.append((0.5 * (a + b), b - a + 1, max(depths)))
    return runs


def _depth_up(mask, start_row, x, limit=40):
    depth = 0
    y = start_row
    while y >= 0 and mask[y, x] and depth < limit:
        depth += 1
        y -= 1
    return depth


def _depth_right(mask, start_col, y, limit=40):
    depth = 0
    x = start_col
    while x < mask.shape[1] and mask[y, x] and depth < limit:
        depth += 1
        x += 1
    return depth


def label_centres(mask: np.ndarray, panel):
    """x-axis label centroids below the panel frame (left to right).

    Digits of one label are merged; the decimal point and stray marks shorter
    than 20 px are dropped (the labels are ~36 px tall).
    """
    y0 = panel["frame_bottom"] + 6
    y1 = min(panel["frame_bottom"] + 60, mask.shape[0])
    xl = max(panel["x0"] - 45, 0)
    xr = min(panel["x1"] + 45, mask.shape[1])
    block = mask[y0:y1, xl:xr]
    labels, _ = ndimage.label(block)
    boxes = []
    for sl in ndimage.find_objects(labels):
        if sl is None:
            continue
        h = sl[0].stop - sl[0].start
        w = sl[1].stop - sl[1].start
        if h < 20 or w < 3:
            continue
        boxes.append((sl[1].start + xl, sl[1].stop + xl))
    boxes.sort()
    groups = []
    for box in boxes:
        if groups and box[0] - groups[-1][1] < 22:
            groups[-1] = (groups[-1][0], max(groups[-1][1], box[1]))
        else:
            groups.append(box)
    return [0.5 * (a + b) for a, b in groups]


def linear_fit(pixels, values, label):
    if len(pixels) < 2:
        raise ValueError(f"{label}: need >= 2 calibration points, "
                         f"got {len(pixels)}")
    matrix = np.vstack([np.array(pixels, dtype=float),
                        np.ones(len(pixels))]).T
    slope, intercept = np.linalg.lstsq(
        matrix, np.array(values, dtype=float), rcond=None)[0]
    residuals = [slope * px + intercept - value
                 for px, value in zip(pixels, values)]
    return {"slope": float(slope), "intercept": float(intercept),
            "pixels": [float(p) for p in pixels], "values": list(values),
            "max_abs_residual": float(max(abs(r) for r in residuals))}


def axis_calibration(mask, panel, axis: str, values):
    """Calibrate an axis from its major ticks + the label centroids.

    The labels give approximate pixel anchors for the labelled values; the
    exact pixels come from the nearest narrow major tick run. The returned
    label calibration (fit to the label centroids themselves) is the
    independent cross-check.
    """
    runs = tick_runs(mask, panel, axis)
    majors = [r for r in runs if r[1] <= 7 and r[2] >= 12]
    if axis == "x":
        centres = label_centres(mask, panel)
        if len(centres) != len(values):
            raise ValueError(f"expected {len(values)} x labels, "
                             f"got {len(centres)}")
        anchors = []
        used = set()
        unmatched = []
        for value, centre in zip(values, centres):
            best = None
            for index, run in enumerate(majors):
                if index in used:
                    continue
                distance = abs(run[0] - centre)
                if distance <= LABEL_MATCH_PX and (best is None or distance < best[0]):
                    best = (distance, index, run)
            if best is None:
                # a labelled value can sit on the panel edge (k figure: 0.04),
                # where the frame replaces the tick; the label still anchors
                # the value, but only the remaining ticks enter the fit.
                if not (centre < panel["x0"] + LABEL_MATCH_PX
                        or centre > panel["x1"] - LABEL_MATCH_PX):
                    raise ValueError(f"no major tick near label {value} "
                                     f"(centre {centre:.1f} px)")
                unmatched.append((centre, value))
                continue
            used.add(best[1])
            anchors.append((best[2][0], value))
        if len(anchors) < 2:
            raise ValueError(f"need >= 2 tick anchors, got {len(anchors)}")
        anchors.sort()
        tick_cal = linear_fit([a[0] for a in anchors],
                              [a[1] for a in anchors], "x tick")
        label_cal = linear_fit(centres, values, "x label")
        # every label must agree with the tick calibration to a few pixels
        for centre, value in zip(centres, values):
            residual_px = pixel_of(tick_cal, value) - centre
            if abs(residual_px) > 6.0:
                raise ValueError(f"label {value} disagrees with the ticks by "
                                 f"{residual_px:.1f} px")
        return tick_cal, label_cal
    # y axis: the labelled values are the only major ticks; require exactly the
    # expected count and a consistent spacing, otherwise refuse to guess.
    majors.sort(key=lambda r: r[0])
    if len(majors) != len(values):
        raise ValueError(f"expected {len(values)} y major ticks, "
                         f"got {len(majors)}: {[round(m[0], 1) for m in majors]}")
    tick_cal = linear_fit([m[0] for m in majors], values, "y tick")
    if tick_cal["max_abs_residual"] > 2.0:
        raise ValueError(f"y tick fit residual "
                         f"{tick_cal['max_abs_residual']:.2f} px")
    return tick_cal, None


def apply_cal(cal: dict, pixel: float) -> float:
    return float(cal["slope"] * pixel + cal["intercept"])


def pixel_of(cal: dict, value: float) -> float:
    return (value - cal["intercept"]) / cal["slope"]


# ---------------------------------------------------------------------------
# detectors
# ---------------------------------------------------------------------------
def ring_template():
    """Two normalised kernels: the marker annulus and the marker interior."""
    size = int(2 * R_OUT) + 9
    centre = (size - 1) / 2
    yy, xx = np.mgrid[0:size, 0:size]
    radius = np.hypot(xx - centre, yy - centre)
    annulus = ((radius >= R_IN) & (radius <= R_OUT)).astype(np.float64)
    disk = (radius <= R_IN - 1.5).astype(np.float64)
    return annulus / annulus.sum(), disk / disk.sum()


def detect_filter(mask: np.ndarray, panels, threshold=FILTER_THRESHOLD):
    """Ring matched filter: a dark annulus around a white interior.

    score = dark fraction on the annulus + 0.5 * white fraction on the disk,
    so a marker crossed by a curve still scores well (the interior stays
    mostly white), while a curve alone scores low (no white disk).
    """
    annulus, disk = ring_template()
    dark = mask.astype(np.float64)
    ring_fraction = fftconvolve(dark, annulus[::-1, ::-1], mode="same")
    interior_dark = fftconvolve(dark, disk[::-1, ::-1], mode="same")
    score = ring_fraction + 0.5 * (1.0 - interior_dark)
    detections = []
    for panel in panels:
        block = score[panel["y0"] + 6: panel["y1"] - 6,
                      panel["x0"] + 6: panel["x1"] - 6].copy()
        while True:
            index = np.unravel_index(np.argmax(block), block.shape)
            value = float(block[index])
            if value < threshold:
                break
            y, x = index
            detections.append({
                "panel": (panel["row"], panel["col"]),
                "xpix": float(x + panel["x0"] + 6),
                "ypix": float(y + panel["y0"] + 6),
                "score": float(value),
                "detector": "filter",
            })
            block[max(0, y - 16): y + 17, max(0, x - 16): x + 17] = -1e9
    return detections


def detect_holes(mask: np.ndarray, panels):
    """Marker interiors: enclosed white regions of marker size.

    Only clean interiors (area >= 200 px) count as marker centres; when the
    curve crosses a marker the interior is split and the pieces are rejected
    here — the matched filter detects those markers instead.
    """
    filled = ndimage.binary_fill_holes(mask)
    holes = filled & ~mask
    labels, _ = ndimage.label(holes)
    detections = []
    for index, sl in enumerate(ndimage.find_objects(labels)):
        if sl is None:
            continue
        h = sl[0].stop - sl[0].start
        w = sl[1].stop - sl[1].start
        if not (11 <= w <= 34 and 11 <= h <= 34):
            continue
        region = labels[sl] == index + 1
        area = int(region.sum())
        if not (200 <= area <= 560):
            continue
        fill = area / (math.pi * (min(w, h) / 2.0) ** 2)
        if not (0.5 <= fill <= 1.35):
            continue
        cy, cx = ndimage.center_of_mass(region)
        xpix, ypix = float(sl[1].start + cx), float(sl[0].start + cy)
        panel = None
        for candidate in panels:
            if (candidate["x0"] < xpix < candidate["x1"]
                    and candidate["y0"] < ypix < candidate["y1"]):
                panel = (candidate["row"], candidate["col"])
                break
        if panel is None:
            continue
        detections.append({"panel": panel, "xpix": float(xpix),
                           "ypix": float(ypix), "area": area,
                           "fill": float(fill), "detector": "hole"})
    return detections


def merge_detections(holes, filters):
    """Combine both detectors.

    The matched filter is the primary detector (it also finds markers crossed
    by a curve); a clean-interior hole that matches a filter detection within
    `MATCH_TOL_PX` confirms it (the centre becomes the mean of the two), and an
    unmatched hole is still kept. Everything else is a filter-only detection.
    """
    merged = []
    matched = []
    used_holes = set()
    for det in filters:
        x, y = det["xpix"], det["ypix"]
        best = None
        for hole in holes:
            if hole["panel"] != det["panel"]:
                continue
            distance = math.hypot(hole["xpix"] - x, hole["ypix"] - y)
            if best is None or distance < best[0]:
                best = (distance, hole)
        if best is not None and best[0] <= MATCH_TOL_PX:
            used_holes.add(id(best[1]))
            x = 0.5 * (x + best[1]["xpix"])
            y = 0.5 * (y + best[1]["ypix"])
            matched.append((best[1], det, best[0]))
        merged.append({"panel": det["panel"], "xpix": x, "ypix": y,
                       "detector": "filter"})
    extras = [hole for hole in holes if id(hole) not in used_holes]
    for hole in extras:
        if any(math.hypot(hole["xpix"] - m["xpix"], hole["ypix"] - m["ypix"])
               <= MATCH_TOL_PX for m in merged):
            continue
        merged.append({"panel": hole["panel"], "xpix": hole["xpix"],
                       "ypix": hole["ypix"], "detector": "hole"})
    return merged, matched, extras


def drop_isolated(detections, max_gap=50.0):
    """Remove detections that are not part of the marker chain.

    The marker series is a dense chain (~30 px spacing along the profile), so
    a detection with no neighbour within `max_gap` px is not a data point —
    in these figures it is the counter of the italic station label (e.g. the
    bowl of the "R" in "1R") caught by the ring filter.
    """
    kept, dropped = [], []
    for det in detections:
        nearest = min(
            (math.hypot(det["xpix"] - other["xpix"], det["ypix"] - other["ypix"])
             for other in detections if other is not det),
            default=float("inf"))
        if nearest <= max_gap:
            det = dict(det, nearest_neighbour_px=nearest)
            kept.append(det)
        else:
            dropped.append(dict(det, nearest_neighbour_px=None
                                if nearest == float("inf") else nearest))
    return kept, dropped


# ---------------------------------------------------------------------------
# per-figure pipeline
# ---------------------------------------------------------------------------
def digitize_figure(mask: np.ndarray, kind: str, overlay_dir: Path | None = None):
    spec = FIGURES[kind]
    panels = find_panels(mask)
    holes = detect_holes(mask, panels)
    filters = detect_filter(mask, panels)
    merged, matched, _ = merge_detections(holes, filters)
    merged, dropped = drop_isolated(merged)

    results = {}
    report = {
        "figure": kind,
        "detectors": {
            "hole": len(holes), "filter": len(filters),
            "matched": len(matched),
            "hole_only": len([d for d in merged if d["detector"] == "hole"]),
            "dropped_isolated": len(dropped),
        },
        "panels": [],
    }
    if dropped:
        report["dropped"] = [
            {"panel": list(d["panel"]), "xpix": round(d["xpix"], 1),
             "ypix": round(d["ypix"], 1)} for d in dropped
        ]
    if matched:
        distances = np.array([m[2] for m in matched])
        report["detector_agreement_px"] = {
            "matched": int(distances.size),
            "median": float(np.median(distances)),
            "max": float(distances.max()),
            "p95": float(np.percentile(distances, 95)),
        }

    for panel in panels:
        key = (panel["row"], panel["col"])
        station = spec["stations"][panel["row"] * 5 + panel["col"]]
        cal_x, cal_x_label = axis_calibration(mask, panel, "x", spec["x_values"])
        cal_y, _ = axis_calibration(mask, panel, "y", spec["y_values"])

        points = []
        for det in merged:
            if det["panel"] != key:
                continue
            points.append((apply_cal(cal_y, det["ypix"]),
                           apply_cal(cal_x, det["xpix"]),
                           det["detector"]))
        points.sort(key=lambda item: item[0])
        if not points:
            raise SystemExit(f"panel {station}: no markers detected")
        results[station] = points

        panel_report = {
            "station": station,
            "n_points": len(points),
            "n_hole": sum(1 for p in points if p[2] == "hole"),
            "n_filter": sum(1 for p in points if p[2] == "filter"),
            "z_over_D_range": [round(points[0][0], 5),
                               round(points[-1][0], 5)],
            "calibration": {"x_tick": cal_x, "x_label": cal_x_label,
                            "y_tick": cal_y},
            "points": [[round(z, 5), round(v, 6), source]
                       for z, v, source in points],
        }
        if cal_x_label is not None:
            # both calibrations evaluated at the panel centre and at the
            # extremes of the digitized profile
            for name, px in (("centre", 0.5 * (panel["x0"] + panel["x1"])),
                             ("left", panel["x0"]),
                             ("right", panel["x1"])):
                panel_report.setdefault("calibration_delta", {})[name] = (
                    apply_cal(cal_x, px) - apply_cal(cal_x_label, px))
        report["panels"].append(panel_report)

        if overlay_dir is not None:
            write_overlay(mask, panel, points, cal_x, cal_y,
                          overlay_dir / f"{kind}_{station}.png")
    return results, report


def write_overlay(mask, panel, points, cal_x, cal_y, path: Path):
    y0, y1 = panel["y0"] - 6, panel["y1"] + 6
    x0, x1 = panel["x0"] - 6, panel["x1"] + 6
    rgb = np.stack([mask[y0:y1, x0:x1]] * 3, axis=-1)
    image = (255 - rgb * 255).astype(np.uint8)
    for z_over_d, value, _ in points:
        xpix = round(pixel_of(cal_x, value)) - x0
        ypix = round(pixel_of(cal_y, z_over_d)) - y0
        for dx, dy in ((-1, 0), (1, 0), (0, -1), (0, 1), (0, 0)):
            xx, yy = xpix + dx, ypix + dy
            if 0 <= yy < image.shape[0] and 0 <= xx < image.shape[1]:
                image[yy, xx] = (255, 0, 0)
    path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(image).save(path)


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------
def parse_named(values, flag):
    out = {}
    for item in values:
        if "=" not in item:
            raise SystemExit(f"{flag} expects KIND=PATH, got {item!r}")
        kind, path = item.split("=", 1)
        if kind not in FIGURES:
            raise SystemExit(f"unknown figure kind {kind!r} (use u or k)")
        out[kind] = Path(path)
    return out


def extract_pdf(pdf: Path, dest: Path) -> Path:
    try:
        import fitz  # PyMuPDF
    except ImportError as exc:  # pragma: no cover - tooling path
        raise SystemExit(f"--pdf requires PyMuPDF: {exc}")
    document = fitz.open(pdf)
    images = document[0].get_images(full=True)
    if not images:
        raise SystemExit(f"no embedded image in {pdf}")
    info = document.extract_image(images[0][0])
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(info["image"])
    return dest


def write_profiles(out_dir: Path, kind: str, results: dict):
    spec = FIGURES[kind]
    written = []
    for station, points in sorted(results.items()):
        path = out_dir / f"{kind}_{station}.csv"
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.writer(stream)
            writer.writerow(["z_over_D", spec["field"]])
            for z_over_d, value, _ in points:
                writer.writerow([f"{z_over_d:.5f}", f"{value:.6f}"])
        written.append(path)
    return written


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--figure", action="append", default=[], metavar="KIND=PATH",
                        help="already-extracted figure image (u=... k=...)")
    parser.add_argument("--pdf", action="append", default=[], metavar="KIND=PATH",
                        help="arXiv source figure PDF to extract (u=... k=...)")
    parser.add_argument("--work-dir", type=Path, default=None,
                        help="where --pdf images are extracted (default: --out)")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT,
                        help="reference CSV output directory")
    parser.add_argument("--report", type=Path, default=None,
                        help="JSON cross-check report (default: <out>/digitization.json)")
    parser.add_argument("--overlay", type=Path, default=None,
                        help="write one QC overlay image per panel into this directory")
    args = parser.parse_args(argv)

    figures = dict(parse_named(args.figure, "--figure"))
    pdfs = parse_named(args.pdf, "--pdf")
    if not figures and not pdfs:
        parser.error("give --figure and/or --pdf (at least one per figure: u, k)")

    args.out.mkdir(parents=True, exist_ok=True)
    report = {"figures": {}, "outputs": {}}
    work = args.work_dir
    if work is None:
        work = Path(tempfile.mkdtemp(prefix="nacelle-digitize-"))
        print(f"extracted figures: {work}")
    for kind, pdf in pdfs.items():
        dest = work / f"{kind}_nacelle.png"
        extract_pdf(pdf, dest)
        figures[kind] = dest
        report["figures"][kind] = {
            "source_pdf": pdf.name,
            "source_pdf_sha256": sha256(pdf),
            "extracted_image": dest.name,
            "extracted_image_sha256": sha256(dest),
        }

    for kind, image_path in sorted(figures.items()):
        report["figures"].setdefault(kind, {})
        report["figures"][kind]["image"] = image_path.name
        report["figures"][kind]["image_sha256"] = sha256(image_path)
        mask = load_mask(image_path)
        results, figure_report = digitize_figure(mask, kind, args.overlay)
        report.setdefault("digitization", {})[kind] = figure_report
        for path in write_profiles(args.out, kind, results):
            report["outputs"][path.name] = {
                "sha256": sha256(path),
                "rows": sum(1 for _ in path.open(encoding="utf-8")) - 1,
            }
        print(f"{kind}: {len(results)} station profiles -> {args.out}")

    report_path = args.report or (args.out / "digitization.json")
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(f"report: {report_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
