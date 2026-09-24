#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Committed 0.25-rev D/32 proxy harness for the rotational-augmentation fix.

Prepares and (inside an authorized Slurm allocation) runs three variants at 7
and 13 m/s on the development queue, then evaluates the fail-loud success
criteria:

    control                  rotational augmentation off, root effects on
    augmentation-on          rotational augmentation on,  root effects on
    augmentation-on-root-off rotational augmentation on,  root effects off

The variants share every other key; only the `rotationalAugmentation` switch
and the Glauert `rootEffects` setting change. Each variant runs in a
self-contained package copy (`system/fvOptions`, its own `constant/polyMesh`
link and a fresh `postProcessing/`) and the variants are chained **serially**,
never submitted as a Slurm array (the development queue has `MaxSubmit=1`).

The harness is a local / authorized-development-queue tool only. `--submit`
refuses any production queue, and no code path cancels or modifies a
production job: the production campaign stays prepared-only.

Usage:
    proxyRotationalAugmentation.py --check|--dry-run
    proxyRotationalAugmentation.py --prepare [--force]
    proxyRotationalAugmentation.py --run
    proxyRotationalAugmentation.py --evaluate
    proxyRotationalAugmentation.py --submit

`--check` renders the matrix in memory and prints the configuration without
submitting anything. `--prepare` renders the self-contained variant copies;
`--run` executes `decomposePar` + `mpirun` for each variant serially (inside a
Slurm allocation); `--evaluate` applies the control gate and the fail-loud
criteria to the collected outputs and exits non-zero on a miss.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

SCRIPTS = Path(__file__).resolve().parent
PACKAGE = SCRIPTS.parent
if str(PACKAGE / "tools") not in sys.path:
    sys.path.insert(0, str(PACKAGE / "tools"))

import generate_case  # noqa: E402
from case_config import (  # noqa: E402
    DEFAULT_CONFIG,
    ROOT_CUTOUT_RADIUS,
    load_config,
)

DEFAULT_WORK_ROOT = PACKAGE / "runs"
EXPERIMENT_CSV = PACKAGE / "data" / "experiment" / "sequence_H_performance.csv"

#: The committed proxy configuration: 0.25 revolutions on the D/32 coarse mesh
#: at 7 and 13 m/s, 48 ranks, on the authorized development queue.
SPEEDS = ("7", "13")
MESH = "coarse"
MESH_D_OVER = 32
REVOLUTIONS = 0.25
RUN_RANKS = 48
DEV_QUEUE = "sequana_cpu_dev"
SEQUENCE = "H"
MODEL = "alm"

#: Any queue containing one of these markers is a production queue and is
#: refused by `--submit` (prepared-only boundary).
PRODUCTION_QUEUE_MARKERS = ("long", "prod", "production")

VARIANTS: tuple[dict[str, Any], ...] = (
    {"name": "control", "augmentation": False, "root": True},
    {"name": "augmentation-on", "augmentation": True, "root": True},
    {"name": "augmentation-on-root-off", "augmentation": True, "root": False},
)

#: The Du-Selig constants are the paper defaults; no free-parameter calibration.
AUGMENTATION = {"model": "DuSelig", "a": 1, "b": 1, "d": 1}

#: Control-reproduction gate: the control U13 integrated `cp` must lie within
#: 15 % of the converged baseline `-0.0411` ([−0.0473, −0.0349]) before any
#: other variant is interpreted (DIAGNOSIS-2026-09-21.md).
CONTROL_CP_BASELINE = -0.0411
CONTROL_CP_TOLERANCE = 0.15
CONTROL_CP_RANGE = (-0.0473, -0.0349)

#: Secondary signal: the U13 mid-span `c_ref_t` must be clearly positive, above
#: the short-window drift tolerance (0.01), versus the control's `-0.063`.
MIDSPAN_STATION = 0.47
MIDSPAN_C_REF_T_MIN = 0.02
MIDSPAN_C_REF_T_CONTROL_REF = -0.063
DRIFT_TOLERANCE = 0.01

#: Primary signal: the U7 augmentation-on + root-off ablation must shrink the
#: converged -16 % power deficit toward or inside the ±15 % band.
U7_DEFICIT_BASELINE = -0.16
BAND = 0.15

TURBINE_CSV = Path("postProcessing") / "turbines" / "0" / "turbine.csv"
ELEMENT_DIR = Path("postProcessing") / "actuatorLineElements" / "0"


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


def variant_dir(work_root: Path, speed: str, variant: dict[str, Any]) -> Path:
    return work_root / f"proxy-U{speed}-{variant['name']}"


def base_run_dir(work_root: Path, speed: str) -> Path:
    return work_root / f"{MODEL}-U{speed}-{MESH}-s0"


def augmentation_block(active: bool) -> dict[str, Any]:
    return {**AUGMENTATION, "active": bool(active)}


def render_variant(cfg, speed, case_dir: Path, variant: dict[str, Any]) -> str:
    """Render the ALM `fvOptions` twin for one proxy variant."""
    return generate_case.render_fv_options(
        cfg,
        speed,
        MESH,
        case_dir,
        generate_case.ALM_ELEMENT,
        sequence=SEQUENCE,
        rotational_augmentation=augmentation_block(variant["augmentation"]),
        root_effects=variant["root"],
    )


def _strip_block(text: str, header: str) -> str:
    """Remove a ``header { ... }`` block with brace matching (line-based)."""
    lines = text.splitlines()
    out: list[str] = []
    depth = 0
    skipping = False
    for line in lines:
        if not skipping and line.strip() == header:
            skipping = True
            depth = 0
        if skipping:
            depth += line.count("{") - line.count("}")
            if depth <= 0 and "{" in line:
                skipping = False
            continue
        out.append(line)
    return "\n".join(out)


def toggle_stripped(text: str) -> list[str]:
    """`fvOptions` lines with the two toggle blocks removed.

    Two variants must be identical after this strip: the only differences
    allowed are the `rotationalAugmentation` switch and the `rootEffects`
    setting.
    """
    stripped = _strip_block(text, "rotationalAugmentation")
    return [
        line
        for line in stripped.splitlines()
        if not line.strip().startswith("rootEffects ")
    ]


def _extract_block(text: str, header: str) -> str:
    """Return the ``header { ... }`` block with brace matching (line-based)."""
    out: list[str] = []
    depth = 0
    collecting = False
    for line in text.splitlines():
        if not collecting and line.strip() == header:
            collecting = True
            depth = 0
        if collecting:
            out.append(line)
            depth += line.count("{") - line.count("}")
            if "{" in line and depth <= 0:
                break
    return "\n".join(out)


def assert_rendered_toggle(text: str, variant: dict[str, Any]) -> None:
    """The rendered `fvOptions` must carry exactly the variant's toggles."""
    expected_active = "active on;" if variant["augmentation"] else "active off;"
    # The rotationalAugmentation block is the only one that renders the
    # `active on|off;` + `model DuSelig;` pair.
    block = _extract_block(text, "rotationalAugmentation")
    if expected_active not in block or "model DuSelig;" not in block:
        raise AssertionError(
            f"variant {variant['name']!r}: rotationalAugmentation block does "
            f"not render {expected_active!r}"
        )
    expected_root = "rootEffects on;" if variant["root"] else "rootEffects off;"
    if expected_root not in text:
        raise AssertionError(
            f"variant {variant['name']!r}: expected {expected_root!r}"
        )


def variant_matrix(cfg, case_dir: Path | None = None) -> dict[tuple[str, str], str]:
    case_dir = case_dir if case_dir is not None else PACKAGE / "case"
    return {
        (speed, variant["name"]): render_variant(cfg, speed, case_dir, variant)
        for speed in SPEEDS
        for variant in VARIANTS
    }


def assert_variants_share_keys(cfg, case_dir: Path | None = None) -> None:
    """Fail when variants differ in anything but the two toggles."""
    matrix = variant_matrix(cfg, case_dir)
    for speed in SPEEDS:
        stripped = {
            variant["name"]: toggle_stripped(matrix[(speed, variant["name"])])
            for variant in VARIANTS
        }
        reference = stripped["control"]
        for name, lines in stripped.items():
            if lines != reference:
                raise AssertionError(
                    f"variant {name!r} at U{speed} differs from the control in a "
                    "key other than the augmentation/root toggles"
                )
        # The toggles themselves must actually be present and distinct.
        for variant in VARIANTS:
            assert_rendered_toggle(matrix[(speed, variant["name"])], variant)
        root_off = matrix[(speed, "augmentation-on-root-off")]
        if "rootEffects off;" not in root_off:
            raise AssertionError(
                f"the U{speed} root-off ablation does not render rootEffects off"
            )


def configuration(cfg, work_root: Path) -> dict[str, Any]:
    return {
        "case": str(PACKAGE / "case"),
        "work_root": str(work_root),
        "mesh": MESH,
        "mesh_cell_size_D": MESH_D_OVER,
        "speeds_m_s": list(SPEEDS),
        "revolutions": REVOLUTIONS,
        "ranks": RUN_RANKS,
        "queue": DEV_QUEUE,
        "model": MODEL,
        "sequence": SEQUENCE,
        "variants": [
            {
                "name": variant["name"],
                "rotationalAugmentation": "on" if variant["augmentation"] else "off",
                "rootEffects": "on" if variant["root"] else "off",
            }
            for variant in VARIANTS
        ],
        "control_baseline_cp": CONTROL_CP_BASELINE,
        "control_acceptance_range": list(CONTROL_CP_RANGE),
        "midspan_c_ref_t_min": MIDSPAN_C_REF_T_MIN,
        "u7_deficit_baseline": U7_DEFICIT_BASELINE,
        "band": BAND,
        "prepared_only": True,
        "production_submission": False,
    }


def execution_plan() -> dict[str, Any]:
    """The serial chain: one step per (speed, variant), never a Slurm array."""
    return {
        "mode": "serial",
        "array": False,
        "max_submit": 1,
        "steps": [
            {"speed": speed, "variant": variant["name"]}
            for speed in SPEEDS
            for variant in VARIANTS
        ],
    }


def assert_dev_queue(queue: str) -> None:
    lowered = queue.lower()
    if any(marker in lowered for marker in PRODUCTION_QUEUE_MARKERS):
        raise ValueError(
            f"refusing production queue {queue!r}: the proxy is prepared-only"
        )
    if queue != DEV_QUEUE:
        raise ValueError(
            f"refusing queue {queue!r}: the proxy runs only on {DEV_QUEUE!r}"
        )


def build_submit_command(job_script: Path, queue: str) -> list[str]:
    """`sbatch` command for the serial chain. Never a Slurm array."""
    assert_dev_queue(queue)
    return [
        "sbatch",
        f"--partition={queue}",
        f"--ntasks={RUN_RANKS}",
        str(job_script),
    ]


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


def _run(command: list[str], cwd: Path | None = None, check: bool = True):
    return subprocess.run(command, cwd=str(cwd) if cwd else None, check=check)


def prepare_base(work_root: Path, speed: str, ranks: int, force: bool) -> Path:
    """Prepare the committed 0.25-rev base case through `runPhaseVI.sh`."""
    base = base_run_dir(work_root, speed)
    if base.is_dir() and not force:
        return base
    if base.is_dir():
        shutil.rmtree(base)
    command = [
        str(SCRIPTS / "runPhaseVI.sh"),
        "-m", MODEL,
        "-u", speed,
        "-mesh", MESH,
        "--stage0",
        "--ranks", str(ranks),
    ]
    _run(command)
    if not base.is_dir():
        raise RuntimeError(f"runPhaseVI.sh did not prepare {base}")
    return base


def prepare(work_root: Path, cfg, ranks: int, force: bool = False) -> list[Path]:
    """Render the self-contained variant copies; never submits a job.

    Idempotent: an existing variant package is kept unless ``force`` is set, so
    a chunked serial run does not wipe the results of an earlier chunk.
    """
    prepared: list[Path] = []
    for speed in SPEEDS:
        base = prepare_base(work_root, speed, ranks, force)
        for variant in VARIANTS:
            vdir = variant_dir(work_root, speed, variant)
            if vdir.is_dir() and not force:
                prepared.append(vdir)
                continue
            if vdir.exists():
                shutil.rmtree(vdir)
            shutil.copytree(base, vdir, copy_function=os.link)
            post = vdir / "postProcessing"
            if post.exists():
                shutil.rmtree(post)
            post.mkdir(parents=True, exist_ok=True)
            fv_options = render_variant(cfg, speed, vdir, variant)
            # The package copy hardlinks the base files; unlink the twin first
            # so the variant's write does not silently rewrite the shared inode
            # (which would make every variant render the last write).
            fv_path = vdir / "system" / "fvOptions"
            if fv_path.exists():
                fv_path.unlink()
            fv_path.write_text(fv_options, encoding="utf-8")
            assert_rendered_toggle(fv_path.read_text(encoding="utf-8"), variant)
            manifest = {
                "speed_m_s": float(speed),
                "variant": variant["name"],
                "rotationalAugmentation": "on" if variant["augmentation"] else "off",
                "rootEffects": "on" if variant["root"] else "off",
                "mesh": MESH,
                "revolutions": REVOLUTIONS,
                "ranks": ranks,
            }
            (vdir / "proxy_variant.json").write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            prepared.append(vdir)
    return prepared


def parse_step(text: str) -> dict[str, str]:
    """Parse a ``SPEED:VARIANT`` step selector (dev-queue chunking)."""
    if ":" not in text:
        raise ValueError(f"step {text!r} must be SPEED:VARIANT")
    speed, _, name = text.partition(":")
    if speed not in SPEEDS:
        raise ValueError(f"step {text!r}: speed must be one of {list(SPEEDS)}")
    if name not in {variant["name"] for variant in VARIANTS}:
        raise ValueError(
            f"step {text!r}: variant must be one of "
            f"{[variant['name'] for variant in VARIANTS]}"
        )
    return {"speed": speed, "variant": name}


def selected_steps(steps: list[dict[str, str]] | None) -> list[dict[str, str]]:
    """The steps to run: the explicit chunk, or the full serial chain."""
    if steps:
        return list(steps)
    return execution_plan()["steps"]


def run_variants(
    work_root: Path,
    ranks: int | None = None,
    steps: list[dict[str, str]] | None = None,
) -> int:
    """Run each selected variant serially. Requires loaded OpenFOAM."""
    if ranks is None:
        ranks = int(os.environ.get("SLURM_NTASKS", RUN_RANKS))
    for tool in ("decomposePar", "mpirun"):
        if shutil.which(tool) is None:
            raise RuntimeError(f"{tool} not found: load the OpenFOAM environment")
    by_name = {variant["name"]: variant for variant in VARIANTS}
    failure = 0
    for step in selected_steps(steps):
        speed = step["speed"]
        variant = by_name[step["variant"]]
        vdir = variant_dir(work_root, speed, variant)
        if not (vdir / "system" / "fvOptions").is_file():
            raise RuntimeError(f"{vdir} is not prepared; run --prepare first")
        post = vdir / "postProcessing"
        if post.exists():
            shutil.rmtree(post)
        print(f"[proxy] {speed} m/s {variant['name']}: decomposePar", flush=True)
        _run(["decomposePar", "-force"], cwd=vdir)
        print(f"[proxy] {speed} m/s {variant['name']}: mpirun -np {ranks}",
              flush=True)
        result = _run(
            ["mpirun", "-np", str(ranks), "pimpleFoam", "-parallel"],
            cwd=vdir,
            check=False,
        )
        if result.returncode != 0:
            failure = 1
            print(
                f"[proxy] {speed} m/s {variant['name']}: solver failed "
                f"rc={result.returncode}",
                file=sys.stderr,
                flush=True,
            )
        # Keep the variant directory small: postProcessing is what the
        # evaluation reads, the decomposed mesh is not needed again.
        for processor in sorted(vdir.glob("processor*")):
            shutil.rmtree(processor, ignore_errors=True)
    return failure


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


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


def mean(values):
    values = [value for value in values if value is not None]
    if not values:
        return None
    return sum(values) / len(values)


def r_over_r(root_dist: float, radius: float) -> float:
    span = 1.0 - ROOT_CUTOUT_RADIUS / radius
    return ROOT_CUTOUT_RADIUS / radius + root_dist * span


def read_midspan_c_ref_t(run_dir: Path, station: float, radius: float):
    element_dir = run_dir / ELEMENT_DIR
    files = sorted(path for path in element_dir.glob("*.csv") if path.is_file())
    if not files:
        raise FileNotFoundError(str(element_dir / "*.csv"))
    profile = []
    for path in files:
        rows = read_rows(path)
        if not rows:
            continue
        root_dist = mean([float(row["root_dist"]) for row in rows])
        c_ref_t = mean([float(row["c_ref_t"]) for row in rows])
        profile.append((r_over_r(root_dist, radius), c_ref_t))
    profile.sort()
    if not profile:
        raise ValueError(f"{element_dir}: no element rows")
    if station <= profile[0][0]:
        return profile[0][1]
    if station >= profile[-1][0]:
        return profile[-1][1]
    for (x0, t0), (x1, t1) in zip(profile, profile[1:]):
        if x0 <= station <= x1:
            ratio = 0.0 if x1 == x0 else (station - x0) / (x1 - x0)
            return t0 + ratio * (t1 - t0)
    return None


def load_measured(speeds=SPEEDS) -> dict[str, dict[str, float]]:
    rows = read_rows(EXPERIMENT_CSV)
    measured: dict[str, dict[str, float]] = {}
    for speed in speeds:
        for row in rows:
            if math.isclose(float(row["wind_speed_m_s"]), float(speed), abs_tol=1e-9):
                measured[speed] = {
                    "rotpow_kw": float(row["rotpow_kw"]),
                    "lsstqcor_nm": float(row["lsstqcor_nm"]),
                    "rho_kg_m3": float(row["rho_kg_m3"]),
                }
                break
        else:
            raise KeyError(f"{EXPERIMENT_CSV}: no measured row for {speed} m/s")
    return measured


def build_metrics(work_root: Path, cfg, measured) -> dict[str, Any]:
    radius = float(cfg["turbine"]["radius"])
    area = math.pi * radius**2
    metrics: dict[str, Any] = {}
    for speed in SPEEDS:
        metrics[speed] = {}
        rho = measured[speed]["rho_kg_m3"]
        q_dyn = 0.5 * rho * area * float(speed) ** 2
        for variant in VARIANTS:
            vdir = variant_dir(work_root, speed, variant)
            rows = read_rows(vdir / TURBINE_CSV)
            if not rows:
                raise ValueError(f"{vdir / TURBINE_CSV}: no data rows")
            cp = mean([float(row["cp"]) for row in rows])
            ct = mean([float(row["ct"]) for row in rows])
            cd = mean([float(row["cd"]) for row in rows])
            midspan = read_midspan_c_ref_t(vdir, MIDSPAN_STATION, radius)
            metrics[speed][variant["name"]] = {
                "cp": cp,
                "ct": ct,
                "cd": cd,
                "power_kw": cp * q_dyn * float(speed) / 1.0e3,
                "torque_nm": ct * q_dyn * radius,
                "midspan_c_ref_t": midspan,
                "samples": len(rows),
            }
    return metrics


def check_gates(metrics, measured) -> tuple[list[str], dict[str, Any]]:
    """Evaluate the fail-loud criteria. Returns (failures, report)."""
    failures: list[str] = []
    report: dict[str, Any] = {"primary": {}, "secondary": {}, "control_gate": {}}

    control13 = metrics["13"]["control"]
    low, high = CONTROL_CP_RANGE
    control_ok = low <= control13["cp"] <= high
    report["control_gate"] = {
        "baseline_cp": CONTROL_CP_BASELINE,
        "acceptance_range": [low, high],
        "control_cp": control13["cp"],
        "pass": control_ok,
    }
    if not control_ok:
        failures.append(
            "control gate (U13 integrated cp): control cp "
            f"{control13['cp']:.5f} outside [{low:.4f}, {high:.4f}] "
            f"(converged baseline {CONTROL_CP_BASELINE:.4f})"
        )

    aug13 = metrics["13"]["augmentation-on"]
    report["primary"]["U13_augmentation_on"] = {
        "cp": aug13["cp"],
        "ct": aug13["ct"],
        "torque_nm": aug13["torque_nm"],
        "cp_positive": aug13["cp"] > 0.0,
        "ct_positive": aug13["ct"] > 0.0,
    }
    if not aug13["cp"] > 0.0:
        failures.append(
            f"U13 augmentation-on primary signal: cp is not positive "
            f"(cp={aug13['cp']:.5f})"
        )
    if not aug13["ct"] > 0.0:
        failures.append(
            f"U13 augmentation-on primary signal: ct is not positive "
            f"(ct={aug13['ct']:.5f})"
        )

    report["secondary"]["U13_augmentation_on"] = {
        "midspan_c_ref_t": aug13["midspan_c_ref_t"],
        "threshold": MIDSPAN_C_REF_T_MIN,
        "control_c_ref_t": control13["midspan_c_ref_t"],
        "control_reference": MIDSPAN_C_REF_T_CONTROL_REF,
        "pass": aug13["midspan_c_ref_t"] >= MIDSPAN_C_REF_T_MIN,
    }
    if not aug13["midspan_c_ref_t"] >= MIDSPAN_C_REF_T_MIN:
        failures.append(
            "U13 augmentation-on secondary signal: mid-span c_ref_t "
            f"{aug13['midspan_c_ref_t']:.5f} < {MIDSPAN_C_REF_T_MIN:.2f} "
            f"(control {control13['midspan_c_ref_t']:.5f}, reference "
            f"{MIDSPAN_C_REF_T_CONTROL_REF:.3f})"
        )

    measured7 = measured["7"]["rotpow_kw"]
    control7 = metrics["7"]["control"]
    ablation7 = metrics["7"]["augmentation-on-root-off"]
    deficit_control = (control7["power_kw"] - measured7) / measured7
    deficit_ablation = (ablation7["power_kw"] - measured7) / measured7
    shrunk = deficit_ablation > deficit_control
    inside = abs(deficit_ablation) <= BAND
    report["primary"]["U7_root_off_ablation"] = {
        "measured_power_kw": measured7,
        "control_power_kw": control7["power_kw"],
        "ablation_power_kw": ablation7["power_kw"],
        "control_deficit_frac": deficit_control,
        "ablation_deficit_frac": deficit_ablation,
        "converged_baseline_frac": U7_DEFICIT_BASELINE,
        "deficit_shrunk": shrunk,
        "inside_band": inside,
        "band": BAND,
    }
    if not shrunk:
        failures.append(
            "U7 augmentation-on + root-off primary signal: the power deficit "
            f"did not shrink (control {deficit_control:+.3%}, ablation "
            f"{deficit_ablation:+.3%}; baseline {U7_DEFICIT_BASELINE:+.0%})"
        )

    report["failures"] = failures
    report["pass"] = not failures
    return failures, report


def evaluate(work_root: Path, cfg, output: Path | None = None) -> int:
    measured = load_measured()
    metrics = build_metrics(work_root, cfg, measured)
    failures, report = check_gates(metrics, measured)
    report["metrics"] = metrics
    report["measured"] = measured
    report["configuration"] = configuration(cfg, work_root)
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(payload, encoding="utf-8")
    print(payload)
    for failure in failures:
        print(f"FAIL: {failure}", file=sys.stderr)
    return 0 if not failures else 1


# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--work-root", type=Path, default=DEFAULT_WORK_ROOT)
    parser.add_argument("--check", "--dry-run", dest="check", action="store_true",
                        help="render the matrix, print the configuration and "
                             "submit nothing")
    parser.add_argument("--prepare", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--evaluate", action="store_true")
    parser.add_argument("--submit", action="store_true",
                        help="submit the serial chain to the development queue")
    parser.add_argument("--queue", default=DEV_QUEUE)
    parser.add_argument("--ranks", type=int, default=RUN_RANKS)
    parser.add_argument("--step", action="append", default=None, metavar="SPEED:VARIANT",
                        help="run only this prepared step (repeatable, for "
                             "chunking across the 20-minute dev queue); the "
                             "default runs the full serial chain")
    parser.add_argument("--force", action="store_true",
                        help="re-prepare the base case even if it exists")
    parser.add_argument("--out", type=Path, default=None,
                        help="write the evaluation report to this path")
    args = parser.parse_args(argv)

    work_root = args.work_root.resolve()
    cfg = load_config(args.config)

    if args.check:
        assert_variants_share_keys(cfg)
        print(json.dumps(configuration(cfg, work_root), indent=2, sort_keys=True))
        print("check: variant matrix shares all keys but the two toggles")
        print("check: submitting nothing")
        return 0

    if args.submit:
        # The committed wrapper is the serial chain driver; the Python harness
        # never builds an array. Production queues are refused.
        job_script = SCRIPTS / "proxyRotationalAugmentation.sh"
        command = build_submit_command(job_script, args.queue)
        print(" ".join(command))
        return subprocess.run(command, check=False).returncode

    if args.prepare or args.run:
        steps = [parse_step(step) for step in args.step] if args.step else None
        if args.run:
            prepared = prepare(work_root, cfg, args.ranks, force=args.force)
            print(f"prepared {len(prepared)} variants under {work_root}")
            return run_variants(work_root, args.ranks, steps)
        prepared = prepare(work_root, cfg, args.ranks, force=args.force)
        print(f"prepared {len(prepared)} variants under {work_root}")
        return 0

    if args.evaluate:
        return evaluate(work_root, cfg, args.out)

    parser.error("choose one of --check/--dry-run, --prepare, --run, "
                 "--evaluate, --submit")


if __name__ == "__main__":
    raise SystemExit(main())
