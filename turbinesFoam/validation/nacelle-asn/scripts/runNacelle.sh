#!/bin/sh
# SPDX-License-Identifier: GPL-3.0-or-later
# Prepare (and optionally run or submit) one periodic-nacelle validation run.
#
# Usage:
#   runNacelle.sh [--mesh coarse|medium] [--solver wale|urans] [--ranks N]
#                 [--stage0] [--end-time T] [--restart] [--run] [--submit]
#
# The case is rendered into runs/<id>/ from config/case.yaml (the committed
# case/ skeleton is never modified), a real 0/ time directory is installed from
# 0.org/, and any stale log.* is removed before running: OpenFOAM's
# runApplication (and any `[ -f log.x ]` guard) silently skips a step whose log
# already exists. Time directories are kept, so --restart resumes from the
# latest written time (falling back to startTime when nothing was written).
#
# --stage0 caps the run at a short stability window (2 s by default, --end-time
# overrides) for the authorized development queue. --run executes
# blockMesh/topoSet/checkMesh/decomposePar/mpirun inside the current
# allocation. --submit hands stage0 to Slurm; the production job is prepared
# only and can never be submitted by this script.
#
# Exit codes:
#   2  unsupported input (mesh, solver, flags, --ranks vs SLURM_NTASKS)
#   3  environment, mesh, decomposition or solver failure
#   4  generated case is stale
#   5  queue gate (long-queue authorization missing, or a production submission)
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(CDPATH= cd -- "$here/.." && pwd)

usage() {
    sed -n '2,26p' "$0" | sed 's/^# \{0,1\}//'
}

mesh=coarse
solver=wale
ranks_override=""
end_time=""
stage0=0
restart=0
run=0
submit=0
while [ $# -gt 0 ]; do
    case "$1" in
        --mesh) mesh="$2"; shift ;;
        --solver) solver="$2"; shift ;;
        --ranks) ranks_override="$2"; shift ;;
        --end-time) end_time="$2"; shift ;;
        --stage0) stage0=1 ;;
        --restart) restart=1 ;;
        --run) run=1 ;;
        --submit) submit=1 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "ERROR: unsupported argument: $1" >&2; usage >&2; exit 2 ;;
    esac
    shift
done

case "$mesh" in
    coarse|medium) ;;
    *) echo "ERROR: unsupported mesh: $mesh (choose coarse or medium)" >&2; exit 2 ;;
esac
case "$solver" in
    wale|urans) ;;
    *) echo "ERROR: unsupported solver: $solver (choose wale or urans)" >&2; exit 2 ;;
esac
if [ -n "$ranks_override" ]; then
    case "$ranks_override" in
        *[!0-9]*) echo "ERROR: --ranks must be a positive integer" >&2; exit 2 ;;
    esac
    if [ "$ranks_override" -le 0 ]; then
        echo "ERROR: --ranks must be a positive integer" >&2
        exit 2
    fi
fi
if [ -n "$end_time" ]; then
    case "$end_time" in
        *[!0-9.]*|*.*.*) echo "ERROR: --end-time must be a positive number" >&2; exit 2 ;;
    esac
fi
if [ "$run" -eq 1 ] && [ "$submit" -eq 1 ]; then
    echo "ERROR: --run and --submit are mutually exclusive" >&2
    exit 2
fi
if [ "$submit" -eq 1 ] && [ "$stage0" -ne 1 ]; then
    echo "ERROR: production is prepared only: this script submits stage0" >&2
    echo "  (sbatch scripts/slurm/stage0.slurm) and nothing else." >&2
    exit 5
fi

# The module environment is not enough to run anything: source the install's
# own bashrc (see the package README and the repository AGENTS.md).
version="${WM_PROJECT_VERSION:-}"
case "$version" in
    v2506|v2412) ;;
    *)
        echo "ERROR: load OpenFOAM v2506 or v2412 first (found '${version:-unset}')." >&2
        echo "  module load openfoam/v2506_openmpi-4.1.4_gnu" >&2
        echo "  source \$WM_PROJECT_DIR/etc/bashrc" >&2
        exit 3
        ;;
esac
for command in blockMesh topoSet checkMesh pimpleFoam decomposePar mpirun python3; do
    if ! command -v "$command" >/dev/null 2>&1; then
        echo "ERROR: required command not found: $command" >&2
        exit 3
    fi
done
library="${FOAM_USER_LIBBIN:-}/libturbinesFoam.so"
if [ ! -f "$library" ]; then
    echo "ERROR: libturbinesFoam.so not found in FOAM_USER_LIBBIN." >&2
    echo "  Build it with: cd turbinesFoam && ./Allwmake" >&2
    exit 3
fi

if ! python3 "$root/tools/generate_case.py" --check; then
    echo "ERROR: generated case/ is stale; re-render with tools/generate_case.py" >&2
    exit 4
fi

# The long queue is prepared only. Running inside a non-development allocation
# requires the explicit authorization the production job itself also checks.
partition="${SLURM_JOB_PARTITION:-${SLURM_PARTITION:-}}"
if [ -n "$partition" ] && [ "$partition" != "sequana_cpu_dev" ] \
        && [ "${NACELLE_LONG_QUEUE_AUTHORIZED:-0}" != "1" ]; then
    echo "ERROR: refusing to run on '$partition': the long queue is prepared" >&2
    echo "  only. Export NACELLE_LONG_QUEUE_AUTHORIZED=1 only after an explicit" >&2
    echo "  authorization to launch the production averaging." >&2
    exit 5
fi

config_ranks=$(python3 -c "import sys; sys.path.insert(0, '$root/tools'); import generate_case; print(generate_case.load_config()['decomposition']['number_of_subdomains'])")
ranks="$config_ranks"
if [ -n "$ranks_override" ]; then
    ranks="$ranks_override"
fi
# Inside a Slurm allocation the rendered decomposeParDict and `mpirun -np` must
# not ask for more ranks than sbatch allocated.
if [ -n "${SLURM_NTASKS:-}" ]; then
    case "$SLURM_NTASKS" in
        *[!0-9]*)
            echo "ERROR: SLURM_NTASKS is not a positive integer: $SLURM_NTASKS" >&2
            exit 2
            ;;
    esac
    if [ "$ranks" -ne "$SLURM_NTASKS" ]; then
        echo "ERROR: ranks $ranks does not match SLURM_NTASKS=$SLURM_NTASKS." >&2
        echo "  Re-run with --ranks $SLURM_NTASKS or fix the allocation." >&2
        exit 2
    fi
fi

if [ "$stage0" -eq 1 ] && [ -z "$end_time" ]; then
    end_time=2
fi

run_id="nacelle-$mesh"
if [ "$solver" = "urans" ]; then
    run_id="$run_id-urans"
fi
if [ "$stage0" -eq 1 ]; then
    run_id="$run_id-s0"
fi
run_dir="$root/runs/$run_id"

render_args="--mesh $mesh --solver $solver --case-dir $run_dir --ranks $ranks"
if [ -n "$end_time" ]; then
    render_args="$render_args --end-time $end_time"
fi
if [ "$restart" -eq 1 ]; then
    written_time=""
    if [ -d "$run_dir" ]; then
        # Any positive numeric time directory, including fractional times below
        # 1 s that a `[1-9]*` glob would miss; `processor*/<time>` directories
        # are covered by the depth-2 search.
        written_time=$(find "$run_dir" -maxdepth 2 -type d -name '[0-9]*' \
            2>/dev/null | awk -F/ '$NF + 0 > 0 {print; exit}')
    fi
    if [ -n "$written_time" ]; then
        render_args="$render_args --start-from latestTime"
    else
        echo "WARNING: --restart requested but no written time exists; using startTime" >&2
    fi
fi

echo "Preparing $run_id in $run_dir"
# shellcheck disable=SC2086
python3 "$root/tools/generate_case.py" $render_args

# Initial condition: the committed case keeps 0.org so that `--check` stays
# clean; the run directory gets a real 0/ time directory.
rm -rf "$run_dir/0"
cp -r "$run_dir/0.org" "$run_dir/0"

python3 - "$run_dir" "$root" "$mesh" "$solver" "$ranks" "$stage0" \
    "$end_time" "$restart" <<'PY'
import datetime
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

run_dir, root, mesh, solver, ranks, stage0, end_time, restart = sys.argv[1:9]
run_dir = Path(run_dir)
root = Path(root)


def git_commit():
    try:
        return subprocess.check_output(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            stderr=subprocess.DEVNULL,
        ).decode().strip()
    except Exception:
        return "unknown"


config = root / "config" / "case.yaml"
payload = {
    "git_commit": git_commit(),
    "openfoam_version": os.environ.get("WM_PROJECT_VERSION", "unknown"),
    "config_sha256": hashlib.sha256(config.read_bytes()).hexdigest(),
    "mesh": mesh,
    "solver": solver,
    "ranks": int(ranks),
    "stage0": stage0 == "1",
    "end_time_s": float(end_time) if end_time else None,
    "restart": restart == "1",
    "run_dir": str(run_dir),
    "start_time_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
}
(run_dir / "run.json").write_text(
    json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
print(f"Wrote {run_dir / 'run.json'}")
PY

if [ "$run" -eq 1 ]; then
    # A pre-existing log makes OpenFOAM's runApplication-style guards skip the
    # step silently; remove them so every requested step really executes.
    rm -f "$run_dir"/log.*
    echo "Removed stale logs (if any) before running"
    cd "$run_dir"
    for step in blockMesh topoSet checkMesh; do
        if ! "$step" > "log.$step" 2>&1; then
            tail -30 "log.$step" >&2
            echo "ERROR: $step failed; inspect $run_dir/log.$step" >&2
            exit 3
        fi
    done
    if ! decomposePar -force > log.decomposePar 2>&1; then
        tail -30 log.decomposePar >&2
        echo "ERROR: decomposePar failed" >&2
        exit 3
    fi
    echo "Running pimpleFoam on $ranks ranks"
    if ! mpirun -np "$ranks" pimpleFoam -parallel > log.pimpleFoam 2>&1; then
        tail -40 log.pimpleFoam >&2
        echo "ERROR: pimpleFoam failed; inspect $run_dir/log.pimpleFoam" >&2
        exit 3
    fi
    echo "Solver finished; log: $run_dir/log.pimpleFoam"
    echo "Compare with: python3 '$here/compareNacelle.py' --run-dir '$run_dir'"
elif [ "$submit" -eq 1 ]; then
    echo "Submitting the stage0 job (development queue)"
    sbatch "$here/slurm/stage0.slurm"
else
    suggestion="$0 --mesh $mesh --solver $solver --ranks $ranks"
    if [ "$stage0" -eq 1 ]; then
        suggestion="$suggestion --stage0"
    fi
    echo "Prepared $run_id. Run the solver with:"
    echo "  $suggestion --run"
fi
