#!/bin/sh
# SPDX-License-Identifier: GPL-3.0-or-later
# Stage 0 of the NREL Phase VI validation: mesh + checkMesh for D/32 and
# D/48, then the <= 0.3 revolution stability runs of ALM and ASM at 7 m/s on
# the authorized development queue. Stages 1-3 are declared in
# config/case.yaml but MUST NOT be submitted.
#
# Usage: stage0.sh [--domain long|squat] [--skip-fine-mesh] [--mesh-only]
#                  [--runs-only] [--execute] [--submit]
#
#   (no flag)    print the stage plan and the prepared commands; nothing runs
#   --execute    run the Stage 0 payload (mesh coarse + fine, ALM + ASM runs);
#                intended to run on rank 0 of scripts/slurm/stage0.slurm
#   --mesh-only  with --execute, run only the mesh steps (coarse, plus fine
#                unless --skip-fine-mesh); used to split Stage 0 when the whole
#                payload does not fit in the 20 min development-queue job
#   --runs-only  with --execute, run only the ALM/ASM stability runs; requires
#                an already validated coarse mesh with the `turbine` cellSet
#   --submit     sbatch scripts/slurm/stage0.slurm on the development queue
#
# Exit codes as runPhaseVI.sh: 2 unsupported, 3 environment/mesh/solver,
# 4 stale case, 5 authorization gate.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(CDPATH= cd -- "$here/.." && pwd)

domain=long
skip_fine=0
mesh_only=0
runs_only=0
execute=0
submit=0
while [ $# -gt 0 ]; do
    case "$1" in
        --domain) domain="$2"; shift ;;
        --skip-fine-mesh) skip_fine=1 ;;
        --mesh-only) mesh_only=1 ;;
        --runs-only) runs_only=1 ;;
        --execute) execute=1 ;;
        --submit) submit=1 ;;
        -h|--help) sed -n '2,21p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "ERROR: unsupported argument: $1" >&2; exit 2 ;;
    esac
    shift
done

if [ "$mesh_only" -eq 1 ] && [ "$runs_only" -eq 1 ]; then
    echo "ERROR: --mesh-only and --runs-only are mutually exclusive" >&2
    exit 2
fi

echo "Phase VI staged plan (config/case.yaml):"
for stage in stage0 stage1 stage2 stage3; do
    printf '%s: ' "$stage"
    python3 "$root/tools/case_config.py" --stage "$stage" | tr '\n' ' '
    echo
done

if [ "$submit" -eq 1 ]; then
    echo "Submitting Stage 0 to the development queue (sequana_cpu_dev, 20 min)"
    # Submit from the package directory: the Slurm job resolves the package
    # from $SLURM_SUBMIT_DIR (the spooled script cannot locate it).
    cd "$root"
    sbatch "$here/slurm/stage0.slurm"
    exit 0
fi

if [ "$execute" -eq 0 ]; then
    echo "Prepared Stage 0 commands (not executed):"
    if [ "$runs_only" -eq 0 ]; then
        echo "  scripts/mesh.sh coarse --domain $domain"
        if [ "$skip_fine" -eq 0 ]; then
            echo "  scripts/mesh.sh fine --domain $domain"
        fi
    fi
    if [ "$mesh_only" -eq 0 ]; then
        echo "  scripts/runPhaseVI.sh -m alm -u 7 -mesh coarse --domain $domain --stage0 --run"
        echo "  scripts/runPhaseVI.sh -m asm -u 7 -mesh coarse --domain $domain --stage0 --run"
    fi
    echo "Run with: scripts/stage0.sh --submit   (development queue only)"
    echo "Split under the 20 min job limit:"
    echo "  scripts/stage0.sh --mesh-only --execute"
    echo "  scripts/stage0.sh --runs-only --execute"
    exit 0
fi

"$here/check_environment.sh" || exit $?

if [ "$runs_only" -eq 1 ]; then
    coarse_dir="$root/runs/mesh-coarse"
    if [ "$domain" = "squat" ]; then
        coarse_dir="$coarse_dir-squat"
    fi
    if [ ! -f "$coarse_dir/log.checkMesh" ] || ! grep -q "Mesh OK" "$coarse_dir/log.checkMesh"; then
        echo "ERROR: --runs-only requires a validated coarse mesh" >&2
        echo "  missing or rejected: $coarse_dir/log.checkMesh" >&2
        exit 3
    fi
    if [ ! -f "$coarse_dir/constant/polyMesh/sets/turbine" ]; then
        echo "ERROR: --runs-only requires the 'turbine' cellSet in $coarse_dir" >&2
        echo "  run the mesh step first: scripts/stage0.sh --mesh-only --execute" >&2
        exit 3
    fi
fi

if [ "$runs_only" -eq 0 ]; then
    "$here/mesh.sh" coarse --domain "$domain"
    if [ "$skip_fine" -eq 0 ]; then
        "$here/mesh.sh" fine --domain "$domain"
    fi
fi

if [ "$mesh_only" -eq 0 ]; then
    "$here/runPhaseVI.sh" -m alm -u 7 -mesh coarse --domain "$domain" --stage0 --run
    "$here/runPhaseVI.sh" -m asm -u 7 -mesh coarse --domain "$domain" --stage0 --run
fi

if [ "$mesh_only" -eq 1 ]; then
    echo "Stage 0 mesh step complete: D/32 and D/48 validated"
elif [ "$runs_only" -eq 1 ]; then
    echo "Stage 0 stability step complete: ALM and ASM 7 m/s D/32 finished"
else
    echo "Stage 0 complete: mesh accepted for D/32 and D/48; ALM and ASM stability runs finished"
fi
