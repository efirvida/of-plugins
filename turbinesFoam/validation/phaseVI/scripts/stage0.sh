#!/bin/sh
# Stage 0 of the NREL Phase VI validation: mesh + checkMesh for D/32 and
# D/48, then the <= 0.3 revolution stability runs of ALM and ASM at 7 m/s on
# the authorized development queue. Stages 1-3 are declared in
# config/case.yaml but MUST NOT be submitted.
#
# Usage: stage0.sh [--domain long|squat] [--skip-fine-mesh] [--execute] [--submit]
#
#   (no flag)   print the stage plan and the prepared commands; nothing runs
#   --execute   run the Stage 0 payload (mesh coarse + fine, ALM + ASM runs);
#               intended to run on rank 0 of scripts/slurm/stage0.slurm
#   --submit    sbatch scripts/slurm/stage0.slurm on the development queue
#
# Exit codes as runPhaseVI.sh: 2 unsupported, 3 environment/mesh/solver,
# 4 stale case, 5 authorization gate.
set -eu

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(CDPATH= cd -- "$here/.." && pwd)

domain=long
skip_fine=0
execute=0
submit=0
while [ $# -gt 0 ]; do
    case "$1" in
        --domain) domain="$2"; shift ;;
        --skip-fine-mesh) skip_fine=1 ;;
        --execute) execute=1 ;;
        --submit) submit=1 ;;
        -h|--help) sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "ERROR: unsupported argument: $1" >&2; exit 2 ;;
    esac
    shift
done

echo "Phase VI staged plan (config/case.yaml):"
for stage in stage0 stage1 stage2 stage3; do
    printf '%s: ' "$stage"
    python3 "$root/tools/case_config.py" --stage "$stage" | tr '\n' ' '
    echo
done

if [ "$submit" -eq 1 ]; then
    echo "Submitting Stage 0 to the development queue (sequana_cpu_dev, 20 min)"
    sbatch "$here/slurm/stage0.slurm"
    exit 0
fi

if [ "$execute" -eq 0 ]; then
    echo "Prepared Stage 0 commands (not executed):"
    echo "  scripts/mesh.sh coarse --domain $domain"
    if [ "$skip_fine" -eq 0 ]; then
        echo "  scripts/mesh.sh fine --domain $domain"
    fi
    echo "  scripts/runPhaseVI.sh -m alm -u 7 -mesh coarse --domain $domain --stage0 --run"
    echo "  scripts/runPhaseVI.sh -m asm -u 7 -mesh coarse --domain $domain --stage0 --run"
    echo "Run with: scripts/stage0.sh --submit   (development queue only)"
    exit 0
fi

"$here/check_environment.sh" || exit $?

"$here/mesh.sh" coarse --domain "$domain"
if [ "$skip_fine" -eq 0 ]; then
    "$here/mesh.sh" fine --domain "$domain"
fi

"$here/runPhaseVI.sh" -m alm -u 7 -mesh coarse --domain "$domain" --stage0 --run
"$here/runPhaseVI.sh" -m asm -u 7 -mesh coarse --domain "$domain" --stage0 --run

echo "Stage 0 complete: mesh accepted for D/32 and D/48; ALM and ASM stability runs finished"
