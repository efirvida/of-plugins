#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# PREPARED ONLY - do not submit until the long-queue authorization is granted.
#
# Stage 1+2 production array of the NREL Phase VI validation: one array task per
# (model, speed, mesh) combination, 48 ranks, at most 24 h wall time, restart
# from the latest written time. 20 m/s is deliberately absent (renderable but
# not staged) so no blanket array can pick it up. Stage 3 (optional IDDES,
# nChordwise sweep, D/64) is intentionally not in this array.
#
# Guarded by PHASEVI_LONG_QUEUE_AUTHORIZED=1; runPhaseVI.sh --submit enforces
# the same gate.
#SBATCH -p sequana_cpu
#SBATCH --ntasks=48
#SBATCH --nodes=1
#SBATCH -J phaseVI-prod
#SBATCH --time=24:00:00
#SBATCH --array=0-13
#SBATCH -e %j.err
#SBATCH -o %j.out

set -eu

if [ "${PHASEVI_LONG_QUEUE_AUTHORIZED:-0}" != "1" ]; then
    echo "ERROR: long-queue authorization not granted; this file is prepared only." >&2
    exit 5
fi

module purge
module load openfoam/v2506_openmpi-4.1.4_gnu gcc

# The module exports Int32 paths; set the Int64 identity first and re-source the
# OpenFOAM bashrc so PATH points at the installed linux64GccDPInt64Opt binaries.
export WM_ARCH=linux64
export WM_COMPILER=Gcc
export WM_COMPILE_OPTION=Opt
export WM_PRECISION_OPTION=DP
export WM_LABEL_SIZE=64
export WM_OPTIONS=linux64GccDPInt64Opt
# etc/bashrc is not written for strict mode: it dereferences optional variables
# (WM_PROJECT_SITE on SDumont) and can return non-zero, so `set -eu` aborts the
# job while sourcing it. Disable both while sourcing, then restore them.
set +eu
. "$WM_PROJECT_DIR/etc/bashrc"
set -eu

export FOAM_USER_LIBBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$WM_PROJECT_DIR/platforms/$WM_OPTIONS/lib:$WM_PROJECT_DIR/platforms/$WM_OPTIONS/lib/sys-openmpi:$LD_LIBRARY_PATH

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Stage 1: 7 m/s headline D/32 then D/48, both models.
# Stage 2: {10,13,15,25} m/s at D/32 both models + Sequence S 7 m/s repeat.
combos=(
    "alm:7:coarse:H"
    "asm:7:coarse:H"
    "alm:7:fine:H"
    "asm:7:fine:H"
    "alm:10:coarse:H"
    "asm:10:coarse:H"
    "alm:13:coarse:H"
    "asm:13:coarse:H"
    "alm:15:coarse:H"
    "asm:15:coarse:H"
    "alm:25:coarse:H"
    "asm:25:coarse:H"
    "alm:7:coarse:S"
    "asm:7:coarse:S"
)

combo="${combos[${SLURM_ARRAY_TASK_ID:?array task id missing}]}"
IFS=':' read -r model speed mesh sequence <<< "$combo"

echo "Array task ${SLURM_ARRAY_TASK_ID}: $combo on $(hostname)"

"${script_dir}/../runPhaseVI.sh" -m "$model" -u "$speed" -mesh "$mesh" \
    -s "$sequence" --restart --run

echo "Array task ${SLURM_ARRAY_TASK_ID} finished"
