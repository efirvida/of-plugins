#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Stage 3 D/64 sensitivity of the NREL Phase VI validation on the long queue:
# one array task per model (ALM, ASM), URANS 7 m/s on the ultra (D/64) mesh.
# Eight nodes / 384 ranks per task, 96 h: at the Stage 1 allocation (48 ranks)
# one D/64 run would need roughly 194 h, which fits comfortably in a single
# 96 h task at 384 ranks (8x). --restart lets a requeued array task continue
# from the latest written time if the wall time is ever exhausted.
#
# Prepared only: refuses to run without PHASEVI_LONG_QUEUE_AUTHORIZED=1.
# Submit from the validation package directory (the package is resolved from
# $SLURM_SUBMIT_DIR because Slurm runs a spooled copy of this script):
#   cd turbinesFoam/validation/phaseVI
#   sbatch scripts/slurm/stage3-d64.slurm
#SBATCH -p sequana_cpu
#SBATCH --ntasks=384
#SBATCH --nodes=8
#SBATCH -J phaseVI-stage3-d64
#SBATCH --time=96:00:00
#SBATCH --array=0-1
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

# Slurm executes a spooled copy of this script (/var/spool/slurmd/...), so
# BASH_SOURCE does not point into the repository; resolve the package from the
# submission directory instead (same fix as stage0.slurm/production.slurm).
pkg_dir="${PHASEVI_PKG_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
if [ ! -x "${pkg_dir}/scripts/runPhaseVI.sh" ]; then
    echo "ERROR: ${pkg_dir}/scripts/runPhaseVI.sh not found; submit from the package directory or set PHASEVI_PKG_DIR." >&2
    exit 3
fi

combos=(
    "alm:7:ultra"
    "asm:7:ultra"
)

combo="${combos[${SLURM_ARRAY_TASK_ID:?array task id missing}]}"
IFS=':' read -r model speed mesh <<< "$combo"

echo "Array task ${SLURM_ARRAY_TASK_ID}: $combo on $(hostname)"

"${pkg_dir}/scripts/runPhaseVI.sh" -m "$model" -u "$speed" -mesh "$mesh" \
    --ranks 384 --restart --run

echo "Array task ${SLURM_ARRAY_TASK_ID} finished"
