#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Stage 3 of the NREL Phase VI validation on the long queue: the optional
# IDDES confirmation (ALM and ASM on D/32 and D/48) plus the ASM nChordwise
# sweep on D/48 (1 and 3; the configured default 5 is the Stage 1 ASM run),
# all at 7 m/s, restarting from the latest written time. One array task per
# combination, six tasks, four nodes / 192 ranks each (48-core sequana nodes):
# the D/48 and IDDES cases need far more than the 48-rank Stage 1 allocation.
#
# Prepared only: refuses to run without PHASEVI_LONG_QUEUE_AUTHORIZED=1.
# Submit from the validation package directory (the package is resolved from
# $SLURM_SUBMIT_DIR because Slurm runs a spooled copy of this script):
#   cd turbinesFoam/validation/phaseVI
#   sbatch scripts/slurm/stage3.slurm
#SBATCH -p sequana_cpu
#SBATCH --ntasks=192
#SBATCH --nodes=4
#SBATCH -J phaseVI-stage3
#SBATCH --time=96:00:00
#SBATCH --array=0-5
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

# Array tasks: IDDES ALM/ASM on D/32 and D/48 (--solver iddes), then the ASM
# nChordwise 1 and 3 runs on D/48 (URANS, the Stage 1 baseline solver). Every
# task renders a 192-rank decomposeParDict and runs on the 192-rank allocation.
combos=(
    "alm:7:coarse:iddes:"
    "asm:7:coarse:iddes:"
    "alm:7:fine:iddes:"
    "asm:7:fine:iddes:"
    "asm:7:fine:urans:1"
    "asm:7:fine:urans:3"
)

combo="${combos[${SLURM_ARRAY_TASK_ID:?array task id missing}]}"
IFS=':' read -r model speed mesh solver nchordwise <<< "$combo"

echo "Array task ${SLURM_ARRAY_TASK_ID}: $combo on $(hostname)"

set -- -m "$model" -u "$speed" -mesh "$mesh" --solver "$solver" \
    --ranks 192 --restart --run
if [ -n "$nchordwise" ]; then
    set -- "$@" --nchordwise "$nchordwise"
fi

"${pkg_dir}/scripts/runPhaseVI.sh" "$@"

echo "Array task ${SLURM_ARRAY_TASK_ID} finished"
