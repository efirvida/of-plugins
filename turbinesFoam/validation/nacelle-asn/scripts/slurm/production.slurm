#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Production array of the periodic-nacelle validation: the full wash-out plus
# time-averaging run of both run grids (coarse 153x80x80, medium 115x151x151),
# 48 ranks, restarting from the latest written time.
#
# PREPARED ONLY: this job MUST NOT be submitted before an explicit
# long-queue authorization. `executes: false` in config/case.yaml,
# runNacelle.sh refuses a non-development partition unless
# NACELLE_LONG_QUEUE_AUTHORIZED=1, and this script checks the same gate.
#
#   cd turbinesFoam/validation/nacelle-asn
#   sbatch scripts/slurm/production.slurm        # only after authorization
#
# After the array finishes, compare each grid with:
#   python3 scripts/compareNacelle.py --run-dir runs/nacelle-coarse --mesh coarse
#   python3 scripts/compareNacelle.py --run-dir runs/nacelle-medium --mesh medium
#SBATCH -p sequana_cpu_long
#SBATCH --ntasks=48
#SBATCH --nodes=1
#SBATCH -J nacelle-prod
#SBATCH --time=96:00:00
#SBATCH --array=0-1
#SBATCH -e %j.err
#SBATCH -o %j.out

set -eu

if [ "${NACELLE_LONG_QUEUE_AUTHORIZED:-0}" != "1" ]; then
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
# submission directory instead (same fix as stage0.slurm).
pkg_dir="${NACELLE_PKG_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
if [ ! -x "${pkg_dir}/scripts/runNacelle.sh" ]; then
    echo "ERROR: ${pkg_dir}/scripts/runNacelle.sh not found; submit from the package directory or set NACELLE_PKG_DIR." >&2
    exit 3
fi

grids=(coarse medium)
grid="${grids[${SLURM_ARRAY_TASK_ID:?array task id missing}]}"

echo "Array task ${SLURM_ARRAY_TASK_ID}: ${grid} grid on $(hostname)"

"${pkg_dir}/scripts/runNacelle.sh" --mesh "$grid" --run --restart --ranks "${SLURM_NTASKS}"

echo "Array task ${SLURM_ARRAY_TASK_ID} finished"
