#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Stage 0 of the NREL Phase VI validation: mesh + checkMesh D/32 and D/48 plus
# the <= 0.3 revolution 7 m/s stability runs of ALM and ASM, on the authorized
# development queue (sequana_cpu_dev, 20 min). This is the only job of this
# change that may be submitted; Stages 1-3 stay prepared-only.
#
# Submit from the validation package directory (the package is resolved from
# $SLURM_SUBMIT_DIR because Slurm runs a spooled copy of this script):
#   cd turbinesFoam/validation/phaseVI
#   sbatch scripts/slurm/stage0.slurm [--mesh-only|--runs-only]
# or, equivalently, scripts/stage0.sh --submit.
#SBATCH -p sequana_cpu_dev
#SBATCH --ntasks=48
#SBATCH --nodes=1
#SBATCH -J phaseVI-stage0
#SBATCH --time=00:20:00
#SBATCH -e %j.err
#SBATCH -o %j.out

set -eu

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

# Slurm spools the submitted script into /var/spool/slurmd/<job>/ and executes
# that copy, so BASH_SOURCE[0] points outside the repository and cannot locate
# the validation package. $SLURM_SUBMIT_DIR is the directory sbatch was invoked
# from (the package dir when submitting from turbinesFoam/validation/phaseVI);
# PHASEVI_PKG_DIR overrides it for direct execution. Extra arguments are
# forwarded to stage0.sh (e.g. --mesh-only / --runs-only to split Stage 0 under
# the 20 min development-queue limit).
pkg_dir="${PHASEVI_PKG_DIR:-${SLURM_SUBMIT_DIR:-$PWD}}"
stage0="$pkg_dir/scripts/stage0.sh"
if [ ! -x "$stage0" ]; then
    echo "ERROR: Stage 0 driver not found: $stage0" >&2
    echo "  Submit from the validation package directory:" >&2
    echo "    cd turbinesFoam/validation/phaseVI" >&2
    echo "    sbatch scripts/slurm/stage0.slurm" >&2
    exit 3
fi

echo "Stage 0 on $(hostname); submit dir: ${SLURM_SUBMIT_DIR:-unset}"
echo "Validation package: $pkg_dir"

if [ "${SLURM_PROCID:-0}" -eq 0 ]; then
    cd "$pkg_dir"
    "$stage0" --execute "$@"
fi

wait
echo "Stage 0 job finished"
