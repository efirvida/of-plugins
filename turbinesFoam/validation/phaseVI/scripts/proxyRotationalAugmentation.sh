#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Serial driver for the committed 0.25-rev D/32 rotational-augmentation proxy.
# The development queue has a 20-minute time limit and MaxSubmit=1, so the
# harness runs the variants serially inside one job. Submit one chunk per job
# and chain the chunks (`SPEED:VARIANT` arguments), for example:
#
#   sbatch proxyRotationalAugmentation.sh 13:control 13:augmentation-on
#   sbatch proxyRotationalAugmentation.sh 13:augmentation-on-root-off 7:control
#   sbatch proxyRotationalAugmentation.sh 7:augmentation-on 7:augmentation-on-root-off
#
# The variants are never submitted as a Slurm array. Submitted only through the
# authorized harness path, which refuses any production queue. The production
# campaign stays prepared-only.
#
# Usage:
#   sbatch validation/phaseVI/scripts/proxyRotationalAugmentation.sh [SPEED:VARIANT ...]
#SBATCH -p sequana_cpu_dev
#SBATCH --ntasks=48
#SBATCH --nodes=1
#SBATCH -J phaseVI-proxy-ra
#SBATCH --time=00:20:00
#SBATCH -e %j.err
#SBATCH -o %j.out

set -eu

module purge
module load openfoam/v2506_openmpi-4.1.4_gnu gcc

# The module does not export the wmake identity variables and its
# LD_LIBRARY_PATH points at a dead Int32 directory; set them explicitly.
export WM_ARCH=linux64 WM_COMPILER=Gcc WM_COMPILE_OPTION=Opt
export WM_PRECISION_OPTION=DP WM_LABEL_SIZE=64
export WM_OPTIONS=linux64GccDPInt64Opt

set +eu
. "$WM_PROJECT_DIR/etc/bashrc"
set -eu

export FOAM_USER_LIBBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$WM_PROJECT_DIR/platforms/$WM_OPTIONS/lib:$WM_PROJECT_DIR/platforms/$WM_OPTIONS/lib/sys-openmpi:$LD_LIBRARY_PATH

# Resolve the harness directory. Under sbatch `$0` is Slurm's spooled copy, so
# fall back to `PROXY_SCRIPTS_DIR` (exported by the harness `--submit` path) or
# to the submitted script path reported by `scontrol`.
here=${PROXY_SCRIPTS_DIR:-$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)}
if [ ! -f "$here/proxyRotationalAugmentation.py" ]; then
    submitted=$(scontrol show job "$SLURM_JOB_ID" 2>/dev/null \
        | sed -n 's/^[[:space:]]*Command=//p')
    if [ -n "$submitted" ]; then
        here=$(CDPATH= cd -- "$(dirname -- "$submitted")" && pwd)
    fi
fi
if [ ! -f "$here/proxyRotationalAugmentation.py" ]; then
    echo "cannot locate proxyRotationalAugmentation.py;" \
         "set PROXY_SCRIPTS_DIR to its directory" >&2
    exit 1
fi

step_args=""
for step in "$@"; do
    step_args="$step_args --step $step"
done

# shellcheck disable=SC2086
python3 "$here/proxyRotationalAugmentation.py" --run --ranks 48 $step_args
