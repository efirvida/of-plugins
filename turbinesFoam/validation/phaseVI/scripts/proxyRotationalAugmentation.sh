#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Serial driver for the committed 0.25-rev D/32 rotational-augmentation proxy.
# The development queue has MaxSubmit=1, so all three variants at both speeds
# run inside this single job, chained serially; it is never submitted as an
# array. `proxyRotationalAugmentation.py --run` prepares each self-contained
# variant copy and executes `decomposePar` + `mpirun` for it in order.
#
# Submitted only through the authorized harness path, which refuses any
# production queue. The production campaign stays prepared-only.
#
# Usage:
#   sbatch validation/phaseVI/scripts/proxyRotationalAugmentation.sh
#SBATCH -p sequana_cpu_dev
#SBATCH --ntasks=48
#SBATCH --nodes=1
#SBATCH -J phaseVI-proxy-ra
#SBATCH --time=01:30:00
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

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)

python3 "$here/proxyRotationalAugmentation.py" --run --ranks 48
