#!/bin/bash
# Stage 0 of the NREL Phase VI validation: mesh + checkMesh D/32 and D/48 plus
# the <= 0.3 revolution 7 m/s stability runs of ALM and ASM, on the authorized
# development queue (sequana_cpu_dev, 20 min). This is the only job of this
# change that may be submitted; Stages 1-3 stay prepared-only.
#
# Submit from anywhere with:
#   sbatch turbinesFoam/validation/phaseVI/scripts/slurm/stage0.slurm
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
. "$WM_PROJECT_DIR/etc/bashrc"

export FOAM_USER_LIBBIN=$WM_PROJECT_USER_DIR/platforms/$WM_OPTIONS/lib
export LD_LIBRARY_PATH=$FOAM_USER_LIBBIN:$WM_PROJECT_DIR/platforms/$WM_OPTIONS/lib:$WM_PROJECT_DIR/platforms/$WM_OPTIONS/lib/sys-openmpi:$LD_LIBRARY_PATH

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
job_dir="$(cd "$script_dir/../.." && pwd)"

echo "Stage 0 on $(hostname); submit dir: ${SLURM_SUBMIT_DIR:-unset}"
echo "Validation package: $job_dir"

if [ "${SLURM_PROCID:-0}" -eq 0 ]; then
    cd "$job_dir"
    "${script_dir}/../stage0.sh" --execute
fi

wait
echo "Stage 0 job finished"
