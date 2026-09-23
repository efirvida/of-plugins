#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# ASM-mesh array of the NREL Phase VI validation (S2 three-way comparison:
# ASM over the imported blade surface). One array task per (model, speed, mesh)
# combination, 48 ranks, at most 24 h wall time per task, restarting from the
# latest written time. Task 0 is the D/32 measurement gate; task 1 is the D/48
# headline run.
#
# PREPARED ONLY. This file refuses to run without
# PHASEVI_LONG_QUEUE_AUTHORIZED=1 and is never submitted by this change.
# runPhaseVI.sh --submit refuses -m asm-mesh and points here.
#
# Gate procedure (design sections 6.5 and 8.4): task 0 (asm-mesh:7:coarse:H,
# the D/32 mesh) is the performance measurement gate. It must be executed only
# under explicit authorization and its per-addSup instrumentation reviewed
# BEFORE any D/48 preparation or submission. The instrumentation is emitted by
# bladeSurfaceSource: one Info line per distribute() on the master rank and the
# per-call rows in postProcessing/bladeSurface/<owner>.surface_distribution.csv
# (time,nodes,candidates,mean_candidates,max_candidates,seconds); it reports
# the surface node count, the bounded candidate entries visited and the wall
# seconds per addSup. The gate decision must be recorded, choosing one of:
#   1. proceed   - the measured per-step cost is within the practical budget;
#   2. harden    - the measured cost exceeds the budget: harden the candidate
#                  query (W4.4, mesh-tree/stencil) and re-measure task 0;
#   3. restrict  - keep the campaign scope (mesh/speed ladder) that the user
#                  authorizes; the scope choice stays with the user
#                  (design section 14 Q2).
# The D/32 task fits in one 24 h allocation; the D/48 headline task may need
# several `--restart` requeues. The measurement is a performance gate only:
# the sub-grid caveat and the kernel/width confound framing are unchanged.
# The existing phaseVI-prod, phaseVI-stage3 and phaseVI-stage3-d64 arrays are
# read-only baselines: not submitted, cancelled or modified.
#
# PREPARED ONLY: this file must NOT be submitted by this change. The D/32
# measurement is prepared here and executed only under explicit authorization;
# no automated step of this change calls sbatch.
#
# Submit from the validation package directory (the package is resolved from
# $SLURM_SUBMIT_DIR because Slurm runs a spooled copy of this script):
#   cd turbinesFoam/validation/phaseVI
#   sbatch scripts/slurm/asm-mesh.slurm
#SBATCH -p sequana_cpu
#SBATCH --ntasks=48
#SBATCH --nodes=1
#SBATCH -J phaseVI-asm-mesh
#SBATCH --time=24:00:00
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

# Task 0: D/32 (coarse) measurement gate - execute and review before task 1.
# Task 1: D/48 (fine) headline run.
combos=(
    "asm-mesh:7:coarse:H"
    "asm-mesh:7:fine:H"
)

combo="${combos[${SLURM_ARRAY_TASK_ID:?array task id missing}]}"
IFS=':' read -r model speed mesh sequence <<< "$combo"

echo "Array task ${SLURM_ARRAY_TASK_ID}: $combo on $(hostname)"

"${pkg_dir}/scripts/runPhaseVI.sh" -m "$model" -u "$speed" -mesh "$mesh" \
    -s "$sequence" --restart --run

echo "Array task ${SLURM_ARRAY_TASK_ID} finished"
