#!/bin/bash -l
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -t 0-02:00:0

# Usage (build/bin/dipole must already exist -- this script does NOT build):
#   sbatch run_lhapdf_scale_variation_oberon.sh
#   SCALE_FACTORS="0.5 2.0" Y_VALS="-3.0 -2.5 2.5 3.0" sbatch run_lhapdf_scale_variation_oberon.sh
#
# SLURM wrapper around run_lhapdf_scale_variation.sh, same role as
# run_many_oberon.sh plays for run_many_Pb.sh -- this one just adds the
# SBATCH header and CORES handling that plain script doesn't have. Both
# scale factors run sequentially in one job (each factor is already a
# single-member run, not a 100+-member sweep, so there's no need for an
# array job here).
#
# NOTE: -c 20 above requests 20 CPUs per task so that SLURM_CPUS_PER_TASK is
# set correctly; without it, CORES falls back to "nproc - 2" (the whole
# node), which mismatches the real allocation if fewer CPUs were granted
# (see run_lhapdf_members_oberon.sh's note -- this bit us once there).
#
# To backfill extra rapidities into an already-completed scale-variation
# set (run_many_Pb.sh only (re)writes the y-tagged files it's asked for, so
# this is safe and won't touch/duplicate already-finished y values):
#   Y_VALS="-3.0 -2.5 2.5 3.0" sbatch run_lhapdf_scale_variation_oberon.sh

set -euo pipefail

# See run_lhapdf_members_oberon.sh's note: the compute node needs the same
# toolchain module as the login/build node, or build/bin/dipole's dynamic
# libstdc++ won't load.
module load GCCcore/13.3.0

CORES=${CORES:-${SLURM_CPUS_PER_TASK:-$(( $(nproc) - 2 ))}}
export CORES

echo "Starting on $(hostname) at $(date)"
echo "CORES=$CORES SCALE_FACTORS=${SCALE_FACTORS:-0.5 2.0} Y_VALS=${Y_VALS:-<run_many_Pb.sh default>}"

bash run_lhapdf_scale_variation.sh

echo "Finished at $(date)"
