#!/bin/bash -l
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -t 0-02:00:0

# quick usage (build/bin/dipole must already exist -- this script does NOT build):
#   sbatch run_HymnD_scale_variation_oberon.sh
#   SCALE_FACTORS="0.5 2.0" Y_VALS="-3.0 -2.5 2.5 3.0" sbatch run_HymnD_scale_variation_oberon.sh
#
# This is just the Oberon version of run_HymnD_scale_variation.sh. Same idea,
# just wrapped in Slurm so we can submit it without sitting on the login node.
#
# Both scale choices run in one job, which is fine because this is not a huge
# sweep. The main thing is to request a sensible number of CPUs so the loops
# don't get stuck running serially on whatever random number Slurm gives you.
#
# If we want extra rapidities later, we can just rerun with a different Y_VALS.
# run_many_Pb.sh only writes the y-values it is told to write, so it is safe.

set -euo pipefail

# See run_HymnD_members_oberon.sh's note: the compute node needs the same
# toolchain module as the login/build node, or build/bin/dipole's dynamic
# libstdc++ won't load.
module load GCCcore/13.3.0

CORES=${CORES:-${SLURM_CPUS_PER_TASK:-$(( $(nproc) - 2 ))}}
export CORES

echo "Starting on $(hostname) at $(date)"
echo "CORES=$CORES SCALE_FACTORS=${SCALE_FACTORS:-0.5 2.0} Y_VALS=${Y_VALS:-<run_many_Pb.sh default>}"

bash run_HymnD_scale_variation.sh

echo "Finished at $(date)"
