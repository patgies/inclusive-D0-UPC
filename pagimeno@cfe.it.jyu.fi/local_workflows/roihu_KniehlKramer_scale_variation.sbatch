#!/bin/bash
#SBATCH --job-name=kk-scale-band      # Job name
#SBATCH --account=lappi               # CSC billing project
#SBATCH --partition=small             # Single-node partition (up to 384 cores/node on Roihu)
#SBATCH --time=08:00:00               # Max. duration of the job -- tune once you've timed one factor
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64            # Parallel dipole runs; raise up to 384 if the queue allows
#SBATCH --mem-per-cpu=1000M
#SBATCH --output=slurm-%j.out

# Run the Kniehl & Kramer factorisation-scale band (factors 0.5, 1.0, 2.0) on
# Roihu. Wraps the existing local_workflows/run_KniehlKramer_scale_variation.sh
# (itself a loop over run_many_Pb.sh) with the module environment and SLURM
# core count.
#
# Usage: sbatch local_workflows/roihu_KniehlKramer_scale_variation.sbatch
#
# Assumes the repo has already been built on Roihu:
#   module load spack/x86_64/v2026_03/Core/gcc/15.2.0
#   module load spack/x86_64/v2026_03/gcc/15.2.0/gsl/2.8
#   module load spack/x86_64/v2026_03/Core/cmake/3.31.11
#   mkdir build && cd build && cmake .. && make -j
#
# Run from the repo root (sbatch inherits the submission directory as $SLURM_SUBMIT_DIR).

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

source /usr/share/lmod/lmod/init/bash
export MODULEPATH=/appl/modulefiles
module load spack/x86_64/v2026_03/Core/gcc/15.2.0
module load spack/x86_64/v2026_03/gcc/15.2.0/gsl/2.8

export CORES="$SLURM_CPUS_PER_TASK"
export OUTBASE=${OUTBASE:-/scratch/lappi/patricia/D0_kk_band}
export SCALE_FACTORS=${SCALE_FACTORS:-"0.5 1.0 2.0"}

mkdir -p "$OUTBASE"

# run_KniehlKramer_scale_variation.sh calls "bash run_many_Pb.sh" (no
# directory prefix), and run_many_Pb.sh in turn opens build/bin/dipole,
# data/Pb/mve, QCDnum/timelike.wgt and inputs/Gamma_AA.dat all as paths
# relative to the repo root. So everything has to run with CWD at the repo
# root, and run_many_Pb.sh (which only physically lives in local_workflows/)
# needs to be reachable there too -- hence the symlink.
ln -sfn local_workflows/run_many_Pb.sh run_many_Pb.sh

bash local_workflows/run_KniehlKramer_scale_variation.sh
