#!/bin/bash -l
#SBATCH --job-name=D0_incl_workflow
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=1000M
#SBATCH --output=slurm-%x-%j.out

# Runs one workflow script as a job. Submit from the repo root, after build_roihu.sh, e.g.
#   FRAG_TYPE=KniehlKramer OUTBASE=/scratch/lappi/patricia/D0_kk_band \
#       sbatch --job-name=kk_scale run_scripts/roihu/run_workflow_roihu.sbatch run_scripts/run_scale_variation.sh

set -e

if [[ $# -lt 1 ]]; then
	echo "Usage: sbatch [options] run_workflow_roihu.sbatch <workflow_script>" >&2
	exit 1
fi

cd "$SLURM_SUBMIT_DIR"
source run_scripts/roihu/roihu_modules.sh
load_roihu_modules
export CORES="$SLURM_CPUS_PER_TASK"

echo "Job $SLURM_JOB_ID ($SLURM_JOB_NAME): $1, CHANNEL=${CHANNEL:-default} FRAG_TYPE=${FRAG_TYPE:-default} OUTDIR=${OUTDIR:-default} OUTBASE=${OUTBASE:-default}, $(date)"
bash "$1"
echo "Job $SLURM_JOB_ID done, $(date)"
