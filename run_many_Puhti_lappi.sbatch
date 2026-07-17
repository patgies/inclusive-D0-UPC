#!/bin/bash -l
#SBATCH --job-name=incl_spectrum
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --output=logs/incl_%j.out
#SBATCH --error=logs/incl_%j.err

set -euo pipefail

# Configurable env vars (override when submitting)
OUTDIR=${OUTDIR:-/scratch/lappi/pagimeno/incl/proton}
# Use SLURM_CPUS_PER_TASK set by sbatch; fallback to 40
CORES=${SLURM_CPUS_PER_TASK:-${CORES:-40}}
PT_MIN=${PT_MIN:-0.1}
PT_STEP=${PT_STEP:-0.5}
PT_MAX=${PT_MAX:-12.0}
Y_VALS=${Y_VALS:-"0.0 1.0 2.0"}
REPO_DIP_DIR=${REPO_DIP_DIR:-data/Pb/mve}

mkdir -p logs
echo "Job $SLURM_JOB_ID starting on $(hostname)"
echo "OUTDIR=$OUTDIR CORES=$CORES PT=[$PT_MIN..$PT_MAX] PT_STEP=$PT_STEP Y_VALS=$Y_VALS"

mkdir -p "$OUTDIR"

# Use node-local storage (SLURM_TMPDIR) for dipole files
tmpdir=${SLURM_TMPDIR:-$(mktemp -d)}
mkdir -p "$tmpdir/mve"

echo "Copying dipole files to node-local storage: $tmpdir/mve"
rsync -a "${REPO_DIP_DIR}/" "$tmpdir/mve/"

function _copy_back() {
    echo "Copying outputs back to persistent OUTDIR: $OUTDIR"
    mkdir -p "$OUTDIR/files"
    rsync -a "$tmpdir/" "$OUTDIR/"
}
trap _copy_back EXIT

# Run from repository root
cd "$PWD"

# Invoke orchestrator using node-local dipole files
DIPOLE_DIR="$tmpdir/mve" \
OUTDIR="$OUTDIR" \
CORES="$CORES" \
PT_MIN="$PT_MIN" PT_STEP="$PT_STEP" PT_MAX="$PT_MAX" \
Y_VALS="$Y_VALS" \
bash run_many_Pb.sh

echo "Job $SLURM_JOB_ID finished"
