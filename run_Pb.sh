#SBATCH --job-name=incl_spectrum_Pb
#SBATCH --account=lappi
#SBATCH --partition=small
#SBATCH --time=3-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --output=logs/incl_Pb_%j.out
#SBATCH --error=logs/incl_Pb_%j.err

# Usage (build/bin/dipole must already exist -- this script does NOT build):
#   sbatch run_puhti_Pb.sh
#   sbatch --export=ALL,OUTDIR=/scratch/otherproject/me run_puhti_Pb.sh
#
# Puhti wrapper around run_many_Pb.sh: sweeps pD0 against every Glauber-
# sampled dipole file in data/Pb/mve/ (17 impact parameters b), i.e.
# ~17 x 60 = ~1020 dipole runs at ~14s/run single-threaded (measured
# locally for the proton case; Pb dipole files are a similar size so
# comparable per-run cost is expected) -- roughly 4 CPU-hours total. At 40
# cores that's well under 10 minutes wall time; the 3-day time limit is
# unused headroom, not an estimate.
#
# Output (one file per rapidity, "b  pD0  dsigma_dy" columns, meant to be
# Simpson-integrated over b in Python afterwards) goes under $OUTDIR/files
# rather than the submission dir, since /scratch has no home-style quota
# (same rationale as run_puhti.sh); SLURM logs stay local under logs/.

set -e

module purge
module load gcc/11.3.0
module load gsl/2.7

mkdir -p logs
export OUTDIR=${OUTDIR:-/scratch/lappi/pagimeno/incl/Pb}
export CORES=${SLURM_CPUS_PER_TASK:-40}

./run_many_Pb.sh
