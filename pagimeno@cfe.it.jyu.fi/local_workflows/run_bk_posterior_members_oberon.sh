#!/bin/bash -l
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -t 0-02:00:0
#SBATCH --array=0-99

# quick usage (build/bin/dipole must already exist -- this script does NOT build):
#   sbatch run_bk_posterior_members_oberon.sh
#   OUTBASE=~/incl/bk_posterior CORES=20 sbatch run_bk_posterior_members_oberon.sh
#
# This is basically the same pattern as the HymnD member scan, but now each job
# uses a different BK posterior sample. The fragmentation is held fixed while
# the dipole input changes from one member to the next.
#
# So the idea is: we vary the dipole amplitude and see how much the result moves.
# That is the BK uncertainty part. This is not the same as the fragmentation
# replica scan, but it is the same general "submit one job per member" setup.
#
# The -c 20 thing matters again because if Slurm gives a strange CPU count, the
# background jobs can end up serializing and the whole thing slows down a lot.
#
# To rerun a few members only:
#   sbatch --array=3,17,42 run_bk_posterior_members_oberon.sh
# To be nicer to the cluster:
#   sbatch --array=0-99%10 run_bk_posterior_members_oberon.sh
#
# Output per member goes under $OUTBASE/member_<NNNN>/files, same layout as the
# other member-based jobs.

set -euo pipefail

# See run_HymnD_members_oberon.sh's note: the compute node needs the same
# toolchain module as the login/build node, or build/bin/dipole's dynamic
# libstdc++ won't load.
module load GCCcore/13.3.0

BK_DIR=${BK_DIR:-bk/bk_posterior}
OUTBASE=${OUTBASE:-$PWD/out/bk_posterior}
CORES=${CORES:-${SLURM_CPUS_PER_TASK:-$(( $(nproc) - 2 ))}}
PT_MIN=${PT_MIN:-0.1}
PT_STEP=${PT_STEP:-0.2}
PT_MAX=${PT_MAX:-12.0}
Y_VALS=${Y_VALS:-"-2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0"}
FRAG_TYPE=${FRAG_TYPE:-LHAPDF}
LHAPDF_FILE=${LHAPDF_FILE:-inputs/prompt-D0-1-109/prompt-D0-1-109_0000.dat}

member_tag=$(printf '%04d' "$SLURM_ARRAY_TASK_ID")
member_dir="$BK_DIR/member_${member_tag}"
if [[ ! -d "$member_dir" ]]; then
	echo "member $member_tag: $member_dir not found, skipping"
	exit 0
fi

mkdir -p logs "$OUTBASE"
echo "Starting member $member_tag on $(hostname) at $(date)"
echo "OUTBASE=$OUTBASE CORES=$CORES PT=[$PT_MIN..$PT_MAX] step=$PT_STEP Y_VALS=$Y_VALS"

OUTDIR="$OUTBASE/member_${member_tag}" \
CORES="$CORES" \
PT_MIN="$PT_MIN" PT_STEP="$PT_STEP" PT_MAX="$PT_MAX" \
Y_VALS="$Y_VALS" \
DIPOLE_DIR="$member_dir" \
DIPOLE_X0=0.01 \
FRAG_TYPE="$FRAG_TYPE" \
LHAPDF_FILE="$LHAPDF_FILE" \
bash run_many_Pb.sh

echo "Finished member $member_tag at $(date)"
