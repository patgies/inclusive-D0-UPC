#!/bin/bash -l
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -t 0-02:00:0
#SBATCH --array=0-99

# Usage (build/bin/dipole must already exist -- this script does NOT build):
#   sbatch run_bk_posterior_members_oberon.sh
#   OUTBASE=~/incl/bk_posterior CORES=20 sbatch run_bk_posterior_members_oberon.sh
#
# NOTE: -c 20 above requests 20 CPUs per array task so that
# SLURM_CPUS_PER_TASK is set correctly; without it, CORES falls back to
# "nproc - 2", which counts the whole node and not just what SLURM actually
# granted this task -- if that mismatches the real allocation, the CORES
# background jobs serialize on however many CPUs were actually granted and
# the run silently takes far longer than expected (see
# run_lhapdf_members_oberon.sh's note -- this bit us once there).
#
# SLURM-array counterpart to run_bk_posterior_members.sh: one array task per
# BK posterior-sample dipole member (100 of them, index 0-99), each running
# that member's full pT/y/b grid with the fragmentation scheme held fixed
# (FRAG_TYPE, default LHAPDF at its central member/scale -- SCALE_FACTOR
# defaults to 1.0 in main.cpp) -- i.e. the fragmentation function is fixed
# and only the dipole amplitude varies, the opposite of what
# run_lhapdf_members_oberon.sh does. This isolates the BK-initial-condition
# uncertainty as a third source (alongside the LHAPDF scale/replica ones)
# to combine in quadrature for the LHAPDF curve on the CMS-comparison plot.
#
# To get the same band for BCFY or Kniehl-Kramer instead, just override
# FRAG_TYPE (LHAPDF_FILE is only read when FRAG_TYPE=LHAPDF, so it's simply
# ignored otherwise); output filenames are tagged by scheme
# (D0_incl_<FRAG_TYPE>_...), so this is safe to run into the same OUTBASE
# as the LHAPDF sweep without clobbering it:
#   FRAG_TYPE=BCFY sbatch run_bk_posterior_members_oberon.sh
#   FRAG_TYPE=KniehlKramer sbatch run_bk_posterior_members_oberon.sh
# Note: cms_comparison.py's compute_bk_posterior_theory_points() only reads
# the LHAPDF-tagged files today -- it would need generalizing to a
# frag_type argument to actually plot a BCFY/KK band from this output.
#
# To resubmit only missing/failed members:
#   sbatch --array=3,17,42 run_bk_posterior_members_oberon.sh
# To cap concurrency (cluster-fairness):
#   sbatch --array=0-99%10 run_bk_posterior_members_oberon.sh
#
# Output per member goes under $OUTBASE/member_<NNNN>/files, matching the
# same layout run_lhapdf_members_oberon.sh uses (so the same kind of
# glob-and-percentile combination in cms_comparison.py can reuse the
# pattern for this band).

set -euo pipefail

# See run_lhapdf_members_oberon.sh's note: the compute node needs the same
# toolchain module as the login/build node, or build/bin/dipole's dynamic
# libstdc++ won't load.
module load GCCcore/13.3.0

BK_DIR=${BK_DIR:-data/Pb/bk_posterior}
OUTBASE=${OUTBASE:-$PWD/out/bk_posterior}
CORES=${CORES:-${SLURM_CPUS_PER_TASK:-$(( $(nproc) - 2 ))}}
PT_MIN=${PT_MIN:-0.1}
PT_STEP=${PT_STEP:-0.2}
PT_MAX=${PT_MAX:-12.0}
Y_VALS=${Y_VALS:-"-2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0"}
FRAG_TYPE=${FRAG_TYPE:-LHAPDF}
LHAPDF_FILE=${LHAPDF_FILE:-data/prompt-D0-1-109/prompt-D0-1-109_0000.dat}

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
