#!/bin/bash -l
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -t 0-02:00:0
#SBATCH --array=0-101

# Usage (build/bin/dipole must already exist -- this script does NOT build):
#   sbatch run_lhapdf_members_oberon.sh
#   OUTBASE=~/incl/lhapdf CORES=20 sbatch run_lhapdf_members_oberon.sh
#
# NOTE: -c 20 above requests 20 CPUs per array task so that
# SLURM_CPUS_PER_TASK is set correctly; without it, CORES falls back to
# "nproc - 2", which counts the whole node and not just what SLURM actually
# granted this task -- if that mismatches the real allocation (e.g. only 1
# CPU granted), the CORES background jobs serialize on however many CPUs
# were actually granted and the run silently takes far longer than
# expected (this bit us once: a run expected to take ~10min took ~1h).
# Adjust -c to match how many cores you actually want per array task.
#
# SLURM-array counterpart to run_lhapdf_members.sh: instead of looping over
# all 101 LHAPDF replica members serially in one job (101x a normal
# run_many_Pb.sh job -- easily a week+, well past any reasonable time limit),
# this submits one array task per member, each running a single member's
# full pT/y/b grid with the same ~2h budget as run_many_oberon.sh.
#
# To resubmit only missing/failed members, pass an explicit list, e.g.:
#   sbatch --array=3,17,42 run_lhapdf_members_oberon.sh
# To cap how many array tasks run concurrently (cluster-fairness), append
# "%N" to --array, e.g. --array=0-101%10 for at most 10 at once.
#
# To backfill extra rapidities into an already-completed set of members
# (run_many_Pb.sh only (re)writes the y-tagged files it's asked for, so this
# is safe and won't touch/duplicate already-finished y values):
#   Y_VALS="-2.0 -1.5 -1.0 -0.5" sbatch run_lhapdf_members_oberon.sh
#
# Output per member goes under $OUTBASE/member_<NNNN>/files, matching what
# cross_section.py's load_lhapdf_band() expects (via the OUTBASE/member_*
# glob) to build the mean +/- std-dev error band across replicas.

set -euo pipefail

LHAPDF_DIR=${LHAPDF_DIR:-data/prompt-D0-1-109}
LHAPDF_SET=${LHAPDF_SET:-prompt-D0-1-109}
OUTBASE=${OUTBASE:-$PWD/out/lhapdf}
CORES=${CORES:-${SLURM_CPUS_PER_TASK:-$(( $(nproc) - 2 ))}}
PT_MIN=${PT_MIN:-0.1}
PT_STEP=${PT_STEP:-0.2}
PT_MAX=${PT_MAX:-12.0}
Y_VALS=${Y_VALS:-"-2.0 -1.5 -1.0 -0.5 0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0"}
DIPOLE_DIR=${DIPOLE_DIR:-data/Pb/mve}

member_tag=$(printf '%04d' "$SLURM_ARRAY_TASK_ID")
member_file="$LHAPDF_DIR/${LHAPDF_SET}_${member_tag}.dat"
if [[ ! -f "$member_file" ]]; then
	echo "member $member_tag: $member_file not found, skipping"
	exit 0
fi

mkdir -p logs "$OUTBASE"
echo "Starting member $member_tag on $(hostname) at $(date)"
echo "OUTBASE=$OUTBASE CORES=$CORES PT=[$PT_MIN..$PT_MAX] step=$PT_STEP Y_VALS=$Y_VALS"

OUTDIR="$OUTBASE/member_${member_tag}" \
CORES="$CORES" \
PT_MIN="$PT_MIN" PT_STEP="$PT_STEP" PT_MAX="$PT_MAX" \
Y_VALS="$Y_VALS" \
DIPOLE_DIR="$DIPOLE_DIR" \
FRAG_TYPE=LHAPDF \
LHAPDF_FILE="$member_file" \
bash run_many_Pb.sh

echo "Finished member $member_tag at $(date)"
