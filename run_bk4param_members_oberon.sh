#!/bin/bash -l
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -t 0-02:00:0
#SBATCH --array=0-99

# Usage (build/bin/dipole must already exist -- this script does NOT build):
#   sbatch run_bk4param_members_oberon.sh
#   OUTBASE=~/incl/bk4param sbatch run_bk4param_members_oberon.sh
#
# SLURM-array counterpart to run_bk4param_members.sh: one array task per BK
# 4-parameter member, each running that member's full pD0/y grid (no
# Glauber b-loop -- proton target, see run_bk4param_members.sh).
#
# NOTE: -c 20 above requests 20 CPUs per array task so that
# SLURM_CPUS_PER_TASK is set correctly; run_many_Pb.sh/run_bk4param_members.sh
# fall back to "nproc - 2" otherwise, which counts the whole node and not
# just what SLURM actually granted this task -- if that mismatches the real
# allocation, the CORES background jobs serialize on however many CPUs were
# actually granted and the run silently takes far longer than expected.
# Adjust -c to match how many cores you actually want per array task.
#
# To resubmit only missing/failed members, pass an explicit list, e.g.:
#   sbatch --array=3,17,42 run_bk4param_members_oberon.sh
# To cap how many array tasks run concurrently (cluster-fairness), append
# "%N" to --array, e.g. --array=0-99%10 for at most 10 at once.
#
# Output per member goes under $OUTBASE/member_<NNNN>, matching what
# cross_section.py's load_bk4param_band() expects (via the
# OUTBASE/member_* glob) to build the mean +/- std-dev error band across
# the 100 Bayesian-sampled BK initial-condition members.

set -euo pipefail

# build/bin/dipole is dynamically linked against libstdc++ (see
# src/CMakeLists.txt), so the compute node needs the same toolchain module
# loaded as the login/build node -- otherwise its default /lib64/libstdc++.so.6
# is too old (missing e.g. GLIBCXX_3.4.32) and the binary fails to even start.
module load GCCcore/13.3.0

DIPOLE_DIR=${DIPOLE_DIR:-data/bk4param/mve}
OUTBASE=${OUTBASE:-$PWD/out/bk4param}
CORES=${CORES:-${SLURM_CPUS_PER_TASK:-$(( $(nproc) - 2 ))}}
PT_MIN=${PT_MIN:-0.1}
PT_STEP=${PT_STEP:-0.2}
PT_MAX=${PT_MAX:-12.0}
Y_VALS=${Y_VALS:-"0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0"}

member_tag=$(printf '%04d' "$SLURM_ARRAY_TASK_ID")
member_file="$DIPOLE_DIR/member_${member_tag}.dat"
if [[ ! -f "$member_file" ]]; then
	echo "member $member_tag: $member_file not found, skipping"
	exit 0
fi

mkdir -p logs "$OUTBASE"
echo "Starting member $member_tag on $(hostname) at $(date)"
echo "OUTBASE=$OUTBASE CORES=$CORES PT=[$PT_MIN..$PT_MAX] step=$PT_STEP Y_VALS=$Y_VALS"

OUTBASE="$OUTBASE" \
DIPOLE_DIR="$DIPOLE_DIR" \
CORES="$CORES" \
PT_MIN="$PT_MIN" PT_STEP="$PT_STEP" PT_MAX="$PT_MAX" \
Y_VALS="$Y_VALS" \
MEMBERS="$SLURM_ARRAY_TASK_ID" \
bash run_bk4param_members.sh

echo "Finished member $member_tag at $(date)"
