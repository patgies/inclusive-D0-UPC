#!/bin/bash -l
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -t 0-02:00:0
#SBATCH --array=0-99

# quick usage (build/bin/dipole must already exist -- this script does NOT build):
#   sbatch run_bk4param_members_oberon.sh
#   OUTBASE=~/incl/bk4param sbatch run_bk4param_members_oberon.sh
#
# This is the Oberon version of the BK 4-parameter member scan. Same general
# idea as the others: one array task per member, each one runs the proton-case
# calculation for that specific dipole sample.
#
# The important part is again the CPU count. If Slurm gives something weird,
# the loop ends up fighting itself and the whole job takes much longer than it
# should. So we ask for a fixed number of cores and keep it consistent.
#
# To rerun only a few members:
#   sbatch --array=3,17,42 run_bk4param_members_oberon.sh
# To not swamp the cluster:
#   sbatch --array=0-99%10 run_bk4param_members_oberon.sh
#
# Output goes under $OUTBASE/member_<NNNN>, which is what the Python scripts
# expect when they combine the BK 4-parameter members into an uncertainty band.

set -euo pipefail

# build/bin/dipole is dynamically linked against libstdc++ (see
# src/CMakeLists.txt), so the compute node needs the same toolchain module
# loaded as the login/build node -- otherwise its default /lib64/libstdc++.so.6
# is too old (missing e.g. GLIBCXX_3.4.32) and the binary fails to even start.
module load GCCcore/13.3.0

DIPOLE_DIR=${DIPOLE_DIR:-bk/bk4param/mve}
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
