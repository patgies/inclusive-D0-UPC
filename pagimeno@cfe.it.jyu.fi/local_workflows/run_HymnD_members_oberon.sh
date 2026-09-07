#!/bin/bash -l
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -t 0-02:00:0
#SBATCH --array=0-101

# quick usage (build/bin/dipole must already exist -- this script does NOT build):
#   sbatch run_HymnD_members_oberon.sh
#   OUTBASE=~/incl/lhapdf CORES=20 sbatch run_HymnD_members_oberon.sh
#
# This is basically the cluster version of run_HymnD_members.sh. Instead of
# doing all 101 members in one job, we submit one array task per member so the
# whole thing does not sit on the cluster for a week.
#
# The important bit is the -c 20 line: without it, Slurm can give you a weird
# number of CPUs, and then the background jobs end up serializing on too many
# or too few cores and the run gets slower than it should. We learned that the
# hard way.
#
# To rerun only a few members:
#   sbatch --array=3,17,42 run_HymnD_members_oberon.sh
# To keep it a bit more polite to the cluster:
#   sbatch --array=0-101%10 run_HymnD_members_oberon.sh
#
# You can also add extra rapidities later if needed; run_many_Pb.sh just writes
# the ones you tell it to, so this is safe and won't overwrite the ones already
# done.
#
# Each member writes to $OUTBASE/member_<NNNN>/files, which is exactly the form
# that cross_section.py expects when it builds the replica band.

set -euo pipefail

# build/bin/dipole is dynamically linked against libstdc++ (see
# src/CMakeLists.txt), so the compute node needs the same toolchain module
# loaded as the login/build node -- otherwise its default /lib64/libstdc++.so.6
# is too old (missing e.g. GLIBCXX_3.4.32) and the binary fails to even start.
module load GCCcore/13.3.0

LHAPDF_DIR=${LHAPDF_DIR:-inputs/prompt-D0-1-109}
LHAPDF_SET=${LHAPDF_SET:-prompt-D0-1-109}
OUTBASE=${OUTBASE:-$PWD/out/HymnD}
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
