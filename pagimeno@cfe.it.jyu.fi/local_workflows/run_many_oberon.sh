#!/bin/bash -l
#SBATCH -n 1
#SBATCH -t 0-02:00:0

# quick usage (build/bin/dipole must already exist -- this script does NOT build):
#   sbatch run_many_oberon.sh
#   OUTDIR=~/incl/Pb CORES=20 sbatch run_many_oberon.sh
#
# This is basically the Oberon wrapper for run_many_Pb.sh. Same idea as the
# other cluster wrappers, just with the Slurm header and a few environment
# variables set so it doesn't completely explode when we submit it.
#
# Output goes under $OUTDIR/files. It's the raw b/pD0/y grid, so later Python
# does the actual b-integration and plotting. Nothing fancy, just the data dump.

set -euo pipefail

# build/bin/dipole is dynamically linked against libstdc++ (see
# src/CMakeLists.txt), so the compute node needs the same toolchain module
# loaded as the login/build node -- otherwise its default /lib64/libstdc++.so.6
# is too old (missing e.g. GLIBCXX_3.4.32) and the binary fails to even start.



OUTDIR=${OUTDIR:-$PWD/out}
# Leave 2 cores free by default instead of claiming the whole box -- override
# with CORES=N, or let SLURM_CPUS_PER_TASK drive it if -c was passed to sbatch.
CORES=${CORES:-${SLURM_CPUS_PER_TASK:-$(( $(nproc) - 2 ))}}
PT_MIN=${PT_MIN:-0.1}
PT_STEP=${PT_STEP:-0.2}
PT_MAX=${PT_MAX:-12.0}
Y_VALS=${Y_VALS:-"0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0"}
DIPOLE_DIR=${DIPOLE_DIR:-data/Pb/mve}

mkdir -p logs "$OUTDIR"
echo "Starting on $(hostname) at $(date)"
echo "OUTDIR=$OUTDIR CORES=$CORES PT=[$PT_MIN..$PT_MAX] step=$PT_STEP Y_VALS=$Y_VALS DIPOLE_DIR=$DIPOLE_DIR"

OUTDIR="$OUTDIR" \
CORES="$CORES" \
PT_MIN="$PT_MIN" PT_STEP="$PT_STEP" PT_MAX="$PT_MAX" \
Y_VALS="$Y_VALS" \
DIPOLE_DIR="$DIPOLE_DIR" \
bash run_many_Pb.sh

echo "Finished at $(date)"
