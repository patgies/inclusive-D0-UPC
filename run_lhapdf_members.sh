#!/bin/bash

# Usage: ./run_lhapdf_members.sh
#
# Runs run_many_Pb.sh once per LHAPDF replica member of the fragmentation
# function set (see LHAPDF_DIR/LHAPDF_SET.info: ErrorType is "replicas", not
# Hessian eigenvectors, so the uncertainty is the standard deviation across
# members, not a Hessian formula). Each member's output lands in its own
# directory: $OUTBASE/member_<NNNN>/files/D0_incl_LHAPDF_<channel>_Pb_y<Y>.dat
#
# cross_section.py picks these up automatically (via the OUTBASE/member_*
# glob) to build the LHAPDF error band; it falls back to the single
# already-existing files/D0_incl_LHAPDF_*.dat (member 0000 only, no errors)
# if no member directories exist yet.
#
# WARNING: this reruns the full pT/y/b grid once per member. With the
# default 101 members this is ~100x the cost of a single run_many_Pb.sh
# call -- easily hours, even with CORES parallelism. Use MEMBERS to try a
# small subset first, or submit on a cluster (see run_many_oberon.sh for the
# SLURM wrapper pattern).
#
#   MEMBERS="0 1 2 3 4" ./run_lhapdf_members.sh   # quick subset test
#   ./run_lhapdf_members.sh                       # full 101-member set

set -euo pipefail

LHAPDF_DIR=${LHAPDF_DIR:-data/prompt-D0-1-109}
LHAPDF_SET=${LHAPDF_SET:-prompt-D0-1-109}
OUTBASE=${OUTBASE:-out/lhapdf}
MEMBERS=${MEMBERS:-$(seq 0 101)}

mkdir -p "$OUTBASE"

for member in $MEMBERS; do
	member_tag=$(printf '%04d' "$member")
	member_file="$LHAPDF_DIR/${LHAPDF_SET}_${member_tag}.dat"
	if [[ ! -f "$member_file" ]]; then
		echo "skipping member $member_tag: $member_file not found"
		continue
	fi

	echo "=== member $member_tag ($(date)) ==="
	OUTDIR="$OUTBASE/member_${member_tag}" \
	FRAG_TYPE=LHAPDF \
	LHAPDF_FILE="$member_file" \
	bash run_many_Pb.sh
done

echo "Finished at $(date)"
