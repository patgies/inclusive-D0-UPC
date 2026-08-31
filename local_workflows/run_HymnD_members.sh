#!/bin/bash

# quick usage: ./run_HymnD_members.sh
#
# Runs run_many_Pb.sh once for each HymnD replica
# member. Each member gives a slightly different fragmentation function, and
# later we combine them to see how big the spread is.
#
# Output goes into a separate folder for each member:
# $OUTBASE/member_<NNNN>/files/D0_incl_HymnD_<channel>_Pb_y<Y>.dat
#
# cross_section.py picks these up automatically via the OUTBASE/member_*
# glob and builds the error band from the spread of the members. If there
# are no member dirs yet, it just falls back to the central file in files/.
#
#   MEMBERS="0 1 2 3 4" ./run_HymnD_members.sh   # quick subset test
#   ./run_HymnD_members.sh                       # full 101-member set

set -euo pipefail

LHAPDF_DIR=${LHAPDF_DIR:-inputs/prompt-D0-1-109}
LHAPDF_SET=${LHAPDF_SET:-prompt-D0-1-109}
OUTBASE=${OUTBASE:-out/HymnD}
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
