#!/bin/bash

# quick usage: ./run_bk_posterior_members.sh
#
# This is the same idea as run_HymnD_members.sh, but instead of varying the
# fragmentation function, we vary the dipole amplitude across BK posterior
# samples. So the fragmentation part stays fixed while the dipole input keeps
# changing.
#
# In other words: one member = one BK posterior sample, and we run the whole
# Pb grid for that sample. This is meant to estimate the BK initial-condition
# uncertainty. Not the most elegant setup, but it does the job.
#
# The member dirs live in data/Pb/bk_posterior/member_<NNNN>, and the output
# gets written in the same general way as the other runs. This one is mostly
# for the BK uncertainty band, not the central curve.
#
# Quick subset test:
#   MEMBERS="0 1 2" ./run_bk_posterior_members.sh
#   ./run_bk_posterior_members.sh

set -euo pipefail

BK_DIR=${BK_DIR:-data/Pb/bk_posterior}
OUTBASE=${OUTBASE:-out/bk_posterior}
MEMBERS=${MEMBERS:-$(seq 0 99)}
FRAG_TYPE=${FRAG_TYPE:-LHAPDF}
LHAPDF_FILE=${LHAPDF_FILE:-data/prompt-D0-1-109/prompt-D0-1-109_0000.dat}

mkdir -p "$OUTBASE"

for member in $MEMBERS; do
	member_tag=$(printf '%04d' "$member")
	member_dir="$BK_DIR/member_${member_tag}"
	if [[ ! -d "$member_dir" ]]; then
		echo "skipping member $member_tag: $member_dir not found"
		continue
	fi

	echo "=== member $member_tag ($(date)) ==="
	OUTDIR="$OUTBASE/member_${member_tag}" \
	DIPOLE_DIR="$member_dir" \
	DIPOLE_X0=0.01 \
	FRAG_TYPE="$FRAG_TYPE" \
	LHAPDF_FILE="$LHAPDF_FILE" \
	bash run_many_Pb.sh
done

echo "Finished at $(date)"
