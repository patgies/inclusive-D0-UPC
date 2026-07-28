#!/bin/bash

# Usage: ./run_lhapdf_scale_variation.sh
#
# Runs run_many_Pb.sh twice with the LHAPDF fragmentation scale
# Q = SCALE_FACTOR * mt0 set to the conventional up/down variation
# (0.5x and 2x) around the central scale (SCALE_FACTOR=1.0, already
# produced by the normal FRAG_TYPE=LHAPDF run into files/), using only the
# central LHAPDF member (member 0000 -- this is a scale-uncertainty
# envelope, not a PDF/FF-fit uncertainty, so it does not need to be redone
# per replica member).
#
# Output lands in $OUTBASE/factor_<0.5|2.0>/files/D0_incl_LHAPDF_<channel>_Pb_y<Y>.dat
#
# cross_section.py / cms_comparison.py combine these with the central
# (SCALE_FACTOR=1.0) run to build a min/max scale-variation envelope,
# separate from the replica (mean +/- std) band built by
# run_lhapdf_members.sh.
#
#   SCALE_FACTORS="0.5 2.0" ./run_lhapdf_scale_variation.sh   # default
#   SCALE_FACTORS="0.25 4.0" ./run_lhapdf_scale_variation.sh  # wider test

set -euo pipefail

LHAPDF_DIR=${LHAPDF_DIR:-data/prompt-D0-1-109}
LHAPDF_SET=${LHAPDF_SET:-prompt-D0-1-109}
OUTBASE=${OUTBASE:-out/lhapdf_scale}
SCALE_FACTORS=${SCALE_FACTORS:-"0.5 2.0"}

member_file="$LHAPDF_DIR/${LHAPDF_SET}_0000.dat"
if [[ ! -f "$member_file" ]]; then
	echo "Error: central member file not found: $member_file" >&2
	exit 1
fi

mkdir -p "$OUTBASE"

for factor in $SCALE_FACTORS; do
	echo "=== scale factor $factor ($(date)) ==="
	OUTDIR="$OUTBASE/factor_${factor}" \
	FRAG_TYPE=LHAPDF \
	LHAPDF_FILE="$member_file" \
	SCALE_FACTOR="$factor" \
	bash run_many_Pb.sh
done

echo "Finished at $(date)"
