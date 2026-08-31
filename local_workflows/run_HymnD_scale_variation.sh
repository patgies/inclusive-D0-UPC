#!/bin/bash

# usage: ./run_HymnD_scale_variation.sh
#
# Standard 0.5x and 2x variations around
# the central scale, but only with the central HymnD member.
#
# So this is not a fit uncertainty; it's just a scale envelope. Output written in
# $OUTBASE/factor_<0.5|2.0>/files/D0_incl_HymnD_<channel>_Pb_y<Y>.dat
#
# Then cross_section.py / cms_comparison.py combine these with the central run
# to get the min/max scale band. This is separate from the replica spread from
# run_HymnD_members.sh.
#
#   SCALE_FACTORS="0.5 2.0" ./run_HymnD_scale_variation.sh   # default
#   SCALE_FACTORS="0.25 4.0" ./run_HymnD_scale_variation.sh  # wider test

set -euo pipefail

LHAPDF_DIR=${LHAPDF_DIR:-inputs/prompt-D0-1-109}
LHAPDF_SET=${LHAPDF_SET:-prompt-D0-1-109}
OUTBASE=${OUTBASE:-out/HymnD_scale}
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
