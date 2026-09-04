#!/bin/bash

# use: ./run_KniehlKramer_scale_variation.sh
#
# Factorisation scale variation for the DGLAP
# evolved Kniehl & Kramer c -> D0 fragmentation function, evaluated at
# Q = scale_factor * mT, mT = sqrt(mc^2 + pT^2) (see src/kk_grid.cpp).
#
# Output written in
# $OUTBASE/factor_<0.5|1.0|2.0>/files/D0_incl_KniehlKramer_<channel>_Pb_y<Y>.dat
#
# cross_section.py / cms_comparison.py combine these to get the min/max
# scale band.
#
#  SCALE_FACTORS="0.5 1.0 2.0" ./run_KniehlKramer_scale_variation.sh  # default


set -euo pipefail

OUTBASE=${OUTBASE:-out/KniehlKramer_scale}
SCALE_FACTORS=${SCALE_FACTORS:-"0.5 1.0 2.0"}

mkdir -p "$OUTBASE"

for factor in $SCALE_FACTORS; do
	echo "=== scale factor $factor ($(date)) ==="
	OUTDIR="$OUTBASE/factor_${factor}" \
	FRAG_TYPE=KniehlKramer \
	SCALE_FACTOR="$factor" \
	bash run_many_Pb.sh
done

echo "Finished at $(date)"
