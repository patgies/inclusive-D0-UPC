#!/bin/bash

# use: ./run_BCFY_scale_variation.sh
#
# Factorisation scale variation for the DGLAP evolved Braaten-Cheung-
# Fleming-Yuan (BCFY) c -> D0 fragmentation function, evaluated at
# Q = scale_factor * mT, mT = sqrt(mc^2 + pT^2) (see src/bcfy_grid.cpp).
#
# Output written in
# $OUTBASE/factor_<0.25|0.5|1.0|2.0|4.0>/files/D0_incl_BCFY_<channel>_Pb_y<Y>.dat
#
# cross_section.py / cms_comparison.py combine these to get the min/max
# scale band.
#
#  SCALE_FACTORS="0.25 1.0 4.0" ./run_BCFY_scale_variation.sh  # default

set -euo pipefail

OUTBASE=${OUTBASE:-out/BCFY_scale}
SCALE_FACTORS=${SCALE_FACTORS:-"0.25 1.0 4.0"}

mkdir -p "$OUTBASE"

for factor in $SCALE_FACTORS; do
	echo "=== scale factor $factor ($(date)) ==="
	OUTDIR="$OUTBASE/factor_${factor}" \
	FRAG_TYPE=BCFY \
	SCALE_FACTOR="$factor" \
	bash run_many_Pb.sh
done

echo "Finished at $(date)"
