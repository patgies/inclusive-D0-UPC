#!/bin/bash

# Usage: ./run_many_Pb.sh
#
# Loop in pD0 and y for every Glauber-sampled dipole file in data/Pb/mve/ (one file per impact
# parameter b, e.g. glauber_mve_0, glauber_mve_2, ..., glauber_mve_32), then
# combines the runs into one file per rapidity y, with columns:
#   b  pD0  dsigma_dyd2pD0
# File name: $OUTDIR/files/D0_incl_<frag>_<channel>_Pb_y<Y>.dat
#
# b is meant to be integrated over afterwards (e.g. Simpson's rule,
# weighted by b, this script only produces the raw per-b, per-pD0, per-y grid.

y_vals=${Y_VALS:-"0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0"}
OUTDIR=${OUTDIR:-.}
CORES=${CORES:-$(( $(nproc) / 2 ))}
PT_MIN=${PT_MIN:-0.1}
PT_STEP=${PT_STEP:-0.2}
PT_MAX=${PT_MAX:-12.0}
DIPOLE_DIR=${DIPOLE_DIR:-data/Pb/mve}
CHANNEL=${CHANNEL:-An0n}
export CHANNEL
FRAG_TYPE=${FRAG_TYPE:-KniehlKramer}
export FRAG_TYPE

mkdir -p "$OUTDIR/files"
tmpdir=$(mktemp -d)

for dfile in "$DIPOLE_DIR"/glauber_mve_*; do
	b=$(basename "$dfile" | sed 's/glauber_mve_//')
	for pt in $(seq $PT_MIN $PT_STEP $PT_MAX)
	do
		for y in $y_vals; do
			ytag=$(echo "$y" | tr -d '.')
			./build/bin/dipole "${pt}" "$dfile" "${y}" > "${tmpdir}/b${b}_pD0_${pt}_y${ytag}.dat" &
			while (( $(jobs -r | wc -l) >= CORES )); do sleep 0.2; done
		done
	done
done
wait


first_tmp=$(ls "${tmpdir}"/b*_pD0_*_y*.dat | head -1)

if [[ ! -s "$first_tmp" ]]; then
	echo "Error: $first_tmp is empty -- the dipole run for that (b, pD0, y) point likely crashed. Aborting instead of silently mislabeling the output." >&2
	exit 1
fi

header=$(grep '^#' "$first_tmp" | grep -v '^# y  dsigma_dy')

if grep -q 'fragmentation.*Kniehl & Kramer' "$first_tmp"; then
	frag_tag="KniehlKramer"
elif grep -q 'fragmentation.*LHAPDF' "$first_tmp"; then
	frag_tag="LHAPDF"
elif grep -q 'fragmentation.*BCFY' "$first_tmp"; then
	frag_tag="BCFY"
else
	echo "Error: could not detect a known fragmentation tag in $first_tmp's header:" >&2
	cat "$first_tmp" >&2
	exit 1
fi

channel_tag=$(grep '^# channel' "$first_tmp" | sed 's/.*: //' | tr -d '() ')
if [[ -z "$channel_tag" ]]; then
	echo "Error: could not detect a channel tag in $first_tmp's header." >&2
	exit 1
fi

for y in $y_vals; do
	ytag=$(echo "$y" | tr -d '.')
	outfile="$OUTDIR/files/D0_incl_${frag_tag}_${channel_tag}_Pb_y${ytag}.dat"

	{
		echo "$header"
		echo "#   dipole file dir       : ${DIPOLE_DIR}/glauber_mve_<b>"
		echo "#   pD0 loop             : ${PT_MIN} to ${PT_MAX} GeV, step ${PT_STEP}"
		echo "#   fixed rapidity y      : ${y}"
		echo "# ============================================================"
		echo "# b  pD0  dsigma_dyd2pD0"
		for dfile in "$DIPOLE_DIR"/glauber_mve_*; do
			b=$(basename "$dfile" | sed 's/glauber_mve_//')
			for pt in $(seq $PT_MIN $PT_STEP $PT_MAX); do
				f="${tmpdir}/b${b}_pD0_${pt}_y${ytag}.dat"
				val=$(awk '$1 !~ /^#/ {print $2}' "$f")
				[[ -n "$val" ]] && echo "$b  $pt  $val"
			done
		done | sort -n -k1,1 -k2,2
	} > "$outfile"
done

rm -rf "${tmpdir}"
