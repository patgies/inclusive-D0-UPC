#!/bin/bash
# D0 cross section for every Glauber sample b_d, pD0 and y, in output/<CHANNEL>/central_values/<FRAG>/.
# Usage: CHANNEL=Xn0n FRAG_TYPE=BCFY ./run_scripts/run_nucleus.sh
# One file per y, with columns b_d  pD0  dsigma_dyd2pD0. The plotting scripts integrate over b_d.

set -e
cd "$(dirname "$0")/.."
source run_scripts/config.sh

require_dir()  { [[ -d "$1" ]] || { echo "Error: $1 does not exist. $2" >&2; exit 1; }; }
# "-1.5" -> "-15": the y tag in file names
tag() { echo "$1" | tr -d '.'; }
# Wait until fewer than CORES jobs are running.
throttle() { while (( $(jobs -r | wc -l) >= CORES )); do sleep 0.2; done; }

TARGET=AA   # always the Pb+Pb photon flux
DIPOLE_DIR=${DIPOLE_DIR:-data/Pb/mve}
OUTDIR=${OUTDIR:-$OUTPUT_ROOT/$channel_tag/central_values/$FRAG_TYPE}

require_dir "$DIPOLE_DIR" "(expected Glauber samples glauber_mve_<b_d>)"
mkdir -p "$OUTDIR"
tmpdir=$(mktemp -d)

echo "Running dipole over $DIPOLE_DIR, pD0 in {$PT_VALS}, y in {$Y_VALS} -> $OUTDIR"
echo "frag_type=$FRAG_TYPE channel=$CHANNEL flux_model=$FLUX_MODEL scale_factor=$SCALE_FACTOR"

for dfile in "$DIPOLE_DIR"/glauber_mve_*; do
	b_d=$(basename "$dfile" | sed 's/glauber_mve_//')
	for pt in $PT_VALS; do
		for y in $Y_VALS; do
			./build/bin/dipole "${pt}" "$dfile" "${y}" > "${tmpdir}/b${b_d}_pD0_${pt}_y$(tag "$y").dat" &
			throttle
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
elif grep -q 'fragmentation.*HymnD' "$first_tmp"; then
	frag_tag="HymnD"
elif grep -q 'fragmentation.*BCFY' "$first_tmp"; then
	frag_tag="BCFY"
else
	echo "Error: could not detect a known fragmentation tag in $first_tmp's header:" >&2
	cat "$first_tmp" >&2
	exit 1
fi

for y in $Y_VALS; do
	ytag=$(tag "$y")
	outfile="$OUTDIR/D0_incl_${frag_tag}_${channel_tag}${g1_tag}_Pb${flux_tag}_y${ytag}.dat"

	{
		echo "$header"
		echo "#   generated             : $(date '+%Y-%m-%d %H:%M:%S %Z')"
		echo "#   flux_model            : ${FLUX_MODEL}"
		echo "#   calls                 : ${CALLS} (VEGAS calls per point)"
		echo "#   dipole file dir       : ${DIPOLE_DIR}/glauber_mve_<b_d>"
		echo "#   pD0 loop             : ${PT_VALS}"
		echo "#   fixed rapidity y      : ${y}"
		echo "# ============================================================"
		echo "# b_d  pD0  dsigma_dyd2pD0"
		for dfile in "$DIPOLE_DIR"/glauber_mve_*; do
			b_d=$(basename "$dfile" | sed 's/glauber_mve_//')
			for pt in $PT_VALS; do
				f="${tmpdir}/b${b_d}_pD0_${pt}_y${ytag}.dat"
				val=$(awk '$1 !~ /^#/ {print $2}' "$f")
				[[ -n "$val" ]] && echo "$b_d  $pt  $val"
			done
		done | sort -n -k1,1 -k2,2
	} > "$outfile"
done

rm -rf "${tmpdir}"

echo "Done. Next: python3 plotting_scripts/cms_comparison.py to integrate over b and plot."
