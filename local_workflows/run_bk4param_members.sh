#!/bin/bash

# quick usage: ./run_bk4param_members.sh
#
# This is the proton-version of the member scan. Instead of the full Pb Glauber
# setup, we just loop over the BK 4-parameter member files and compute the pD0/y
# grid for each one. It's basically the same idea as the other member scripts,
# just without the extra b integration.
#
# The member files are already prepared in data/bk4param/mve/. We just read them
# and run the dipole code for each member. The output lands in
# $OUTBASE/member_<NNNN>/spectrum_y<Y>.dat.
#
# This is mostly for building the uncertainty band from the BK 4-parameter
# sample, and then cross_section.py does the normalization and combination.
#
#   MEMBERS="0 1 2" ./run_bk4param_members.sh   # quick subset test
#   ./run_bk4param_members.sh                   # full 100-member set

set -euo pipefail

BIN=./build/bin/dipole
DIPOLE_DIR=${DIPOLE_DIR:-data/bk4param/mve}
OUTBASE=${OUTBASE:-out/bk4param}
CORES=${CORES:-$(( $(nproc) / 2 ))}
PT_MIN=${PT_MIN:-0.1}
PT_STEP=${PT_STEP:-0.2}
PT_MAX=${PT_MAX:-12.0}
Y_VALS=${Y_VALS:-"0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0"}
CHANNEL=${CHANNEL:-An0n}
export CHANNEL
FRAG_TYPE=${FRAG_TYPE:-KniehlKramer}
export FRAG_TYPE
MEMBERS=${MEMBERS:-$(seq 0 99)}

mkdir -p "$OUTBASE"

for member in $MEMBERS; do
	member_tag=$(printf '%04d' "$member")
	dfile="$DIPOLE_DIR/member_${member_tag}.dat"
	if [[ ! -f "$dfile" ]]; then
		echo "skipping member $member_tag: $dfile not found"
		continue
	fi

	echo "=== member $member_tag ($(date)) ==="
	outdir="$OUTBASE/member_${member_tag}"
	mkdir -p "$outdir"
	tmpdir=$(mktemp -d)

	for y in $Y_VALS; do
		ytag=$(echo "$y" | tr -d '.')
		for pt in $(seq $PT_MIN $PT_STEP $PT_MAX); do
			DIPOLE_FILE="$dfile" "$BIN" "$pt" "$y" > "${tmpdir}/pD0_${pt}_y${ytag}.dat" &
			while (( $(jobs -r | wc -l) >= CORES )); do sleep 0.2; done
		done
	done
	wait

	for y in $Y_VALS; do
		ytag=$(echo "$y" | tr -d '.')
		outfile="$outdir/spectrum_y${ytag}.dat"
		{
			echo "# fixed rapidity y : ${y}"
			echo "# pD0  dsigma_dydpt"
			for pt in $(seq $PT_MIN $PT_STEP $PT_MAX); do
				f="${tmpdir}/pD0_${pt}_y${ytag}.dat"
				val=$(awk '$1 !~ /^#/ {print $2}' "$f")
				[[ -n "$val" ]] && echo "$pt  $val"
			done
		} > "$outfile"
	done

	rm -rf "${tmpdir}"
done

echo "Finished at $(date)"
