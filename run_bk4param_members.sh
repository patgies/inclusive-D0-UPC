#!/bin/bash

# Usage: ./run_bk4param_members.sh
#
# Proton-target analog of run_lhapdf_members.sh: runs the full pD0/y grid
# once per BK 4-parameter Bayesian-sampled dipole-amplitude member.
#
# data/theta_100_mve.dat holds 100 samples of (Q_{s,0}^2, e_c, C^2,
# sigma0/2) from Bayesian sampling of the BK initial condition; row i+1 of
# that file was BK-evolved (externally, via rcbkdipole/rbk) into
# data/bk4param/mve/member_<NNNN>.dat, NNNN = i zero-padded. This script
# does not run that evolution -- it only consumes the already-evolved
# member files.
#
# Unlike a nucleus target (run_many_Pb.sh), a proton target has no Glauber
# b-loop and no external b Simpson integral: see README.md, "For a proton
# target there is no such b and the impact parameter integral results in
# the proton size." So per member this just loops pD0 x y once, like
# run_local.sh, but parallelized across CORES the same way run_many_Pb.sh
# parallelizes its (b, pD0, y) loop.
#
# Each member's output: $OUTBASE/member_<NNNN>/spectrum_y<Y>.dat, columns
# "pD0  dsigma_dydpt" (raw dipole-binary units, not yet normalized).
# cross_section.py's load_bk4param_band() picks these up (via the
# OUTBASE/member_* glob), applies the (2*pi)*pT*factor_A*GEVSQR_TO_MB
# normalization, and combines them into a mean +/- std band, analogous to
# load_lhapdf_band().
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
