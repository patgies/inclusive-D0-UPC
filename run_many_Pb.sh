#!/bin/bash

# Usage: ./run_many_Pb.sh
#   OUTDIR=/scratch/me CORES=40 ./run_many_Pb.sh   # override output dir / concurrency
#   CHANNEL="PL(AnAn)" ./run_many_Pb.sh             # override EMD channel (default An0n)
#
# Sweeps pD0 and y (one dipole process per (pD0, b, y) triple -- the binary
# takes a single required y argument, it does not loop internally) against
# EVERY Glauber-sampled dipole file in data/Pb/mve/ (one file per impact
# parameter b, e.g. glauber_mve_0, glauber_mve_2, ..., glauber_mve_32), then
# combines the runs into one file per rapidity y, with columns:
#   b  pD0  dsigma_dy
# File name: $OUTDIR/files/D0_incl_<frag>_<channel>_Pb_y<Y>.dat
# where <frag> (BCFY/KniehlKramer) and <channel> (e.g. An0n) are read back from
# the run's own header, matching whatever main.cpp was compiled with.
#
# b is meant to be integrated over afterwards in Python (Simpson's rule,
# weighted by b, i.e. the radial measure in 2D polar coordinates) -- this
# script only produces the raw per-b, per-pD0, per-y grid.
#
# Also called (with OUTDIR/CORES set) by run_csc_Pb.sh: build/bin/dipole
# has no internal OpenMP, so parallelism is at the process level, one
# dipole process per (pD0, b, y) triple.

y_vals=${Y_VALS:-"0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0"}
OUTDIR=${OUTDIR:-.}
CORES=${CORES:-$(( $(nproc) / 2 ))}
PT_MIN=${PT_MIN:-0.1}
PT_STEP=${PT_STEP:-0.2}
PT_MAX=${PT_MAX:-12.0}
DIPOLE_DIR=${DIPOLE_DIR:-data/Pb/mve}
CHANNEL=${CHANNEL:-An0n}
export CHANNEL

mkdir -p "$OUTDIR/files"
tmpdir=$(mktemp -d)

for dfile in "$DIPOLE_DIR"/glauber_mve_*; do
	b=$(basename "$dfile" | sed 's/glauber_mve_//')
	for pt in $(seq $PT_MIN $PT_STEP $PT_MAX)
	do
		for y in $y_vals; do
			ytag=$(echo "$y" | tr -d '.')
			./build/bin/dipole "${pt}" "$dfile" "${y}" > "${tmpdir}/b${b}_pD0_${pt}_y${ytag}.dat" &
			# cap concurrency to $CORES -- plain poll loop instead of `wait -n`
			# (bash >=4.3 only; some servers still ship an older bash as /bin/bash)
			while (( $(jobs -r | wc -l) >= CORES )); do sleep 0.2; done
		done
	done
done
wait

# Shared physics parameters are identical across all runs (only pD0, b and y
# vary), so reuse one run's header for every output file.
first_tmp=$(ls "${tmpdir}"/b*_pD0_*_y*.dat | head -1)
header=$(grep '^#' "$first_tmp" | grep -v '^# y  dsigma_dy')

# Derive filename tags from the actual run instead of hardcoding them, so the
# name can't drift out of sync with param.use_kniehl_kramer / param.channel in main.cpp.
frag_tag=$(grep -q 'fragmentation.*Kniehl & Kramer' "$first_tmp" && echo "KniehlKramer" || echo "BCFY")
channel_tag=$(grep '^# channel' "$first_tmp" | sed 's/.*: //' | tr -d '() ')

for y in $y_vals; do
	ytag=$(echo "$y" | tr -d '.')
	outfile="$OUTDIR/files/D0_incl_${frag_tag}_${channel_tag}_Pb_y${ytag}.dat"

	{
		echo "$header"
		echo "#   dipole file dir       : ${DIPOLE_DIR}/glauber_mve_<b>"
		echo "#   pD0 sweep             : ${PT_MIN} to ${PT_MAX} GeV, step ${PT_STEP}"
		echo "#   fixed rapidity y      : ${y}"
		echo "# ============================================================"
		echo "# b  pD0  dsigma_dy"
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
