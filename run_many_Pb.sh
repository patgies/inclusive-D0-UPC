#!/bin/bash

# Usage: ./run_many_Pb.sh
#   OUTDIR=/scratch/me CORES=40 ./run_many_Pb.sh   # override output dir / concurrency
#
# Sweeps pD0 (dipole loops over all rapidities internally per run) against
# EVERY Glauber-sampled dipole file in data/Pb/mve/ (one file per impact
# parameter b, e.g. glauber_mve_0, glauber_mve_2, ..., glauber_mve_32), then
# splits the combined output into one file per rapidity y, with columns:
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
# dipole process per (pD0, b) pair.

y_vals="0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0"
OUTDIR=${OUTDIR:-.}
CORES=${CORES:-$(( $(nproc) / 2 ))}
PT_MIN=0.1
PT_STEP=0.2
PT_MAX=12.0
DIPOLE_DIR=data/Pb/mve

mkdir -p "$OUTDIR/files"
tmpdir=$(mktemp -d)

for dfile in "$DIPOLE_DIR"/glauber_mve_*; do
	b=$(basename "$dfile" | sed 's/glauber_mve_//')
	for pt in $(seq $PT_MIN $PT_STEP $PT_MAX)
	do
		./build/bin/dipole ${pt} "$dfile" > "${tmpdir}/b${b}_pD0_${pt}.dat" &
		# cap concurrency to $CORES
		while (( $(jobs -r | wc -l) >= CORES )); do wait -n; done
	done
done
wait

# Shared physics parameters are identical across all runs (only pD0 and the
# dipole file itself vary), so reuse one run's header for every output file.
first_tmp=$(ls "${tmpdir}"/b*_pD0_*.dat | head -1)
header=$(grep '^#' "$first_tmp" | grep -v '^#   dipole amplitude file' | grep -v '^#   pD0 ' | grep -v '^# y  dsigma_dy')

# Derive filename tags from the actual run instead of hardcoding them, so the
# name can't drift out of sync with param.use_kniehl_kramer / param.channel in main.cpp.
frag_tag=$(grep -q 'fragmentation.*Kniehl & Kramer' "$first_tmp" && echo "KniehlKramer" || echo "BCFY")
channel_tag=$(grep 'EMD channel' "$first_tmp" | sed 's/.*: //' | tr -d '() ')

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
		for f in "${tmpdir}"/b*_pD0_*.dat; do
			tag=$(basename "$f" .dat)                    # b<B>_pD0_<PT>
			b=$(echo "$tag"  | sed -E 's/^b([0-9.]+)_pD0_.*/\1/')
			pt=$(echo "$tag" | sed -E 's/^b[0-9.]+_pD0_//')
			val=$(awk -v y="$y" '$1 !~ /^#/ && ($1+0)==(y+0) {print $2}' "$f")
			[[ -n "$val" ]] && echo "$b  $pt  $val"
		done | sort -n -k1,1 -k2,2
	} > "$outfile"

	
done

rm -rf "${tmpdir}"
