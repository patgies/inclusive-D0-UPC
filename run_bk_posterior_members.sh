#!/bin/bash

# Usage: ./run_bk_posterior_members.sh
#
# Nucleus-target analog of run_lhapdf_members.sh, but for the *other*
# uncertainty source: instead of varying the fragmentation function while
# holding the dipole amplitude fixed, this holds the fragmentation scheme
# fixed (FRAG_TYPE, default LHAPDF at its central member 0000,
# SCALE_FACTOR=1 -- matching the plotted central LHAPDF curve) and varies
# the dipole amplitude across 100 posterior samples of the BK
# initial-condition fit ("LOmvefit": Q_{s,0}^2, e_c, C^2, sigma0/2 -- see
# bk/posteriorsamples_100_LOmvefit.dat), each independently BK-evolved to
# the same rapidities as the central mve set. This isolates the
# BK-initial-condition uncertainty as a third uncertainty source (alongside
# the LHAPDF scale/replica ones) for the LHAPDF curve on the
# CMS-comparison plot.
#
# To get the same band for BCFY or Kniehl-Kramer instead, override
# FRAG_TYPE (LHAPDF_FILE is only read when FRAG_TYPE=LHAPDF, so it's simply
# ignored otherwise); output filenames are tagged by scheme
# (D0_incl_<FRAG_TYPE>_...), so this is safe to run into the same OUTBASE
# as the LHAPDF sweep without clobbering it:
#   FRAG_TYPE=BCFY ./run_bk_posterior_members.sh
#   FRAG_TYPE=KniehlKramer ./run_bk_posterior_members.sh
# Note: cms_comparison.py's compute_bk_posterior_theory_points() only reads
# the LHAPDF-tagged files today -- it would need generalizing to a
# frag_type argument to actually plot a BCFY/KK band from this output.
#
# data/Pb/bk_posterior/member_<NNNN>/glauber_mve_<b> are symlinks into
# bk/bks_Pbtargets_1/bks/<NNNN>/ic_208_<b>.dat (member index NNNN matches
# row NNNN+1 of bk/posteriorsamples_100_LOmvefit.dat, confirmed by
# cross-checking each file's own alphas_scaling against that row's C^2).
#
# DIPOLE_X0=0.01: the rcbk solver run that produced these files left x0 out
# of its header (see src/main.cpp's DIPOLE_X0 handling, and datafile.cpp's
# "Invalid x0!" check) -- 0.01 matches both the existing central mve
# dataset's header and the standard HERA-fit convention this family of
# Bayesian MVe fits uses.
#
#   MEMBERS="0 1 2" ./run_bk_posterior_members.sh   # quick subset test
#   ./run_bk_posterior_members.sh                   # full 100-member set

set -euo pipefail

BK_DIR=${BK_DIR:-data/Pb/bk_posterior}
OUTBASE=${OUTBASE:-out/bk_posterior}
MEMBERS=${MEMBERS:-$(seq 0 99)}
FRAG_TYPE=${FRAG_TYPE:-LHAPDF}
LHAPDF_FILE=${LHAPDF_FILE:-data/prompt-D0-1-109/prompt-D0-1-109_0000.dat}

mkdir -p "$OUTBASE"

for member in $MEMBERS; do
	member_tag=$(printf '%04d' "$member")
	member_dir="$BK_DIR/member_${member_tag}"
	if [[ ! -d "$member_dir" ]]; then
		echo "skipping member $member_tag: $member_dir not found"
		continue
	fi

	echo "=== member $member_tag ($(date)) ==="
	OUTDIR="$OUTBASE/member_${member_tag}" \
	DIPOLE_DIR="$member_dir" \
	DIPOLE_X0=0.01 \
	FRAG_TYPE="$FRAG_TYPE" \
	LHAPDF_FILE="$LHAPDF_FILE" \
	bash run_many_Pb.sh
done

echo "Finished at $(date)"
