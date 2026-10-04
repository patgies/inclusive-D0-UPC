#!/bin/bash
# One run per member of a set, in output/<CHANNEL>/<set>_band/member_NNNN/.
# Usage: MEMBER_SET=HymnD MEMBERS="0 1 2" ./run_scripts/run_members.sh
#   HymnD    : HymnD replicas (input/HymnD/), Pb target
#   bk       : BK posterior samples (data/Pb/bk_posterior/), Pb target
#   bk4param : proton BK samples (input/BK/bk4param/mve/), proton target

# The bk set uses HymnD by default.
[[ "${MEMBER_SET:-HymnD}" == "bk" ]] && : "${FRAG_TYPE:=HymnD}"
set -e
cd "$(dirname "$0")/.."
source run_scripts/config.sh

require_file() { [[ -f "$1" ]] || { echo "Error: $1 does not exist. $2" >&2; exit 1; }; }
require_dir()  { [[ -d "$1" ]] || { echo "Error: $1 does not exist. $2" >&2; exit 1; }; }

MEMBER_SET=${MEMBER_SET:-HymnD}
HYMND_DIR=${HYMND_DIR:-input/HymnD}
HYMND_SET=${HYMND_SET:-prompt-D0-1-109}
BK_DIR=${BK_DIR:-data/Pb/bk_posterior}
BK4PARAM_DIR=${BK4PARAM_DIR:-input/BK/bk4param/mve}

case "$MEMBER_SET" in
	HymnD)
		OUTBASE=${OUTBASE:-$OUTPUT_ROOT/$channel_tag/HymnD_band}
		all_members=$(ls "$HYMND_DIR"/${HYMND_SET}_[0-9][0-9][0-9][0-9].dat 2>/dev/null | sed -E 's/.*_([0-9]{4})\.dat/\1/')
		;;
	bk)
		OUTBASE=${OUTBASE:-$OUTPUT_ROOT/$channel_tag/bk_band}
		all_members=$(ls -d "$BK_DIR"/member_[0-9][0-9][0-9][0-9] 2>/dev/null | sed -E 's/.*member_([0-9]{4})/\1/')
		;;
	bk4param)
		OUTBASE=${OUTBASE:-$OUTPUT_ROOT/$channel_tag/bk4param_band}
		all_members=$(ls "$BK4PARAM_DIR"/member_[0-9][0-9][0-9][0-9].dat 2>/dev/null | sed -E 's/.*member_([0-9]{4})\.dat/\1/')
		;;
	*)
		echo "Error: MEMBER_SET must be HymnD, bk or bk4param, not '$MEMBER_SET'." >&2
		exit 1
		;;
esac

MEMBERS=${MEMBERS:-$all_members}
if [[ -z "$MEMBERS" ]]; then
	echo "Error: no $MEMBER_SET members found (see the header of this script)." >&2
	exit 1
fi

for m in $MEMBERS; do
	member_tag=$(printf '%04d' "$((10#$m))")
	echo "=== $MEMBER_SET member $member_tag ($(date)) ==="
	case "$MEMBER_SET" in
		HymnD)
			member_file="$HYMND_DIR/${HYMND_SET}_${member_tag}.dat"
			require_file "$member_file"
			FRAG_TYPE=HymnD HYMND_FILE="$member_file" OUTDIR="$OUTBASE/member_${member_tag}" \
				./run_scripts/run_nucleus.sh
			;;
		bk)
			member_dir="$BK_DIR/member_${member_tag}"
			require_dir "$member_dir" "Run ./run_scripts/setup_bk_posterior_links.sh first."
			DIPOLE_DIR="$member_dir" DIPOLE_X0=0.01 OUTDIR="$OUTBASE/member_${member_tag}" \
				./run_scripts/run_nucleus.sh
			;;
		bk4param)
			member_file="$BK4PARAM_DIR/member_${member_tag}.dat"
			require_file "$member_file"
			DIPOLE_FILE="$member_file" OUTDIR="$OUTBASE/member_${member_tag}" \
				./run_scripts/run_proton.sh
			;;
	esac
done

echo "Finished at $(date)"
