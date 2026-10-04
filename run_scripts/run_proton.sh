#!/bin/bash
# D0 cross section on a proton for every pD0 and y.
# TARGET=AA (default): proton with the Pb+Pb flux, in output/<CHANNEL>/proton_baseline/<FRAG>/.
# TARGET=pA: p+Pb at 8.16 TeV, in output/pPb/central_values/<FRAG>/.
# Usage: TARGET=pA FRAG_TYPE=HymnD ./run_scripts/run_proton.sh
# One file per y, with columns pD0  dsigma_dyd2pD0 (without the proton normalization).

# Default rapidities for the proton
: "${Y_VALS:=0.0 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0}"
# Defaults for p+Pb
if [[ "${TARGET:-AA}" == "pA" ]]; then
	: "${FLUX_MODEL:=WS}"
	: "${SIGMA_NN:=99}"
	export SIGMA_NN
fi
set -e
cd "$(dirname "$0")/.."
source run_scripts/config.sh

require_file() { [[ -f "$1" ]] || { echo "Error: $1 does not exist. $2" >&2; exit 1; }; }
# "-1.5" -> "-15": the y tag in file names
tag() { echo "$1" | tr -d '.'; }
# Wait until fewer than CORES jobs are running.
throttle() { while (( $(jobs -r | wc -l) >= CORES )); do sleep 0.2; done; }

DIPOLE_FILE=${DIPOLE_FILE:-data/proton/mve.dat}
if [[ "$TARGET" == "pA" ]]; then
	OUTDIR=${OUTDIR:-$OUTPUT_ROOT/pPb/central_values/$FRAG_TYPE}
	# no neutron class in p+Pb; the flux is in the name if it is not WS
	pA_flux_tag=$([[ "$FLUX_MODEL" == "WS" ]] && echo "" || echo "_${FLUX_MODEL}")
	outname="D0_pA_incl_${FRAG_TYPE}${pA_flux_tag}"
	description="p+Pb (TARGET=pA): lead emits, proton target, 8.16 TeV, Gamma_pA [sigma_NN=${SIGMA_NN} mb], no EMD"
else
	OUTDIR=${OUTDIR:-$OUTPUT_ROOT/$channel_tag/proton_baseline/$FRAG_TYPE}
	outname="D0_proton_baseline_incl_${FRAG_TYPE}_${channel_tag}${g1_tag}${flux_tag}"
	description="proton target with the Pb+Pb photon flux (TARGET=AA), ${channel_tag} channel"
fi

require_file "$DIPOLE_FILE"
mkdir -p "$OUTDIR"
tmpdir=$(mktemp -d)

echo "Running dipole (proton target, TARGET=$TARGET, $DIPOLE_FILE) over pD0 in {$PT_VALS}, y in {$Y_VALS} -> $OUTDIR"
echo "frag_type=$FRAG_TYPE channel=$CHANNEL flux_model=$FLUX_MODEL scale_factor=$SCALE_FACTOR"

for y in $Y_VALS; do
	for pt in $PT_VALS; do
		./build/bin/dipole "$pt" "$DIPOLE_FILE" "$y" > "${tmpdir}/pD0_${pt}_y$(tag "$y").dat" &
		throttle
	done
done
wait

for y in $Y_VALS; do
	ytag=$(tag "$y")
	{
		echo "# ${description}"
		echo "# proton target, ${FRAG_TYPE} fragmentation"
		echo "# generated      : $(date '+%Y-%m-%d %H:%M:%S %Z')"
		echo "# flux_model     : ${FLUX_MODEL}"
		echo "# scale_factor   : ${SCALE_FACTOR} (fragmentation scale Q = SCALE_FACTOR * m_T)"
		echo "# calls          : ${CALLS} (VEGAS calls per point)"
		echo "# dipole file    : ${DIPOLE_FILE}"
		echo "# fixed rapidity y : ${y}"
		echo "# pD0  dsigma_dydpt"
		for pt in $PT_VALS; do
			f="${tmpdir}/pD0_${pt}_y${ytag}.dat"
			val=$(awk '$1 !~ /^#/ {print $2}' "$f")
			[[ -n "$val" ]] && echo "$pt  $val"
		done
	} > "$OUTDIR/${outname}_y${ytag}.dat"
done

rm -rf "${tmpdir}"

echo "Done."
