#!/bin/bash

# Run ./build/bin/dipole at a single rapidity, loop in pD0, for one dipole file.
# Writes one file: pD0  dsigma_dyd^2pD

BIN=./build/bin/dipole
DIPOLE_FILE=${DIPOLE_FILE:-./data/proton/mve.dat}
export DIPOLE_FILE
Y=1.0
OUTDIR=./out
CHANNEL=${CHANNEL:-An0n}
export CHANNEL

mkdir -p "$OUTDIR"

outfile="$OUTDIR/spectrum_y${Y}.dat"
echo "# pD0  dsigma_dydpt (y=$Y)" > "$outfile"

for pt in $(seq 0.1 0.2 12.0); do
    echo "pD0=$pt"
    val=$("$BIN" "$pt" "$Y" | grep -v '^#' | tr -s ' ' | cut -d' ' -f2)
    echo "$pt  $val" >> "$outfile"
done
