#!/bin/bash

# Run ./build/bin/dipole at a single rapidity, sweeping pD0, for one dipole file.
# Writes one file: pD0  dsigma_dydptd

BIN=./build/bin/dipole
DIPOLE_FILE=./data/proton/mve.dat
Y=1.0
OUTDIR=./out
CHANNEL=${CHANNEL:-An0n}
export CHANNEL

mkdir -p "$OUTDIR"

outfile="$OUTDIR/spectrum_y${Y}.dat"
echo "# pD0  dsigma_dydpt (y=$Y)" > "$outfile"

for pt in $(seq 0.1 0.2 12.0); do
    echo "pD0=$pt"
    val=$("$BIN" "$pt" "$DIPOLE_FILE" "$Y" | awk '$1 !~ /^#/ {print $2}')
    echo "$pt  $val" >> "$outfile"
done
