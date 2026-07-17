# Inclusive D0 photoproduction

Differential cross section for inclusive D⁰ photoproduction in UPCs, in the CGC formalism.
Pick a UPC event (`Xn0n`, `An0n`, `AnAn`) and a fragmentation function (`KniehlKramer`, `BCFY`).

Based on P.Gimeno-Estivill, T.Lappi, H.Mäntysaari, [2503.16108](https://arxiv.org/abs/2503.16108)

***
## Build
```
mkdir -p build && cd build
cmake .. && make
```
Needs CMake + GSL.

## Run
```
./build/bin/dipole <pD0> <dipole_file> <y>
```
All three args required, no defaults/sweep. `CHANNEL` env var picks the UPC event (default `An0n`).

```
./build/bin/dipole 2.0 ./data/Pb/mve/glauber_mve_10 1.5
```
prints `y  dsigma_dyd^2p` (one line).

Dipole files: `data/proton/mve.dat` (single file) or `data/Pb(Au)/mve/glauber_mve_<b>` (Glauber-sampled, one file per nuclear impact parameter `b`). MVe dipole from [rcbkdipole](https://github.com/hejajama/rcbkdipole).

## Sweeping
No sweep built in — loop over `dipole` calls yourself, e.g. proton, fixed y:
```
for pt in $(seq 0.1 0.2 12.0); do
    val=$(./build/bin/dipole "$pt" ./data/proton/mve.dat 1.0 | awk '$1 !~ /^#/ {print $2}')
    echo "$pt  $val"
done > out/spectrum_y1.0.dat
```
For Pb/Au also loop over the `glauber_mve_*` files and tag each line with its `b`:
```
for y in 0.0 1.0 2.0; do
  for dfile in data/Pb/mve/glauber_mve_*; do
    b=$(basename "$dfile" | sed 's/glauber_mve_//')
    for pt in $(seq 0.1 0.5 12.0); do
        val=$(./build/bin/dipole "$pt" "$dfile" "$y" | awk '$1 !~ /^#/ {print $2}')
        echo "$b  $pt  $val"
    done
  done > out/D0_incl_Pb_y${y}.dat
done
```
That `b` (nuclear impact parameter, from the file name — not the photon-emission `b` already integrated inside the Vegas call) still needs integrating over: see `b_integral.ipynb` (Simpson, weighted by `b`). Proton has no such `b`. `cross_section.ipynb` turns the resulting `pD0` grid into `dsigma/dpD0 dy`.

Each `(pD0, dipole_file, y)` call is independent — parallelize with `xargs -P`/GNU `parallel` if it's slow.

Units: GeV^n throughout.
