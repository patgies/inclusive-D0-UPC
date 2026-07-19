# Inclusive D0 photoproduction

Differential cross section in rapidity and transverse momentum for inclusive D⁰ photoproduction in UPCs in the CGC formalism.
Select a UPC event (`Xn0n`, `An0n`, `PL(AnAn)`) and a fragmentation function ( [`KniehlKramer`](https://arxiv.org/abs/hep-ph/0504058), [`BCFY`](https://arxiv.org/abs/hep-ph/9409316), `LHAPDF`).

Based on P.Gimeno-Estivill, T.Lappi, H.Mäntysaari, *Inclusive D⁰ photoproduction in ultraperipheral collisions*, Phys. Rev. D 111, 114036 (2025) [[doi:10.1103/7741-585p](https://doi.org/10.1103/7741-585p)]

***
## Build
```
mkdir build
cd build
cmake ..
make
```
Requires CMake + GSL (GNU Scientific Library).

## Run
```
./build/bin/dipole <pD0> [<dipole_file>] <y>
```
prints the rapidity and transverse momentum of the $D0$ meson: `y  dsigma_dyd^2pD0` .

`<pD0>` and `<y>` are required. `[<dipole_file>]` can be given as the second positional arg, or omitted and read from the `DIPOLE_FILE` env var instead. See e.g. `run_local.sh` for a loop over pD0: `CHANNEL` env var picks the UPC event: `Xn0n`, `An0n`, or `PL(AnAn)` (default `An0n`).


Dipole files: `data/proton/mve.dat` (single file) or `data/Pb(Au)/mve/glauber_mve_<b>` (Glauber-sampled, one file per nuclear impact parameter `b`). See e.g. `run_many_Pb.sh` .

For Pb(Au), each `run_many_Pb.sh` output row is per impact parameter `b`, so `b` still needs to be integrated over (Simpson's rule, weighted by `b`) to get a `pD0`-only spectrum. See e.g., `cross_section.py` which does this internally and turns the raw `b pD0 dsigma_dypD0` grid into `dsigma/dpD0 dy` directly. For a proton target there is no such `b` and the impact parameter integral results in the proton size.


MVe dipole from [rcbkdipole](https://github.com/hejajama/rcbkdipole).

Units: GeV^n throughout.

## Fragmentation function
`FRAG_TYPE` env var picks the c → D⁰ fragmentation function: `BCFY`, `KniehlKramer` (default), or `LHAPDF`.
```
FRAG_TYPE="LHAPDF" ./build/bin/dipole 2.0 ./data/proton/mve.dat 1.0
```