# Inclusive D0 photoproduction

Differential cross section in rapidity and transverse momentum for inclusive D⁰ photoproduction in UPCs in the CGC formalism.
Select a UPC event (`Xn0n`, `An0n`, `PL(AnAn)`) and a fragmentation function (`KniehlKramer` [[hep-ph/0504058](https://arxiv.org/abs/hep-ph/0504058)], `BCFY` [[hep-ph/9409316](https://arxiv.org/abs/hep-ph/9409316)], `LHAPDF`).

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

`<pD0>` and `<y>` are required. `dipole_file` can be given as the second positional arg, or omitted and read from the `DIPOLE_FILE` env var instead. See e.g. `run_local.sh` for a loop over pD0. `CHANNEL` env var picks the UPC event: `Xn0n`, `An0n`, or `PL(AnAn)` (default `An0n`).


Dipole files: `data/proton/mve.dat` (single file) or `data/Pb(Au)/mve/glauber_mve_<b>` (Glauber-sampled, one file per nuclear impact parameter `b`). See e.g. `run_many_Pb.sh` .

MVe dipole from [rcbkdipole](https://github.com/hejajama/rcbkdipole).

## Fragmentation function
`FRAG_TYPE` env var picks the c → D⁰ fragmentation function: `BCFY`, `KniehlKramer` (default), or `LHAPDF`.
```
FRAG_TYPE="LHAPDF" ./build/bin/dipole 2.0 ./data/proton/mve.dat 1.0
```
`LHAPDF` reads member 0 (central value) of an LHAPDF `lhagrid1`-format grid, evaluated at fixed `Q = m` (charm mass) — no dependency on the LHAPDF library itself. Defaults to `data/prompt-D0-1-109/prompt-D0-1-109_0000.dat`, override with `LHAPDF_FILE`. 



Units: GeV^n throughout.
