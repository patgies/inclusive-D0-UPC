# Inclusive D0 photoproduction


This code provides the differential cross section for inclusive D⁰ photoproduction in ultraperipheral collisions in the Color Glass Condensate formalism. 
The user can select the UPC event: {"Xn0n", "An0n", "AnAn"} and fragmentation function: {KniehlKramer, LHAPDF, BCFY}.

Based on P.Gimeno-Estivill, T.Lappi, H.Mäntysaari, [2503.16108](https://arxiv.org/abs/2503.16108) 


***
## Files

- `src/main.cpp`: sets all physical parameters (charm mass, energy, EMD channel, VEGAS integration box, fragmentation scheme...), loads the dipole amplitude and the EMD photon flux table, and prints `y  dsigma_dy` for the D0 transverse momentum `pD0` given on the command line, either for a single rapidity `y` (if given) or for the hardcoded list `{0.0, 0.5, ..., 4.0}` (if omitted).
- `src/inclusive_spectrum.cpp`: defines the 5-dimensional integrand (fragmentation variable `zh`, photon plus-momentum `qp`, impact parameter `b`, gluon angle `phi`, gluon transverse momentum `l`) combining the photon flux, the hard photon-gluon-splitting amplitude and the fragmentation function, and integrates it with the GSL Monte Carlo Vegas algorithm.
- `src/gamma_aa.hpp` / `src/gamma_aa.cpp`: equivalent photon flux from the nucleus (Weizsäcker-Williams flux times nuclear form factor), electromagnetic-dissociation probability for the `An0n`/`Xn0n`/`PL(AnAn)` channels, and the `Gamma_AA(b)` table loader/interpolator.
- `src/fragmentation.hpp` / `src/fragmentation.cpp`: charm-to-D0 fragmentation functions — BCFY (`Dc_to_D0`), Kniehl & Kramer (`D_kniehl_kramer`), and an LHAPDF-backed fragmentation function (`D_lhapdf`).
- `src/amplitudelib.hpp` / `src/amplitudelib.cpp`: reads a dipole amplitude data file and provides `MakeSkInterpolator`, which builds an interpolator of the momentum-space dipole/gluon TMD `S_k(l)` at fixed `xbj` for use inside the Vegas integrand.
- `src/def.hpp`: the `parameters` struct threaded through the integrand (kinematics, fragmentation settings, photon-flux settings, integration box, precomputed `S_k` interpolator).
- `src/interpolation.hpp` / `src/interpolation.cpp`, `src/fourier.h` / `src/fourier.c`, `src/datafile.hpp` / `src/datafile.cpp`, `src/tools.hpp` / `src/tools.cpp`: supporting interpolation, Hankel/Fourier transform and dipole-datafile-parsing utilities used by `amplitudelib`.

## Input
- The dipole amplitude file is passed as the second command-line argument (see **Running locally**); if omitted it defaults to `./data/proton/mve.dat`. Dipole amplitude data files for Pb and Au (Glauber-sampled per impact parameter) live under `data/Pb/mve`, `data/Pb/mv1`, `data/Au/mve`, `data/Au/mv1`; the proton dipole is under `data/proton`.

The dipole amplitude in momentum space, `S_k(l, xbj)`, is precomputed once per rapidity via `dipole.MakeSkInterpolator(xbj, lmax)` in `main.cpp`, at `xbj = min(0.01, (mt/ss) * exp(-y))`, where `ss` is the collision energy [GeV] and `mt` the transverse mass of the produced D⁰.

Dipole parametrization MVe from [https://github.com/hejajama/rcbkdipole](https://github.com/hejajama/rcbkdipole).

- The equivalent-photon flux table `data/Gamma_AA.dat` (`b  T(b)`) is loaded at startup and is required for the nuclear EMD/flux calculation.
- The user needs to specify the momentum `pD0` of the produced D⁰ meson: see **Running locally**.

## Output
The program prints `y  dsigma_dy` for the given `pD0`, either for a single rapidity (if `y` is given on the command line) or, by default, for each rapidity in the hardcoded list `{0.0, 0.5, ..., 4.0}` in `main.cpp`. The 5D Vegas integral already integrates over the fragmentation variable `zh`, the photon momentum `qp`, the photon-emission impact parameter `b` (within `[bmin, bmax]`) and the gluon transverse momentum/angle `(l, phi)`.

The `b` in the *dipole file name* (e.g. `glauber_mve_10`) is a separate, nuclear-geometry impact parameter: each file is one Glauber sample, and results from different files still need to be combined/integrated over that impact parameter — see `run_many_Pb.sh` (production sweep: all Pb dipole files, all rapidities, regrouped into `b  pD0  dsigma_dy` per rapidity) and the example notebook `b_integral.ipynb` (Simpson integration over that impact parameter).

`cross_section.ipynb` shows how to turn the resulting `pD0`-grid into `dsigma/dpD0 dy`.

Note: this code uses GeV^n units.

## Building
Requires
- CMake
- GSL (GNU Scientific Library)
- LHAPDF (`lhapdf-config` must be on `PATH`), with the fragmentation function set used in `main.cpp` (e.g. `prompt-D0-1-109`) installed

How to compile:
```
mkdir -p build
cd build
cmake ..
make
```

## Running locally

### Direct executable run
The built executable takes a momentum argument, a dipole amplitude file, and an optional rapidity:
```
./build/bin/dipole <pD0> <dipole_file> [y]
```
If `y` is omitted, it prints `y  dsigma_dydpt` for every rapidity in the hardcoded list `{0.0, 0.5, ..., 4.0}`; if given, only that one rapidity is computed.

Examples:
```
./build/bin/dipole 2.0 ./data/Pb/mve/glauber_mve_10
./build/bin/dipole 2.0 ./data/Pb/mve/glauber_mve_10 1.5
```

### Local batch wrapper
`run_local.sh` is a quick local test script: it fixes a single dipole file and rapidity (edit `DIPOLE_FILE`/`Y` at the top of the script) and sweeps `pD0` from 0.1 to 12.0 in steps of 0.2, writing `pD0  dsigma_dydpt` pairs to `out/spectrum_y<Y>.dat`:
```
./run_local.sh
```

## SLURM / cluster execution

`run_Pb.sh` is a SLURM wrapper that loads modules and calls `run_many_Pb.sh`.
It requires `build/bin/dipole` to exist and `data/Pb/mve/` to contain the dipole files.

`run_many_Pb.sh` is a local-style script that sweeps `pD0` over all impact parameter files and writes one output file per rapidity.
