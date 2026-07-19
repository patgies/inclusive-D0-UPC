# Inclusive D0 photoproduction

Differential cross section $d\sigma/dyd^2p_{D^0}$ for inclusive D⁰ photoproduction in UPCs in the CGC formalism.
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
prints the rapidity and differential cross section at that fixed $p_{D^0}$: $y \quad d\sigma/dyd^2p_{D^0}$.

`<pD0>` and `<y>` are required. `[<dipole_file>]` can be given as the second arg, or omitted and read from the `DIPOLE_FILE` env var instead. See e.g. `run_local.sh` for a loop over pD0, there `CHANNEL` env var picks the UPC event: `Xn0n`, `An0n`, or `PL(AnAn)` (default `An0n`).

- Dipole files: `data/proton/mve.dat` single file for proton target or `data/Pb(Au)/mve/glauber_mve_<b>` for nuclei (Glauber-sampled, one file per nuclear impact parameter `b`).
See e.g. `run_many_Pb.sh`,where each output row is per impact parameter $b$, so $b$ still needs to be integrated over (e.g. Simpson's rule, weighted by $2\pi b$) to get a $p_{D^0}$-only spectrum for each rapidity. For a proton target there is no such $b$ and the impact parameter integral results in the proton size. A factor $\alpha_{\textrm{em}}e^{2}_{c}N_{c}/(2\pi)^4$ needs to also be included.

 MVe dipole from [rcbkdipole](https://github.com/hejajama/rcbkdipole).



- The `FRAG_TYPE` env var picks the c → D⁰ fragmentation function: `BCFY`, `KniehlKramer` (default), or `LHAPDF`. Like `CHANNEL`, it's also set in batch scripts `run_local.sh` or `run_many_Pb.sh`.

Units: GeV^n throughout.