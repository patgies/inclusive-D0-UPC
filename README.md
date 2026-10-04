# Inclusive D0 photoproduction

This project computes the inclusive D0 photoproduction cross section `dσ / (dy d²p_D0)` in ultraperipheral collisions (UPCs) in the CGC framework. The code supports different UPC channels, photon fluxes and fragmentation functions, for proton and nuclear targets.

Based on P. Gimeno-Estivill, T. Lappi, and H. Mäntysaari, *Inclusive D⁰ photoproduction in ultraperipheral collisions*, Phys. Rev. D 111, 114036 (2025) [[doi:10.1103/7741-585p](https://doi.org/10.1103/7741-585p)].

Two notes explain the details and how to change them:

- [notes_photon_flux.pdf](notes_photon_flux.pdf): the photon fluxes, the variables, and how to use another flux.
- [notes_FF.pdf](notes_FF.pdf): the fragmentation functions and how to add another one.

## Build

```bash
mkdir build
cd build
cmake ..
make
```

Requirements: CMake and GSL (GNU Scientific Library).

## Run

The C++ sources are in [src](src), the run scripts in [run_scripts](run_scripts) and the plotting scripts in [plotting_scripts](plotting_scripts). The executable is

```bash
./build/bin/dipole <pD0> [<dipole_file>] <y>
./build/bin/dipole 3 data/proton/mve.dat 0
```

with `pD0` the D0 transverse momentum in GeV and `y` its rapidity. It prints `y` and the cross section, without prefactors. If the dipole file is omitted, it is read from `DIPOLE_FILE`.

The run scripts loop over `pD0`, `y` and the dipole files, and can be started from any directory:

```bash
CHANNEL=Xn0n FRAG_TYPE=BCFY ./run_scripts/run_nucleus.sh       # Pb+Pb, one run per Glauber sample b_d
FRAG_TYPE=HymnD ./run_scripts/run_scale_variation.sh           # fragmentation scale 0.5, 1, 2 x mT
MEMBER_SET=HymnD ./run_scripts/run_members.sh                  # HymnD replicas (or MEMBER_SET=bk, bk4param)
./run_scripts/run_proton.sh                                    # proton with the Pb+Pb flux
TARGET=pA ./run_scripts/run_proton.sh                          # p+Pb
Y_VALS="0.0 1.0" PT_VALS="1.1 3.1" CALLS=1e4 ./run_scripts/run_nucleus.sh   # quick test
```

The cluster wrappers are in [run_scripts/roihu](run_scripts/roihu) and `run_scripts/oberon`.

## Settings

The settings are environment variables. The defaults are in [run_scripts/config.sh](run_scripts/config.sh).

| Variable | Default | Values |
|---|---|---|
| `CHANNEL` | `An0n` | neutron class: `An0n`, `Xn0n`, `PL(AnAn)` |
| `FLUX_MODEL` | `EFF` | photon flux: `EFF`, `PL`, `WS`, or `TABLE` with `FLUX_FILE=<file>` |
| `TARGET` | `AA` | `AA` (Pb+Pb, 5.36 TeV) or `pA` (p+Pb, 8.16 TeV) |
| `FRAG_TYPE` | `KniehlKramer` | fragmentation function: `BCFY`, `KniehlKramer`, `HymnD` |
| `SCALE_FACTOR` | `1.0` | fragmentation scale `Q = SCALE_FACTOR * mT` |
| `Y_VALS`, `PT_VALS` | see `config.sh` | rapidities and `pD0` values of the run scripts |
| `CALLS` | `2e5` | VEGAS calls per point |

- Photon flux: `EFF` is the effective flux of Eskola, Guzey, Helenius, Paakkinen and Paukkunen, Phys. Rev. C 110, 054906 (2024) [arXiv:2404.09731]. To use another flux, write it as a file with the columns `z_gamma  f(z_gamma)` and run with `FLUX_MODEL=TABLE FLUX_FILE=<file>`; no code has to be changed. See [notes_photon_flux.pdf](notes_photon_flux.pdf).
- Fragmentation functions: `BCFY` (Braaten, Cheung, Fleming and Yuan, Phys. Rev. D 51, 4819 (1995)) and `KniehlKramer` (Kniehl and Kramer, Phys. Rev. D 74, 037502 (2006)) are evolved with DGLAP using [eko](https://github.com/NNPDF/eko) in [EKO-FF](https://github.com/patgies/EKO-FF). `HymnD` is the set of Epele, Hekhorn, Helenius, Paukkunen and Zurita (arXiv:2609.10327). See [notes_FF.pdf](notes_FF.pdf).

## Folders

- `data/`: the dipole files. Proton: `data/proton/mve.dat`. Nucleus: `data/Pb/mve/glauber_mve_<b_d>`, one Glauber-sampled file per dipole impact parameter `b_d`. MVe parametrization from [rcbkdipole](https://github.com/hejajama/rcbkdipole).
- `input/BCFY_EKO/`, `input/KK_EKO/`, `input/HymnD/`: the fragmentation function grids.
- `input/WS_photon_flux/`: the survival factor `Gamma_AA(b)`, made with [src/make_gamma_aa.py](src/make_gamma_aa.py).
- `input/Starlight_photon_flux/`: the flux tables of P. Paakkinen, and the same in the `TABLE` format, made with [src/make_flux_table.py](src/make_flux_table.py).
- `input/flux_cache/`: the `EFF` flux, saved the first time it is computed. It can be deleted at any time.
- `input/BK/`: the posterior samples of the BK initial condition. `run_scripts/setup_bk_posterior_links.sh` links the Pb samples as `data/Pb/bk_posterior/member_NNNN/`.
- `input/CMS_data/`: the CMS D0 cross sections from HEPData.
- `output/<CHANNEL>/`: the results, in `central_values/<FRAG>/`, `scale_variation/<FRAG>/factor_<factor>/`, `HymnD_band/`, `bk_band/`, `bk4param_band/` and `proton_baseline/<FRAG>/`. The p+Pb results are in `output/pPb/`.

## Uncertainty bands

[plotting_scripts/cms_comparison.py](plotting_scripts/cms_comparison.py) compares the results with the CMS data. The band of the HymnD curve is the fragmentation scale variation, `Q` between `0.5` and `2` times the central scale (`run_scale_variation.sh`).

Two other uncertainties can be computed with `run_members.sh` but are not in the plot: the HymnD replicas (`MEMBER_SET=HymnD`) and the BK initial condition (`MEMBER_SET=bk` for Pb, `bk4param` for the proton), which is about 2% in the CMS bins.

## Differential cross section and normalization

The output of the nuclear runs is per dipole impact parameter `b_d`, in files with the columns `b_d  pD0  dsigma_dyd2pD0`. The plotting scripts integrate over `b_d` (Simpson's rule, multiplied by `2πb_d`) and multiply by the prefactors.

For a proton target, the impact-parameter integral is replaced by the proton normalization `16.36` mb of the MVe dipole parametrization.

Units are GeV throughout.
