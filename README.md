# Inclusive D0 photoproduction

Computes the inclusive D0 photoproduction cross section `dσ / (dy d²p_D0)` in ultraperipheral collisions (UPCs), in the CGC framework. It works for proton and nuclear targets, with several UPC channels, photon fluxes and fragmentation functions.

Based on P. Gimeno-Estivill, T. Lappi, and H. Mäntysaari, *Inclusive D⁰ photoproduction in ultraperipheral collisions*, Phys. Rev. D 111, 114036 (2025) [[doi:10.1103/7741-585p](https://doi.org/10.1103/7741-585p)].

---

## Build

```bash
mkdir build && cd build
cmake ..
make
```

Requirements: CMake and GSL (GNU Scientific Library). The plotting scripts need Python 3 with numpy, scipy and matplotlib.

## Repository layout

- [src](src): C++ sources
- [run_scripts](run_scripts): run scripts, with the default settings in [config.sh](run_scripts/config.sh)
- [plotting_scripts](plotting_scripts): plotting scripts
- `data/`: dipole files, `input/`: other inputs, `output/`: results, `plots/`: figures

## Basic run

```bash
./build/bin/dipole <pD0> [<dipole_file>] <y>
./build/bin/dipole 3 data/proton/mve.dat 0
```

- `<pD0>`: transverse momentum of the D0 meson, in GeV
- `<dipole_file>`: optional; if omitted, it is read from `DIPOLE_FILE`
- `<y>`: rapidity

This prints the rapidity and the differential cross section, without prefactors.

The run scripts loop over `pD0`, `y` and the dipole files. Their settings are environment variables, and the results go to `output/`:

```bash
CHANNEL=An0n FRAG_TYPE=BCFY ./run_scripts/run_nucleus.sh   # Pb+Pb
FRAG_TYPE=BCFY ./run_scripts/run_scale_variation.sh        # fragmentation scale 0.5, 1, 2 x mT
TARGET=AA FRAG_TYPE=BCFY ./run_scripts/run_proton.sh       # proton with the Pb+Pb flux (edge, R_pA)
TARGET=pA FRAG_TYPE=BCFY ./run_scripts/run_proton.sh       # p+Pb
```

The plots go to `plots/`:

```bash
python3 plotting_scripts/PbPb_bins.py              # Pb+Pb in bins of pT and y, with the CMS data
python3 plotting_scripts/PbPb_bins_subtracted.py   # Pb+Pb, inclusive (An0n) minus diffractive (0n0n)
python3 plotting_scripts/pPb_bins.py               # p+Pb
python3 plotting_scripts/RpA.py                    # R_pA, Eq. (40) of arXiv:2503.16108
python3 plotting_scripts/cms_comparison.py         # comparison with the CMS data, rapidity bins
```

`FLUX=Starlight` in front of the two Pb+Pb scripts plots the runs with the STARlight flux table (`output/An0n_Starlight/`).

## Inputs

Dipoles, MVe parametrization from [rcbkdipole](https://github.com/hejajama/rcbkdipole):

- proton: `data/proton/mve.dat`
- nucleus: `data/Pb/mve/glauber_mve_<b_d>`, one Glauber-sampled file per dipole impact parameter `b_d`

In `input/`:

- fragmentation function grids
- photon flux tables
- BK posterior samples
- CMS data
- `diffractive/`: the diffractive cross sections of [diffractive-D0-UPC](https://github.com/patgies/diffractive-D0-UPC), which `PbPb_bins_subtracted.py` subtracts

## Settings

| Variable | Values |
|---|---|
| `CHANNEL` | `An0n` (default), `Xn0n`, `PL(AnAn)` |
| `FLUX_MODEL` | `EFF` (default) [1]; `PL` (point-like); `WS` (Woods-Saxon); `TABLE`, a flux from `FLUX_FILE=<file>` with the columns `z_gamma  f(z_gamma)` |
| `TARGET` | `AA`: Pb+Pb at 5.36 TeV (default); `pA`: p+Pb at 8.16 TeV |
| `FRAG_TYPE` | `KniehlKramer` (default) [2], `BCFY` [3], `HymnD` [4] |
| `SCALE_FACTOR` | fragmentation scale `Q = SCALE_FACTOR * mT` (default 1) |
| `X_MODE` | `charm` (default) or `fixed`, see below |

**`X_MODE`.**
- `charm`: the dipole is evaluated at the `x` of the charm quark, `x = m_T e^{-y} / sqrt(s)` with `p_c = p_D0 / z_h`, as in arXiv:2503.16108.
- `fixed`: the `x` of `z_h = 1` is used for the whole `z_h` integral. It is about 5 times faster and meant for quick tests.

**Fragmentation functions.**
- `BCFY` and `KniehlKramer` are DGLAP-evolved with [eko](https://github.com/NNPDF/eko) in [EKO-FF](https://github.com/patgies/EKO-FF).
- The two parts of `BCFY` (direct D0, and D*0 → D0) are each normalized to one, then multiplied by their weights 0.168 and 0.39, as in Cacciari and Nason, JHEP 09 (2003) 006. So the integral of `D(z_h)` over `z_h` is 0.558 at the initial scale.

The details, and how to use another flux or fragmentation function, are in [notes_photon_flux.pdf](notes_photon_flux.pdf) and [notes_FF.pdf](notes_FF.pdf).

## Uncertainty bands

- The bands in the plots are the fragmentation scale variation, `Q` between `0.5` and `2` times the central scale (`run_scale_variation.sh`).
- The HymnD replicas and the BK initial-condition uncertainty (about 2% in the CMS bins) are computed with `run_scripts/run_members.sh`. `cms_comparison.py` shows them when their runs are in `output/An0n/HymnD_band/` and `output/An0n/bk_band/`.

## Normalization

- The nuclear runs give one result per `b_d`. The plotting scripts integrate over `b_d` (Simpson's rule, times `2πb_d`) and multiply by the prefactors.
- For a proton target, the `b_d` integral is replaced by the proton normalization `16.36` mb of MVe.
- The Glauber files stop at `b_d = 32 GeV⁻¹`, which leaves out the edge of the nucleus (about 20 nucleons). The plotting scripts add it as `A_edge` times the proton cross section with the same flux, fragmentation function and scale, as in arXiv:2503.16108. This needs `TARGET=AA ./run_scripts/run_proton.sh`, with results in `output/<CHANNEL>/proton_baseline/`. `EDGE=0` leaves the edge out.
- Units are GeV throughout.

## References

1. K. J. Eskola, V. Guzey, I. Helenius, P. Paakkinen, and H. Paukkunen, "Spatial resolution of dijet photoproduction in near-encounter ultraperipheral nuclear collisions," Phys. Rev. C 110, 054906 (2024) [arXiv:2404.09731].
2. B. A. Kniehl and G. Kramer, Phys. Rev. D 74 (2006) 037502 [arXiv:hep-ph/0607306].
3. E. Braaten, K.-m. Cheung, S. Fleming, and T.-C. Yuan, Phys. Rev. D 51, 4819–4829 (1995).
4. Epele, Hekhorn, Helenius, Paukkunen, and Zurita, arXiv:2609.10327.
