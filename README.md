# Inclusive D0 photoproduction

This project computes the inclusive D0 photoproduction cross section `dσ / (dy d²p_D0)` in ultraperipheral collisions (UPCs) in the CGC framework.

The code supports different UPC channels, photon fluxes and fragmentation functions, and it can be used for proton and nuclear targets.

Based on P. Gimeno-Estivill, T. Lappi, and H. Mäntysaari, *Inclusive D⁰ photoproduction in ultraperipheral collisions*, Phys. Rev. D 111, 114036 (2025) [[doi:10.1103/7741-585p](https://doi.org/10.1103/7741-585p)].

---

## Build

```bash
mkdir build
cd build
cmake ..
make
```

Requirements:

- CMake
- GSL (GNU Scientific Library)

---

## Basic run

The core calculation lives in the C++ sources under [src](src). The run scripts are in [run_scripts](run_scripts) and the plotting scripts in [plotting_scripts](plotting_scripts). The main executable is built as

```bash
./build/bin/dipole <pD0> [<dipole_file>] <y>
```

Example:

```bash
./build/bin/dipole 3 data/proton/mve.dat 0
```

This prints the rapidity and the differential cross section (without prefactors) at the selected `p_D0` value.

The arguments are:

- `<pD0>`: transverse momentum of the D0 meson in GeV
- `<y>`: rapidity
- `[<dipole_file>]`: optional dipole file; if omitted, it is read from `DIPOLE_FILE`

The run scripts loop over `p_D0`, `y` and the dipole files. Their settings are environment variables, with the defaults in [run_scripts/config.sh](run_scripts/config.sh):

```bash
CHANNEL=An0n FRAG_TYPE=BCFY ./run_scripts/run_nucleus.sh   # Pb+Pb
FRAG_TYPE=BCFY ./run_scripts/run_scale_variation.sh        # fragmentation scale 0.5, 1, 2 x mT
TARGET=pA FRAG_TYPE=BCFY ./run_scripts/run_proton.sh       # p+Pb
```

The results go to `output/`. The plots are made with `plotting_scripts/PbPb_bins.py` (Pb+Pb, with the CMS data) and `plotting_scripts/pPb_bins.py` (p+Pb).

---

## Inputs and targets

The code expects dipole input files of the form

- proton: `data/proton/mve.dat`
- nucleus: `data/Pb/mve/glauber_mve_<b_d>`, one Glauber-sampled file per dipole impact parameter `b_d`

Dipole parametrization MVe from [https://github.com/hejajama/rcbkdipole](https://github.com/hejajama/rcbkdipole).

The other inputs are in `input/`: the fragmentation function grids, the photon flux tables, the BK posterior samples and the CMS data.

---

## UPC channel, photon flux and fragmentation functions

The code supports the following channels, selected through `CHANNEL`:

- `Xn0n`
- `An0n`
- `PL(AnAn)`

The photon flux is selected through `FLUX_MODEL`:

- `EFF` (default): effective flux of K. J. Eskola, V. Guzey, I. Helenius, P. Paakkinen, and H. Paukkunen, "Spatial resolution of dijet photoproduction in near-encounter ultraperipheral nuclear collisions," Phys. Rev. C 110, 054906 (2024) [arXiv:2404.09731].
- `PL`: point-like flux.
- `WS`: Woods-Saxon flux.
- `TABLE`: a flux from a file, `FLUX_FILE=<file>`, with the columns `z_gamma  f(z_gamma)`.

`TARGET=pA` computes p+Pb collisions at 8.16 TeV instead of Pb+Pb at 5.36 TeV.

The fragmentation function is selected through `FRAG_TYPE`, at the scale `Q = SCALE_FACTOR * mT`:

- `BCFY`: E. Braaten, K.-m. Cheung, S. Fleming, and T.-C. Yuan, Phys. Rev. D 51, 4819–4829 (1995).
- `KniehlKramer`: B. A. Kniehl and G. Kramer, Phys. Rev. D 74 (2006) 037502 [arXiv:hep-ph/0607306].
- `HymnD`: Epele, Hekhorn, Helenius, Paukkunen, and Zurita, arXiv:2609.10327.

`BCFY` and `KniehlKramer` are DGLAP-evolved with [eko](https://github.com/NNPDF/eko) in the repository [EKO-FF](https://github.com/patgies/EKO-FF).

The details, and how to use another flux or fragmentation function, are in [notes_photon_flux.pdf](notes_photon_flux.pdf) and [notes_FF.pdf](notes_FF.pdf).

### Uncertainty bands

The bands in the plots are the fragmentation scale variation, `Q` between `0.5` and `2` times the central scale (`run_scale_variation.sh`).

The HymnD replicas and the BK initial-condition uncertainty (about 2% in the CMS bins) can be computed with `run_scripts/run_members.sh`, but they are not included in the plots.

---

## Differential cross section and normalization

The output of the nuclear runs is per dipole impact parameter `b_d`. The plotting scripts integrate over `b_d` (Simpson's rule, multiplied by `2πb_d`) and multiply by the prefactors.

For a proton target, the impact-parameter integral is replaced by the proton normalization `16.36` mb of the MVe dipole parametrization.

Units are GeV throughout.

---
