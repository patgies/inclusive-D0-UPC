# Inclusive D0 photoproduction

Code to compute inclusive D0 photoproduction `dσ / (dy d²p_D0)` in ultraperipheral collisions in the CGC, for proton and nuclear targets.

Based on P. Gimeno-Estivill, T. Lappi, and H. Mäntysaari, *Inclusive D⁰ photoproduction in ultraperipheral collisions*, Phys. Rev. D 111, 114036 (2025) [[doi:10.1103/7741-585p](https://doi.org/10.1103/7741-585p)].

## Build

Needs CMake and GSL. The plotting scripts need Python 3 with numpy, scipy and matplotlib.

```bash
mkdir build && cd build
cmake ..
make
```

## Run

```bash
./build/bin/dipole <pD0> [<dipole_file>] <y>
./build/bin/dipole 3 data/proton/mve.dat 0
```

It prints the rapidity and the cross section (without prefactors). If no dipole file is given, it is read from `DIPOLE_FILE`.

The run scripts are in `run_scripts/`. The settings are environment variables, with the defaults in `run_scripts/config.sh`. Results go to `output/`.

```bash
CHANNEL=An0n FRAG_TYPE=BCFY ./run_scripts/run_nucleus.sh   # Pb+Pb
FRAG_TYPE=BCFY ./run_scripts/run_scale_variation.sh        # scale variation
TARGET=AA FRAG_TYPE=BCFY ./run_scripts/run_proton.sh       # proton with the Pb+Pb flux (edge, R_pA)
TARGET=pA FRAG_TYPE=BCFY ./run_scripts/run_proton.sh       # p+Pb
```

Plots (they go to `plots/`):

```bash
python3 plotting_scripts/PbPb_bins.py              # Pb+Pb with the CMS data
python3 plotting_scripts/PbPb_bins_subtracted.py   # inclusive (An0n) minus diffractive (0n0n)
python3 plotting_scripts/pPb_bins.py               # p+Pb
python3 plotting_scripts/RpA.py                    # R_pA, Eq. (40) of arXiv:2503.16108
python3 plotting_scripts/cms_comparison.py         # CMS comparison in rapidity bins
```

Use `FLUX=Starlight` with the Pb+Pb scripts to plot the runs with the STARlight flux (`output/An0n_Starlight/`).

## Inputs

Dipoles (MVe, from [rcbkdipole](https://github.com/hejajama/rcbkdipole)):

- proton: `data/proton/mve.dat`
- Pb: `data/Pb/mve/glauber_mve_<b_d>`, one file per impact parameter `b_d`

`input/` has the fragmentation functions, the flux tables, the BK posterior samples, the CMS data and the diffractive results from [diffractive-D0-UPC](https://github.com/patgies/diffractive-D0-UPC) (used in `PbPb_bins_subtracted.py`).

## Options

- `CHANNEL`: `An0n` (default), `Xn0n`, `PL(AnAn)`
- `FLUX_MODEL`: `EFF` (default, [arXiv:2404.09731](https://arxiv.org/abs/2404.09731)), `PL`, `WS`, or `TABLE` with `FLUX_FILE=<file>` (columns `z_gamma  f(z_gamma)`)
- `TARGET`: `AA` (Pb+Pb, 5.36 TeV, default) or `pA` (p+Pb, 8.16 TeV)
- `FRAG_TYPE`: `KniehlKramer` (default, [hep-ph/0607306](https://arxiv.org/abs/hep-ph/0607306)), `BCFY`, `HymnD` ([arXiv:2609.10327](https://arxiv.org/abs/2609.10327))
- `SCALE_FACTOR`: fragmentation scale `Q = SCALE_FACTOR * mT`
- `X_MODE`: `charm` (default) evaluates the dipole at the `x` of the charm quark, as in arXiv:2503.16108. `fixed` uses the `x` of `z_h = 1`, which is about 5 times faster, for quick tests.

`BCFY` and `KniehlKramer` are evolved with [eko](https://github.com/NNPDF/eko) in [EKO-FF](https://github.com/patgies/EKO-FF). The two parts of `BCFY` (direct D0 and D*0 → D0) are normalized with the weights 0.168 and 0.39 (Cacciari and Nason, JHEP 09 (2003) 006).

More details on the flux and the fragmentation functions are in [notes_photon_flux.pdf](notes_photon_flux.pdf) and [notes_FF.pdf](notes_FF.pdf).

## Uncertainties

- The bands in the plots are the fragmentation scale variation, `Q` from 0.5 to 2 times `mT`.
- The HymnD replicas and the BK initial condition (about 2%) are run with `run_members.sh`. They are only shown in `cms_comparison.py`.

## Normalization

The nuclear runs are per `b_d`, and the plotting scripts integrate over `b_d`. For the proton, the `b_d` integral is replaced by `16.36` mb. The Pb dipoles stop at `b_d = 32 GeV⁻¹`, so the edge of the nucleus (about 20 nucleons) is added as `A_edge` times the proton result, as in arXiv:2503.16108. This needs `TARGET=AA ./run_scripts/run_proton.sh`. `EDGE=0` turns it off.

Units are GeV.
