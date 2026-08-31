# Inclusive D0 photoproduction

This project computes the inclusive D0 photoproduction cross section `dσ / (dy d p_D0)` in ultraperipheral collisions (UPCs) in the CGC framework.

The code supports different UPC channels and fragmentation functions, and it can be used for proton and nuclear targets.

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

The core calculation lives in the C++ sources under [src](src). The analysis and plotting scripts are kept in [python](python). The main executable is built as

```bash
./build/bin/dipole <pD0> [<dipole_file>] <y>
```

Example:

```bash
./build/bin/dipole 3 data/proton/mve.dat 0
```

This prints the rapidity and differential cross section at the selected `p_D0` value.

The arguments are:

- `<pD0>`: transverse momentum of the D0 meson in GeV
- `<y>`: rapidity
- `[<dipole_file>]`: optional dipole file; if omitted, it can be read from `DIPOLE_FILE`

The simple example is [run_local.sh](run_local.sh).

---

## Inputs and targets

The code expects dipole input files of the form

- proton: `data/proton/mve.dat`
- nucleus: `data/Pb/...` or `data/Au/...` with Glauber-sampled files such as `glauber_mve_<b>`

Dipole parametrization MVe from [https://github.com/hejajama/rcbkdipole](https://github.com/hejajama/rcbkdipole).

---

## UPC channel and fragmentation functions

The code supports the following channels:

- `Xn0n`
- `An0n`
- `PL(AnAn)`

These are controlled through the environment variable `CHANNEL`.

The fragmentation function is selected through `FRAG_TYPE`:

- `BCFY`: E. Braaten, K.-m. Cheung, S. Fleming, and T.-C. Yuan, "Perturbative QCD fragmentation functions as a model for heavy quark fragmentation," Phys. Rev. D 51, 4819–4829 (1995).
- `KniehlKramer`: B. A. Kniehl and G. Kramer, "Charmed-hadron fragmentation functions from CERN LEP1 revisited," Phys. Rev. D 74 (2006) 037502 [arXiv:hep-ph/0607306].
- `HymnD` 

### Uncertainty bands

The comparison plots in [python/cms_comparison.py](python/cms_comparison.py) include the dominant theory uncertainty for the HymnD curve:

- a factorization-scale band from varying `Q` by a factor of `0.5` and `2` around the central scale.

The BK initial-condition uncertainty is only about $2\%$ in the relevant bins and is therefore not included in the displayed public band. The BK posterior samples are still available in the output directories under `data/Pb/bk_posterior/member_*/...`, and the replica set is also kept in `files/HymnD/member_*/...` for reference.

The final CMS-style theory band keeps the dominant factorization-scale variation, while the small BK contribution and the small HymnD replica contribution are omitted. The BK uncertainty is an independent variation of the dipole initial condition: the code samples many posterior configurations of the BK fit and evolves each one separately, while keeping the fragmentation function fixed. 

---

## Differential cross section and normalization

The raw output from the large runs is per impact parameter `b`. To obtain a final `p_D0` spectrum, one still needs to integrate over `b` (for example by Simpson's rule, weighted by `2π b`).

For a proton target, the impact-parameter integral is effectively absorbed into the proton normalization. The code includes the standard overall prefactor associated with the photon flux and the charm charge.

The final observable is built from:

- photon flux,
- dipole amplitude,
- fragmentation function,
- impact-parameter integration,
- and the chosen UPC channel.

Units are GeV-based throughout.

---