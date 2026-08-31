# Inclusive D0 photoproduction

This project computes the inclusive D0 photoproduction cross section in ultraperipheral collisions (UPCs) in the CGC framework.
The main observable is

$$
\frac{d\sigma}{dy\,d^2p_{D^0}}.
$$

The code supports different UPC channels and fragmentation functions, and it can be used for proton and nuclear targets.

Based on P. Gimeno-Estivill, T. Lappi, and H. Mäntysaari, *Inclusive D⁰ photoproduction in ultraperipheral collisions*, Phys. Rev. D 111, 114036 (2025) [[doi:10.1103/7741-585p](https://doi.org/10.1103/7741-585p)].

---

## Project structure

The core calculation lives in the C++ sources under [src](src). The main executable is built as

```bash
./build/bin/dipole
```

The repository also contains figure scripts such as:

- [cms_comparison.py](cms_comparison.py)
- [cross_section.py](cross_section.py)
- [fragmentation_comparison.py](fragmentation_comparison.py)
- [HymnD_pt_spectrum.py](HymnD_pt_spectrum.py)

The heavier batch/cluster workflows are kept in [local_workflows](local_workflows). They are useful for large scans, but they are not required for the basic calculation.

---

## Build

```bash
mkdir -p build
cd build
cmake ..
make -j4
```

Requirements:

- CMake
- GSL (GNU Scientific Library)

---

## Basic run

```bash
./build/bin/dipole <pD0> [<dipole_file>] <y>
```

Example:

```bash
./build/bin/dipole 3 data/proton/mve.dat 0
```

This prints the rapidity and differential cross section at the selected $p_{D^0}$ value.

The arguments are:

- `<pD0>`: transverse momentum of the D0 meson in GeV
- `<y>`: rapidity
- `[<dipole_file>]`: optional dipole file; if omitted, it can be read from `DIPOLE_FILE`

The simple example runner is [run_local.sh](run_local.sh).

---

## Inputs and targets

The code expects dipole input files of the form

- proton: `data/proton/mve.dat`
- nucleus: `data/Pb/...` or `data/Au/...` with Glauber-sampled files such as `glauber_mve_<b>`

The dipole input is taken from the MVE model used in the project. The code also supports different UPC channels and fragmentation functions.

---

## UPC channel and fragmentation functions

The code supports the following channels:

- `Xn0n`
- `An0n`
- `PL(AnAn)`

These are controlled through the environment variable `CHANNEL`.

The fragmentation function is selected through `FRAG_TYPE`:

- `BCFY`
- `KniehlKramer`
- `HymnD` (the public label for the LHAPDF-based fragmentation set used here)

### Uncertainty bands

The comparison plots in [cms_comparison.py](cms_comparison.py) include the dominant theory uncertainty for the HymnD curve:

- a factorization-scale band from varying $Q$ by a factor of $0.5$ and $2$ around the central scale.

The BK initial-condition uncertainty is only about $2\%$ in the relevant bins and is therefore not included in the displayed public band. The BK posterior samples are still available in the output directories under `data/Pb/bk_posterior/member_*/...`, and the replica set is also kept in `files/HymnD/member_*/...` for reference.

The final CMS-style theory band keeps the dominant factorization-scale variation, while the small BK contribution and the small HymnD replica contribution are omitted from the public comparison. The BK uncertainty is an independent variation of the dipole initial condition: the code samples many posterior configurations of the BK fit and evolves each one separately, while keeping the fragmentation function fixed. This probes the uncertainty in the underlying nuclear-gluon evolution rather than in the fragmentation-function fit itself.

---

## Differential cross section and normalization

The raw output from the large runs is per impact parameter $b$. To obtain a final $p_{D^0}$ spectrum, one still needs to integrate over $b$ (for example by Simpson's rule, weighted by $2\pi b$).

For a proton target, the impact-parameter integral is effectively absorbed into the proton normalization. The code includes the standard overall prefactor associated with the photon flux and the charm charge.

In short, the final observable is built from:

- photon flux,
- dipole amplitude,
- fragmentation function,
- impact-parameter integration,
- and the chosen UPC channel.

Units are GeV-based throughout.

---

## Notes 

For a clean public GitHub release, the main code and example workflow remain at the project root, while the heavy local production scripts are kept in [local_workflows](local_workflows). These are convenient for large scans and local batch production, but they are not required for a standard user to build and run the code.

---

## Relevant files

- [CMakeLists.txt](CMakeLists.txt)
- [src](src)
- [run_local.sh](run_local.sh)
- [cross_section.py](cross_section.py)
- [cms_comparison.py](cms_comparison.py)
- [fragmentation_comparison.py](fragmentation_comparison.py)
- [HymnD_pt_spectrum.py](HymnD_pt_spectrum.py)
- [local_workflows](local_workflows)