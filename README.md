# Inclusive D0 photoproduction


This code provides the differential cross section for inclusive D0 photoproduction in ultraperipheral collisions in the "Xn0n" event.

Based on P.Gimeno-Estivill, T.Lappi, H.Mäntysaari, [2503.16108](https://arxiv.org/abs/2503.16108) 


***
## Files

- function.cpp : defines the integrand (photon flux, electromagnetic dissociation projectile (0n), photon splitting+dipole scattering)
- Integral_parton.cpp : computes the partonic cross section for charm photoproduction using MonteCarlo method Vegas.
- integral_D0.cpp : adds the partonic cross section, defines the fragmentation function and computes the differential cross section for D0 photoproduction
- main.cpp: defines the parameters (charm mass, rapidity, collision energy...) and calls integral_D0.cpp


## Input
- User needs to specify the impact parameter X (GeV^(-1)) in `main.cpp`: `AmplitudeLib dipole("./data/Pb/mve/glauber_mve_X");`. 
The dipole amplitudes in momentum space `dipole -> S_k( l, par -> xbj, FUNDAMENTAL, 1.0);` in `function.cpp` are evolved for each impact parameter X independently using the Balitsky-Kovchegov equation from initial Bjorken x `xbj=0.01` to `xbj = (mt / ss) * exp(-(y))` where 
where `ss` is the collision energy (GeV), `mt` is the invariant mass (GeV) and `y` the rapidity of the meson.
Dipole parametrization from [https://github.com/hejajama/rcbkdipole](https://github.com/hejajama/rcbkdipole).

## Output

- D0 cross section differential in momentum, rapidity and impact parameter. The user needs to integrate over the 2D impact parameter: see example code
`b_integral.ipynb`

## Building
Requires
- CMake
- GSL

How to compile:
```
mkdir build
cd build
cmake ..
make
```
Execute:

`./build/bin/dipole <kd>`

Where `<kd>` is the momentum of the produced D⁰ meson, passed as a command-line argument.


### Questions and comments
Please send an email to patricia.p.gimenoestivill@jyu.fi
