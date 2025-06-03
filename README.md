# Inclusive D0 photoproduction


This code provides the differential cross section for inclusive D⁰ photoproduction in ultraperipheral collisions in the "Xn0n" event.

Based on P.Gimeno-Estivill, T.Lappi, H.Mäntysaari, [2503.16108](https://arxiv.org/abs/2503.16108) 


***
## Files

- function.cpp : defines the integrand (photon flux, electromagnetic dissociation projectile (0n), photon splitting+dipole scattering).
- Integral_parton.cpp : computes the partonic cross section for charm photoproduction using MonteCarlo method Vegas.
- integral_D0.cpp : adds the partonic cross section, defines the fragmentation function and computes the differential cross section for D0 photoproduction.
- main.cpp: defines the parameters (charm mass, rapidity, collision energy...) and calls integral_D0.cpp.

## Input
- The user needs to specify the impact parameter X [GeV^-1] in `main.cpp`: `AmplitudeLib dipole("./data/Pb/mve/glauber_mve_X");`. 

The dipole amplitudes in momentum space `dipole -> S_k( l, par -> xbj, FUNDAMENTAL, 1.0);` in `function.cpp` are evolved for each impact parameter X independently using the BK equation from Bjorken  `xbj=0.01` to `xbj = (mt / ss) * exp(-(y))` where 
where `ss` is the collision energy [GeV], `mt` is the invariant mass [GeV] and `y` the rapidity of the meson.

Dipole parametrization MVe from [https://github.com/hejajama/rcbkdipole](https://github.com/hejajama/rcbkdipole).

- The user needs to specify the momentum of the produced D⁰ meson: see _Building_ section.

## Output
- Differential D0 cross section in momentum `kD`, rapidity `y` and impact parameter `b`: ```dsigma/d^2kD dy d^2b [GeV^(-2)]```.

The user needs to integrate over the impact parameter `b`: see example code `b_integral.ipynb`.

Finally, add the corresponding factors to obtain the differential cross section ```dsigma/dkD dy [GeV^(-3)]```: see example code `cross_section.ipynb`.

Note: this code uses GeV^n units.

## Building
Requires
- CMake
- GSL (GNU Scientific Library)

How to compile:
```
mkdir build
cd build
cmake ..
make
```
Execute:

`./build/bin/dipole <kd>`

Where `<kd>` is the momentum [GeV] of the produced D⁰ meson, passed as a command-line argument.


### Questions and comments
Please send an email to patricia.p.gimenoestivill@jyu.fi .
