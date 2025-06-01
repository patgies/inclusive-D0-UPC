# Inclusive D0 photoproduction


This code provides the differential cross section for inclusive D0 photoproduction in ultraperipheral collisions in the "Xn0n" event.


Based on P.Gimeno-Estivill, T.Lappi, H.Mäntysaari, [2503.16108](https://arxiv.org/abs/2503.16108) 

See [https://github.com/hejajama/rcbkdipole](https://github.com/hejajama/rcbkdipole) 
***
## Files

- function.cpp : defines the integrand (photon flux, 0n emission probability, photon splitting into quark pair)
- Integral_parton.cpp : computes the partonic cross section for charm photoproduction using MC las vegas.
- integral_D0.cpp : adds the partonic cross section, defines the fragmentation function and computes the differential cross section for D0 photoproduction
- main.cpp: defines the parameters and calls integral_D0.cpp


## Input
-  `AmplitudeLib dipole("./data/Pb/mve/glauber_mve_2");` Dipole amplitudes: describe the scattering of the charm quark pair with lead nucleus.
From [https://github.com/hejajama/rcbkdipole](https://github.com/hejajama/rcbkdipole). 
These dipole amplitudes are evolved using the Balitsky-Kovchegov equation from initial Bjorken x: xbj=0.01 to `xbj = (mt / ss) * exp(-(y))` where 
where `ss` is the collision energy, `mt` is the invariant mass and `y` the rapidity of the meson. 


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

`./build/bin/dipole ${kd}`

where kd is the momentum of the produced D0 meson. 


### Questions and comments
Please send an email to patricia.p.gimenoestivill@jyu.fi
