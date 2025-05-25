# Inclusive D0


This code provides the cross section for inclusive D0 photoproduction in ultraperipheral collisions. 

Based on [2503.16108](https://arxiv.org/abs/2503.16108) 

Dipoles amplitudes from [https://github.com/hejajama/rcbkdipole](https://github.com/hejajama/rcbkdipole) 
***


## Building
Requires

 - Cmake
 - GSL

How to compile:

```
   mkdir build
   cd build
   cmake ..
   make
```

This generates a library build/lib/libamplitude.a that you can link in your own program.

### Questions and comments
Please send an email to patricia.p.gimenoestivill@jyu.fi
