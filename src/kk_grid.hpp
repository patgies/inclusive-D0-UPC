#ifndef KK_GRID_HPP
#define KK_GRID_HPP

#include <memory>
#include "interpolation.hpp"

// Builds a z-interpolator for the DGLAP-evolved Kniehl & Kramer c -> D0
// fragmentation function at a fixed factorisation scale Q (GeV).

// The charm-quark Peterson at mu0=mc and the bottom-quark
// at mu0=mb (Kniehl & Kramer hep-ph/0607306 LO)
// are evolved once per process and summed linearly for Q >= mb.
std::unique_ptr<Interpolator> MakeKniehlKramerInterpolator(double Q);

// Raw (non-evolved) c/b -> D0 input fragmentation functions, D(x), at their
// own natural starting scales (mc, mb) -- the un-evolved reference curves.
double KKInitialConditionC(double x);
double KKInitialConditionB(double x);

#endif
