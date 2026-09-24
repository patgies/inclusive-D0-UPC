#ifndef KK_GRID_HPP
#define KK_GRID_HPP

#include <memory>
#include "interpolation.hpp"

// Builds a z-interpolator for the DGLAP-evolved Kniehl & Kramer c -> D0
// fragmentation function at a fixed factorisation scale Q (GeV). Charm
// channel only, evolved in another repository with EKO package from mu0=mc; this reads the
// resulting grid, inputs/kk_eko/kk_eko_0000.dat.
std::unique_ptr<Interpolator> MakeKniehlKramerInterpolator(double Q);

// Non-evolved c/b -> D0 input fragmentation functions D(x) at their
// own natural starting scales (mc, mb) 
double KKInitialConditionC(double x);
double KKInitialConditionB(double x);

#endif
