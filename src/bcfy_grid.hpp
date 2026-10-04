#ifndef BCFY_GRID_HPP
#define BCFY_GRID_HPP

#include <memory>
#include "interpolation.hpp"

// BCFY c -> D0 fragmentation function D(z_h) at the scale Q [GeV],
// from the DGLAP-evolved grid input/BCFY_EKO/bcfy_eko_0000.dat.
std::unique_ptr<Interpolator> MakeBCFYInterpolator(double Q);

#endif
