#ifndef KK_GRID_HPP
#define KK_GRID_HPP

#include <memory>
#include "interpolation.hpp"

// Kniehl-Kramer c -> D0 fragmentation function D(z_h) at the scale Q [GeV],
// from the DGLAP-evolved grid input/KK_EKO/kk_eko_0000.dat.
std::unique_ptr<Interpolator> MakeKniehlKramerInterpolator(double Q);

#endif
