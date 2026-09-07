#ifndef BCFY_GRID_HPP
#define BCFY_GRID_HPP

#include <memory>
#include "interpolation.hpp"

// Builds a z-interpolator for the DGLAP-evolved Braaten-Cheung-Fleming-Yuan
// (BCFY) c -> D0 fragmentation function at a fixed factorisation scale Q
// (GeV). Mirrors kk_grid.hpp's MakeKniehlKramerInterpolator.
//
// The pseudoscalar (c -> D0 direct) and vector (c -> D*0, feeding down to
// D0) channels are evolved as two independent DGLAP sets, both starting at
// mu0 = mc, and combined AFTER evolution (see bcfy_grid.cpp and
// testjobs/bcfyD0.cc in QCDnumFF for the underlying physics).
std::unique_ptr<Interpolator> MakeBCFYInterpolator(double Q);

#endif
