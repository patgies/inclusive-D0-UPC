#ifndef BCFY_GRID_HPP
#define BCFY_GRID_HPP

#include <memory>
#include "interpolation.hpp"

// Builds a z-interpolator for the DGLAP-evolved Braaten-Cheung-Fleming-Yuan
// (BCFY) c -> D0 fragmentation function at a fixed factorisation scale Q
// (GeV). The evolution was done in another repository with EKO package, this reads the
// resulting grid, inputs/bcfy_eko/bcfy_eko_0000.dat.
std::unique_ptr<Interpolator> MakeBCFYInterpolator(double Q);

#endif
