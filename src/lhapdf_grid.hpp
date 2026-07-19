#ifndef LHAPDF_GRID_HPP
#define LHAPDF_GRID_HPP

#include <memory>
#include <string>
#include "interpolation.hpp"

// Minimal reader for one member file of an LHAPDF "lhagrid1"-format grid
// (no dependency on the LHAPDF library itself). Builds a z-interpolator for
// one PDG parton flavor at a fixed scale Q, by linearly interpolating
// between the two Q-grid points (within whichever subgrid brackets Q) that
// straddle the requested Q.
//
// member_file : path to the grid file, e.g.
//               "data/prompt-D0-1-109/prompt-D0-1-109_0000.dat"
// pdg_flavor  : parton flavor id as listed in the grid's flavor line
//               (e.g. 4 for charm quark)
// Q           : fixed factorization scale at which to evaluate
//
// LHAPDF grids store x*f(x,Q) (the momentum-weighted convention used for
// both PDFs and fragmentation functions), so the returned interpolator's
// values are divided by z to give D(z,Q) directly, matching the convention
// of BCFY_DP/D_kniehl_kramer.
std::unique_ptr<Interpolator> MakeLHAPDFZInterpolator(
    const std::string& member_file, int pdg_flavor, double Q);

#endif
