#include "bcfy_grid.hpp"
#include "hymnd_grid.hpp"
#include <cstdlib>

namespace {
const char* DefaultGridPath() { return "inputs/bcfy_eko/bcfy_eko_0000.dat"; }
const int kCharmPid = 4;
}

std::unique_ptr<Interpolator> MakeBCFYInterpolator(double Q)
{
    // BCFY_EKO_FILE overrides the default grid path.
    const char* path = std::getenv("BCFY_EKO_FILE");
    return MakeHymnDZInterpolator(path ? path : DefaultGridPath(), kCharmPid, Q);
}
