#include "kk_grid.hpp"
#include "hymnd_grid.hpp"
#include <cmath>
#include <cstdlib>

namespace {
const char* DefaultGridPath() { return "inputs/kk_eko/kk_eko_0000.dat"; }
const int kCharmPid = 4;
}

// Non-evolved Kniehl & Kramer c/b -> D0 input fragmentation functions
// D(x) at their natural starting scales mc/mb.
double KKInitialConditionC(double x) {
  double N = 0.694, eps = 0.101;
  double den = (1-x)*(1-x) + eps*x;
  return N * x * (1-x)*(1-x) / (den*den);
}

double KKInitialConditionB(double x) {
  double N = 81.7, alfa = 1.81, beta = 4.95;
  return N * std::pow(x, alfa) * std::pow(1-x, beta);
}

std::unique_ptr<Interpolator> MakeKniehlKramerInterpolator(double Q)
{
    // KK_EKO_FILE overrides the default grid path.
    const char* path = std::getenv("KK_EKO_FILE");
    return MakeHymnDZInterpolator(path ? path : DefaultGridPath(), kCharmPid, Q);
}
