#include "gamma_aa.hpp"
#include "def.hpp"
#include <memory>
#include <cmath>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_bessel.h>

using namespace std;


// Computes the Bessel function K_nu(x). GSL calls this a "modified
// Bessel function of the second kind". We only ever call it with
// nu = 1, so this is really just "K_1(x)" with a safety check added.
static double besselK(int nu, double x)
{
    gsl_sf_result result;
    int status = gsl_sf_bessel_Kn_e(nu, x, &result);

    // If GSL could not compute the value (for example because x is
    // huge and the true answer would be an extremely tiny number),
    // just treat the result as zero instead of crashing.
    if (status != GSL_SUCCESS) {
        return 0.0;
    }

    return result.val;
}

// Small helper so the rest of this file can just write
// get_GammaAA(b) instead of gamma_aa()(b).
static double get_GammaAA(double b)
{
    GammaAA& spline = gamma_aa();
    double value = spline(b);
    return value;
}


// The number of photons per unit qp and per unit transverse area,
// coming from a nucleus with charge Z, evaluated at impact parameter b.
//
//   qp = the photon's "plus momentum" (related to its energy)
//   b  = distance from the photon's emission point to the nucleus
double flux_density(double qp, double b, void* p)
{
    parameters* par = (parameters*)p;

    // Step 1: turn the plus-momentum qp into an ordinary photon energy.
    double omega = qp / sqrt(2.0);

    // Step 2: express that energy as a fraction "z" of the beam energy.
    // par->ss is sqrt(s), the total collision energy, so ss/2 is
    // (approximately) the energy carried by one beam.
    double beam_energy = par->ss / 2.0;
    double z = omega / beam_energy;

    // Step 3: build the argument that goes inside the Bessel function.
    double eta = z * par->mn * b;

    // If eta is very large, K_1(eta) is essentially zero, and computing
    // it directly could underflow or take a long time. Just return 0.
    if (eta > 50.0) {
        return 0.0;
    }

    // Step 4: compute the overall (dimensionless) prefactor.
    double alpha = par->alpha;
    double Z     = par->Z;
    double pref  = (alpha * Z * Z) / (M_PI * M_PI);

    // Step 5: compute K_1(eta) and combine everything.
    double K1 = besselK(1, eta);

    double eta_over_b = eta / b;
    double flux = (pref / qp) * eta_over_b * eta_over_b * K1 * K1;

    return flux;
}


// The full photon flux used in the cross-section integral: the plain
// photon number density from flux_density(), multiplied by
//   - "gamma", the probability that the emitting nucleus survives
//     without breaking apart, and
//   - "emd_factor", the probability of seeing the requested neutron
//     multiplicity (An0n, Xn0n, ...) from electromagnetic dissociation.
double photon_flux(double b, double qp, void* p)
{
    parameters* par = (parameters*)p;
    string channel = par->channel;

    double flux = flux_density(qp, b, p);

    // Work out the nuclear-breakup survival probability "gamma".
    double gamma;
    if (channel == "PL(AnAn)") {
        double b_cutoff = 14.0 / 0.197327;
        if (b >= b_cutoff) {
            gamma = 1.0;
        } else {
            gamma = 0.0;
        }
    } else {
        if (b < 150.0) {
            gamma = get_GammaAA(b);
        } else {
            gamma = 1.0;
        }
    }

    // Work out the electromagnetic-dissociation (neutron tagging)
    // probability "emd_factor".
    double P_emd    = par->S / (b * b);
    double P_no_emd = std::exp(-P_emd);

    double emd_factor = 1.0;
    if (channel == "An0n") {
        emd_factor = P_no_emd;
    } else if (channel == "Xn0n") {
        emd_factor = P_no_emd * (1.0 - P_no_emd);
    }

    double result = 2.0 * M_PI * b * flux * gamma * emd_factor;
    return result;
}

// These two unique_ptrs replace BOTH of the old constructs:
//   - the function-local "static TASpline my_spline(...)" inside the old
//     inline load_data_and_initialize()
//   - the raw global "GammaAA* g_gamma" that used to live in integrand.cpp
//
// They live in this anonymous namespace, so they're only visible inside
// this file (same idea as "static" at file scope in C). The rest of the
// program reaches them only through load_data_and_initialize() and
// gamma_aa() below.
//
// Lifetime: both unique_ptrs are constructed empty when the program
// starts and stay alive until the program exits (or you explicitly
// reset them) - the same "lives for the rest of the program" lifetime
// the old static/global pair had. The difference is ownership is now
// explicit and there is no leaked "new" with no matching "delete":
// the unique_ptrs free their objects automatically at program exit.
namespace {
    std::unique_ptr<TASpline> ta_spline;
    std::unique_ptr<GammaAA>  gamma_instance;
}

void load_data_and_initialize(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file!");
    }

    std::vector<double> b_values;
    std::vector<double> T_values;
    double b, T;

    while (file >> b >> T) {
        b_values.push_back(b);
        T_values.push_back(T);
    }

    // make_unique constructs the TASpline/GammaAA on the heap and hands
    // ownership to ta_spline/gamma_instance in one step - equivalent to
    // "new TASpline(...)" but with no risk of forgetting to free it later.
    ta_spline = std::make_unique<TASpline>(b_values, T_values);
    gamma_instance = std::make_unique<GammaAA>(*ta_spline);
}

GammaAA& gamma_aa()
{
    if (!gamma_instance) {
        throw std::runtime_error(
            "gamma_aa() called before load_data_and_initialize()");
    }
    return *gamma_instance;
}
