#include "gamma_aa.hpp"
#include "def.hpp"
#include <memory>
#include <cmath>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_bessel.h>

using namespace std;

static double besselK(int nu, double x)
{
    gsl_sf_result result;
    int status = gsl_sf_bessel_Kn_e(nu, x, &result);

    if (status != GSL_SUCCESS) {
        return 0.0;
    }

    return result.val;
}

static double get_GammaAA(double b)
{
    GammaAA& spline = gamma_aa();
    double value = spline(b);
    return value;
}

// Photon flux in the equivalent-photon approximation. The target spatial
// resolution is included through the impact-parameter dependence of the
// nuclear profile embedded in gamma_aa(b), following the discussion in
// K. J. Eskola, V. Guzey, I. Helenius, P. Paakkinen, and H. Paukkunen,
// "Spatial resolution of dijet photoproduction in near-encounter
// ultraperipheral nuclear collisions," Phys. Rev. C 110, 054906 (2024).
double flux_density(double qp, double b, void* p)
{
    parameters* par = (parameters*)p;

    double omega = qp / sqrt(2.0);
    double beam_energy = par->ss / 2.0;
    double z = omega / beam_energy;
    double eta = z * par->mn * b;

    if (eta > 50.0) {
        return 0.0;
    }

    double alpha = par->alpha;
    double Z     = par->Z;
    double pref  = (alpha * Z * Z) / (M_PI * M_PI);

    double K1 = besselK(1, eta);

    double eta_over_b = eta / b;
    double flux = (pref / qp) * eta_over_b * eta_over_b * K1 * K1;

    return flux;
}

double photon_flux(double b, double qp, void* p)
{
    parameters* par = (parameters*)p;
    string channel = par->channel;

    double flux = flux_density(qp, b, p);

    double gamma;
    if (channel == "PL(AnAn)") {
        double b_cutoff = 14.0 / 0.197327;
        if (b >= b_cutoff) {
            gamma = 1.0;
        } else {
            gamma = 0.0;
        }
    } else if (par->gamma_aa_one) {
        gamma = 1.0;
    } else {
        if (b < 150.0) {
            gamma = get_GammaAA(b);
        } else {
            gamma = 1.0;
        }
    }

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
