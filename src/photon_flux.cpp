#include "photon_flux.hpp"
#include "def.hpp"
#include <memory>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_bessel.h>

using namespace std;

// Every flux here is f = dN/dz_gamma, with z_gamma = omega/E_beam = sqrt(2) q+ / sqrt(s_NN),
// as in arXiv:2606.05469 Eqs. 15-22. The factor dz_gamma/dq+ = z_gamma/q+ is applied
// in integrand_inclusive() (cross_section_inclusive.cpp), not here.

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

// Point-like flux (arXiv:2606.05469 Eq. 18):
//   f(z_gamma,r) = alpha Z^2/(pi^2 z_gamma) * (z_gamma mN K1(z_gamma mN r))^2
double flux_density(double z_gamma, double r, void* p)
{
    parameters* par = (parameters*)p;

    double eta = z_gamma * par->mn * r;

    if (eta > 50.0) {
        return 0.0;
    }

    double alpha = par->alpha;
    double Z     = par->Z;
    double pref  = (alpha * Z * Z) / (M_PI * M_PI);

    double K1 = besselK(1, eta);

    double zmK = z_gamma * par->mn * K1;
    double flux = (pref / z_gamma) * zmK * zmK;

    return flux;
}

namespace {
    std::unique_ptr<TASpline>     ta_spline;
    std::unique_ptr<GammaAA>      gamma_instance;
    std::unique_ptr<WSFormFactor> ws_form_factor;
    std::unique_ptr<WSFluxTable>  ws_flux_table;
    double ws_RA = 0.0;       // GeV^-1
    double ws_r_switch = 0.0; // GeV^-1: WS flux below it, point-like flux above it
    const double WS_QMAX = 5.0;   // GeV
    uint64_t gamma_aa_hash = 0;   // identifies the Gamma_AA table in the flux cache
}

// Woods-Saxon flux (arXiv:2606.05469 Eq. 16):
//   f(z_gamma,r) = alpha Z^2/(pi^2 z_gamma) * I(z_gamma,r)^2
// Above 2R_PL = 14.2 fm it is the same as the point-like flux.
double flux_density_WS(double z_gamma, double r, void* p)
{
    parameters* par = (parameters*)p;
    if (!ws_flux_table)
        throw std::runtime_error("flux_density_WS() called before init_ws_form_factor()");
    if (r >= ws_r_switch) {
        return flux_density(z_gamma, r, p);
    }
    double pref = (par->alpha * par->Z * par->Z) / (M_PI * M_PI * z_gamma);
    return pref * ws_flux_table->I2(z_gamma, r);
}

namespace {
    std::unique_ptr<WSThickness> pA_TA;
    double pA_sigma_NN = 0.0;   // GeV^-2
}

void init_pA_flux(double sigma_NN_mb, double RA_fm, double a_fm, double B_mass)
{
    const double hbarc = 0.197327;
    double RA = RA_fm / hbarc, a = a_fm / hbarc;
    double s_max = RA + 12.0 * a;
    pA_TA.reset(new WSThickness(RA, a, B_mass, s_max));
    // 1 mb = 0.1 fm^2 = 0.1/hbarc^2 GeV^-2
    pA_sigma_NN = sigma_NN_mb * 0.1 / (hbarc * hbarc);
}

double GammaPA(double b)
{
    if (!pA_TA)
        throw std::runtime_error("GammaPA() called before init_pA_flux()");
    return std::exp(-pA_sigma_NN * (*pA_TA)(b));
}

// PL and WS: integrand of the b integral (arXiv:2606.05469 Eqs. 15, 20, 22),
//   2 pi b f(z_gamma,r=b) Gamma_AA(b) P_EMD(b)
double photon_flux(double b, double z_gamma, void* p)
{
    parameters* par = (parameters*)p;
    string channel = par->channel;

    // the flux is taken at the centre of the target nucleus, r = b
    double r     = b;
    double flux  = (par->flux_model == "WS") ? flux_density_WS(z_gamma, r, p)
                                              : flux_density(z_gamma, r, p);

    if (par->target == "pA") {
        // p+Pb (arXiv:2606.05469 Eq. 25): Gamma_pA(b) and no EMD factor. PL is zero below bmin.
        double gamma_pA = (par->flux_model == "PL") ? ((b >= par->bmin) ? 1.0 : 0.0)
                                                    : GammaPA(b);
        return 2.0 * M_PI * b * flux * gamma_pA;
    }

    double gamma;
    if (channel == "PL(AnAn)") {
        double b_cutoff = par->bmin;   // 2R_PL = 14.2 fm
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

    if (gamma <= 0.0) return 0.0;   // also avoids inf*0 for the point-like flux at b -> 0

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

void load_data_and_initialize(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file " + filename);
    }

    std::vector<double> b_values;
    std::vector<double> T_values;
    std::string line;

    // columns: b [GeV^-1]  Gamma_AA(b)
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        double b, T;
        if (iss >> b >> T) {
            b_values.push_back(b);
            T_values.push_back(T);
        }
    }

    // Hash of the table, so that the flux cache notices a new Gamma_AA table.
    gamma_aa_hash = 1469598103934665603ULL;
    for (size_t i = 0; i < b_values.size(); i++) {
        double v[2] = {b_values[i], T_values[i]};
        unsigned char bytes[sizeof(v)];
        std::memcpy(bytes, v, sizeof(v));
        for (size_t k = 0; k < sizeof(v); k++) { gamma_aa_hash ^= bytes[k]; gamma_aa_hash *= 1099511628211ULL; }
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

void init_ws_form_factor(double RA_fm, double a_fm, double mn)
{
    const double hbarc = 0.197327;
    ws_RA = RA_fm / hbarc;
    // Switch from WS to PL at 2R_PL = 14.2 fm.
    ws_r_switch = 14.2 / hbarc;
    double a_GeV = a_fm / hbarc;
    ws_form_factor.reset(new WSFormFactor(ws_RA, a_GeV, WS_QMAX));
    ws_flux_table.reset(new WSFluxTable(*ws_form_factor, mn, WS_QMAX,
                                         1e-6, 1.0, 160,           // z_gamma grid
                                         1.0, ws_r_switch, 260));  // r grid
}


// Fluxes already integrated over b: EFF (computed here) and TABLE (read from FLUX_FILE).
// Both are a spline of log f in log z_gamma and use the same file format:
//   column 1: z_gamma,  column 2: f(z_gamma) = dN/dz_gamma.  Lines with '#' are comments.
namespace {
    std::unique_ptr<WSThickness>   eff_TB;
    std::unique_ptr<EffFluxRadial> eff_H;
    gsl_spline* flux_spline = nullptr;
    gsl_interp_accel* flux_acc = nullptr;
    double flux_lz_gamma_min = 0.0, flux_lz_gamma_max = 0.0;
    bool flux_zero_above = false;        // TABLE: f = 0 above the last z_gamma
    std::string flux_source;             // name for the output header
    const double EFF_R_SWITCH = 800.0;   // GeV^-1
    const char* EFF_CACHE_DIR = "input/flux_cache";

    void set_flux_spline(const std::vector<double>& z_gamma, const std::vector<double>& f)
    {
        std::vector<double> lz_gamma(z_gamma.size()), lf(z_gamma.size());
        for (size_t i = 0; i < z_gamma.size(); i++) { lz_gamma[i] = std::log(z_gamma[i]); lf[i] = std::log(f[i]); }
        if (flux_spline) { gsl_spline_free(flux_spline); gsl_interp_accel_free(flux_acc); }
        flux_spline = gsl_spline_alloc(gsl_interp_cspline, z_gamma.size());
        flux_acc = gsl_interp_accel_alloc();
        gsl_spline_init(flux_spline, lz_gamma.data(), lf.data(), z_gamma.size());
        flux_lz_gamma_min = lz_gamma.front(); flux_lz_gamma_max = lz_gamma.back();
    }

    // Reads a flux table. key gets the "# key:" line of the file.
    bool read_flux_table(const std::string& filename, std::vector<double>& z_gamma, std::vector<double>& f, std::string* key)
    {
        std::ifstream file(filename);
        if (!file.is_open()) return false;
        std::string line;
        while (std::getline(file, line)) {
            if (line.rfind("# key:", 0) == 0 && key) *key = line.substr(6);
            if (line.empty() || line[0] == '#') continue;
            std::istringstream iss(line);
            double z_gamma_i, ff;
            if (iss >> z_gamma_i >> ff) { z_gamma.push_back(z_gamma_i); f.push_back(ff); }
        }
        return true;
    }

    bool write_flux_table(const std::string& filename, const std::vector<double>& z_gamma, const std::vector<double>& f,
                          const std::string& description, const std::string& key)
    {
        FILE* out = std::fopen(filename.c_str(), "w");
        if (!out) return false;
        std::fprintf(out, "# Photon flux f(z_gamma) = dN/dz_gamma, z_gamma = omega/E_beam = photon energy / beam energy per nucleon\n");
        std::fprintf(out, "# %s\n", description.c_str());
        std::fprintf(out, "# key:%s\n", key.c_str());
        std::fprintf(out, "# z_gamma  f(z_gamma)  z_gamma*f(z_gamma)\n");
        for (size_t i = 0; i < z_gamma.size(); i++) std::fprintf(out, "%.17g  %.17g  %.17g\n", z_gamma[i], f[i], z_gamma[i] * f[i]);
        return std::fclose(out) == 0;
    }
}

void init_effective_flux(const std::string& channel, void* p, double RA_fm, double a_fm, double B_mass)
{
    parameters* par = (parameters*)p;
    const double hbarc = 0.197327;
    const double ws_r_switch_fm = 14.2;
    int nz_gamma = 100;
    double z_gamma_min = 1e-6, z_gamma_max = 1.0;

    // The table does not depend on pD0, y or the dipole, so it is kept in input/flux_cache/.
    // The key lists what it depends on; if the key of the file is different, compute again.
    // FLUX_NO_CACHE=1 turns the cache off.
    std::string neutron_class = (channel == "An0n" || channel == "Xn0n") ? channel : "AnAn";
    char keybuf[512];
    std::snprintf(keybuf, sizeof(keybuf),
                  " v1 class=%s alpha=%.17g Z=%.17g mn=%.17g S=%.17g RA_fm=%.17g a_fm=%.17g B=%.17g"
                  " ws_switch_fm=%.17g eff_r_switch=%.17g nz=%d gamma_aa=%016llx",
                  neutron_class.c_str(), par->alpha, par->Z, par->mn, par->S, RA_fm, a_fm, B_mass,
                  ws_r_switch_fm, EFF_R_SWITCH, nz_gamma, (unsigned long long)gamma_aa_hash);
    std::string key = keybuf;
    std::string cache_file = std::string(EFF_CACHE_DIR) + "/EFF_" + neutron_class + ".dat";
    flux_source = "EFF (" + neutron_class + ", arXiv:2404.09731 Eq. 4)";
    flux_zero_above = false;

    {
        std::vector<double> z_gamma_cached, fc;
        std::string cached_key;
        if (!getenv("FLUX_NO_CACHE") && read_flux_table(cache_file, z_gamma_cached, fc, &cached_key)
            && cached_key == key && (int)z_gamma_cached.size() == nz_gamma) {
            set_flux_spline(z_gamma_cached, fc);
            return;
        }
    }

    init_ws_form_factor(RA_fm, a_fm, par->mn);

    double RA = RA_fm / hbarc, a = a_fm / hbarc;
    double s_max = RA + 12.0 * a;

    eff_TB.reset(new WSThickness(RA, a, B_mass, s_max));
    eff_H.reset(new EffFluxRadial(*eff_TB, channel, par->S, B_mass, EFF_R_SWITCH));

    // f_eff(z_gamma) = (2 pi/B) int_0^{r_max} r dr f_WS(z_gamma,r) H(r),  r_max = 60/(z_gamma*mn)
    double lz_gamma_min = std::log(z_gamma_min), lz_gamma_max = std::log(z_gamma_max);
    std::vector<double> z_gamma_grid(nz_gamma), logf(nz_gamma), fg(nz_gamma);
    for (int i = 0; i < nz_gamma; i++) {
        double lz_gamma_i = lz_gamma_min + (lz_gamma_max - lz_gamma_min) * i / (nz_gamma - 1.0);
        double z_gamma = std::exp(lz_gamma_i);
        z_gamma_grid[i] = z_gamma;
        double r_max = 60.0 / (z_gamma * par->mn);
        double r0 = 1e-3;
        int Nr = 2000;
        double aa = std::log(r0), cc = std::log(r_max), hh = (cc - aa) / Nr, sum = 0.0;
        for (int k = 0; k <= Nr; k++) {
            double r = std::exp(aa + k * hh);
            double flux = flux_density_WS(z_gamma, r, p);
            double Hval = (*eff_H)(r);
            double w = (k == 0 || k == Nr) ? 1.0 : (k % 2 ? 4.0 : 2.0);
            sum += w * (r * flux * Hval) * r;   // dr = r dln(r)
        }
        double r_integral = sum * hh / 3.0;
        double f_eff = (2.0 * M_PI / B_mass) * r_integral;
        // If f is 0 or NaN (z_gamma = 1), extrapolate from the two points before.
        if (std::isfinite(f_eff) && f_eff > 0.0) logf[i] = std::log(f_eff);
        else logf[i] = (i >= 2) ? 2.0 * logf[i-1] - logf[i-2] : std::log(1e-300);
        fg[i] = std::exp(logf[i]);
    }
    set_flux_spline(z_gamma_grid, fg);

    // Save the table. Write a private file and rename it, since many processes may do this at once.
    if (!getenv("FLUX_NO_CACHE")) {
        mkdir(EFF_CACHE_DIR, 0775);
        std::string tmp = cache_file + ".tmp" + std::to_string((long)getpid());
        bool ok = write_flux_table(tmp, z_gamma_grid, fg, "effective flux " + flux_source
                                   + ", made by init_effective_flux(), safe to delete", key);
        if (!ok || std::rename(tmp.c_str(), cache_file.c_str()) != 0) std::remove(tmp.c_str());
    }
}

void init_table_flux(const std::string& filename)
{
    std::vector<double> z_gamma, f;
    if (!read_flux_table(filename, z_gamma, f, nullptr))
        throw std::runtime_error("FLUX_MODEL=TABLE: could not open FLUX_FILE " + filename);
    if (z_gamma.size() < 3)
        throw std::runtime_error("FLUX_MODEL=TABLE: " + filename + " needs at least 3 lines \"z_gamma  f(z_gamma)\"");
    for (size_t i = 0; i < z_gamma.size(); i++) {
        if (!(z_gamma[i] > 0.0) || !(f[i] > 0.0) || (i > 0 && z_gamma[i] <= z_gamma[i-1]))
            throw std::runtime_error("FLUX_MODEL=TABLE: " + filename + " must have increasing z_gamma > 0 and f(z_gamma) > 0 (data line "
                                     + std::to_string(i + 1) + ")");
    }
    set_flux_spline(z_gamma, f);
    flux_zero_above = true;
    flux_source = "TABLE (" + filename + ")";
}

bool flux_is_integrated(void* p)
{
    parameters* par = (parameters*)p;
    return par->flux_model == "EFF" || par->flux_model == "TABLE";
}

double integrated_photon_flux(double z_gamma)
{
    if (!flux_spline)
        throw std::runtime_error("integrated_photon_flux() called before init_photon_flux()");
    double lz_gamma = std::log(z_gamma);
    if (lz_gamma > flux_lz_gamma_max) {
        if (flux_zero_above) return 0.0;
        lz_gamma = flux_lz_gamma_max;
    }
    if (lz_gamma < flux_lz_gamma_min) lz_gamma = flux_lz_gamma_min;
    return std::exp(gsl_spline_eval(flux_spline, lz_gamma, flux_acc));
}

void require_flux_covers(double z_gamma_lowest, void* p)
{
    parameters* par = (parameters*)p;
    if (par->flux_model != "TABLE") return;
    if (std::log(z_gamma_lowest) < flux_lz_gamma_min) {
        std::ostringstream msg;
        msg << "FLUX_MODEL=TABLE: this point needs the flux down to z_gamma = " << z_gamma_lowest
            << " but the table starts at z_gamma = " << std::exp(flux_lz_gamma_min) << ".";
        throw std::runtime_error(msg.str());
    }
}

void set_collision_parameters(void* p)
{
    parameters* par = (parameters*)p;
    // TARGET=pA: p+Pb. It changes the collision energy and the photon flux.
    par->target  = getenv("TARGET") ? getenv("TARGET") : "AA";
    par->ss      = (par->target == "pA") ? 8160.0 : 5360.0;
    par->alpha   = 1.0/137.0;
    par->Z       = 82.0;
    par->mn      = 0.938;
    par->S       = pow(17.4, 2) / pow(0.197327, 2);
    par->channel = getenv("CHANNEL") ? getenv("CHANNEL") : "An0n";
    par->gamma_aa_one = getenv("GAMMA_AA_ONE") != nullptr;
    par->bmin    = 14.2 / 0.197327;
    par->flux_model = getenv("FLUX_MODEL") ? getenv("FLUX_MODEL")
                                           : (par->target == "pA" ? "WS" : "EFF");
}

void init_photon_flux(void* p)
{
    parameters* par = (parameters*)p;
    const std::string& model = par->flux_model;
    if (par->target != "AA" && par->target != "pA")
        throw std::runtime_error("Unknown TARGET '" + par->target + "'. Expected AA or pA.");
    if (model == "TABLE") {
        if (!getenv("FLUX_FILE"))
            throw std::runtime_error("FLUX_MODEL=TABLE needs FLUX_FILE=<file with lines \"z_gamma  f(z_gamma)\">");
        init_table_flux(getenv("FLUX_FILE"));
        return;
    }
    if (model != "EFF" && model != "PL" && model != "WS")
        throw std::runtime_error("Unknown FLUX_MODEL '" + model + "'. Expected EFF, PL, WS or TABLE.");

    if (par->target == "pA") {
        // p+Pb: WS (default) or PL, with Gamma_pA(b) = exp(-sigma_NN T_A(b)).
        if (model == "EFF")
            throw std::runtime_error("FLUX_MODEL=EFF is not for TARGET=pA (the target is a proton). Use WS, PL or TABLE.");
        double sigma_NN_mb = getenv("SIGMA_NN") ? atof(getenv("SIGMA_NN")) : 99.0;
        init_pA_flux(sigma_NN_mb);
        std::ostringstream info;
        if (model == "WS") {
            init_ws_form_factor(6.49, 0.54, par->mn);
            par->bmin = 1e-3;
            info << "WS, p+Pb (Gamma_pA with sigma_NN = " << sigma_NN_mb << " mb)";
        } else {
            const double R_PL_fm = 7.1;   // b_min = 1.1 R_PL (arXiv:2606.05469 Sec. 6.1)
            par->bmin = 1.1 * R_PL_fm / 0.197327;
            info << "PL, p+Pb (b > 1.1 R_PL)";
        }
        flux_source = info.str();
        return;
    }

    load_data_and_initialize("./input/WS_photon_flux/Gamma_AA.dat");
    if (model == "WS") {
        init_ws_form_factor(6.49, 0.54, par->mn);
        flux_source = "WS (single b, Woods-Saxon form factor)";
    } else if (model == "PL") {
        flux_source = "PL (single b, point-like)";
    } else {
        init_effective_flux(par->channel, p);
    }
}

const std::string& photon_flux_info() { return flux_source; }
