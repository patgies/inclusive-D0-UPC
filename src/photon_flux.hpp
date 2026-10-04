#ifndef PHOTON_FLUX_HPP
#define PHOTON_FLUX_HPP
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include <string>
#include <vector>
#include <stdexcept>
#include <fstream>

// Woods-Saxon charge form factor F_A(q), formula of Maximon & Schrack,
// J. Res. Natl. Bur. Stand. B 70, 85 (1966). F(0) = 1.
class WSFormFactor {
private:
    double c_, a_, q_max_, rho0_;
    static const int N_TERMS = 10;

    static double series_cube(double c, double a)
    {
        double sum = 0.0, sign = 1.0;
        for (int n = 1; n <= N_TERMS; n++) {
            sum += sign * std::exp(-n * c / a) / (n * n * n);
            sign = -sign;
        }
        return sum;
    }

    static double series_F(double q, double c, double a)
    {
        double qa2 = (q * a) * (q * a);
        double sum = 0.0, sign = 1.0;
        for (int n = 1; n <= N_TERMS; n++) {
            double denom = (n * n + qa2);
            sum += sign * n * std::exp(-n * c / a) / (denom * denom);
            sign = -sign;
        }
        return sum;
    }

public:
    // RA and a in GeV^-1, q_max in GeV.
    WSFormFactor(double RA, double a, double q_max) : c_(RA), a_(a), q_max_(q_max)
    {
        double bracket = (4.0 * M_PI * c_ / 3.0) * (M_PI * a_ * M_PI * a_ + c_ * c_)
                        + 8.0 * M_PI * a_ * a_ * a_ * series_cube(c_, a_);
        rho0_ = 1.0 / bracket;
    }

    WSFormFactor(const WSFormFactor&) = delete;
    WSFormFactor& operator=(const WSFormFactor&) = delete;

    double operator()(double q) const
    {
        if (q >= q_max_) return 0.0;      // F is already zero before q_max
        if (q < 1e-6) return 1.0;         // avoids 0/0 in the formula below
        double qa  = q * a_;
        double pqa = M_PI * qa;
        double sh = std::sinh(pqa), ch = std::cosh(pqa);
        double elementary = (4.0 * M_PI * M_PI * rho0_ * a_ * a_ * a_) / (qa * qa * sh * sh)
                           * (pqa * ch * std::sin(q * c_) - q * c_ * std::cos(q * c_) * sh);
        double series = 8.0 * M_PI * rho0_ * a_ * a_ * a_ * series_F(q, c_, a_);
        return elementary + series;
    }
};

#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_sf_result.h>

// Table of the k_perp integral I(z_gamma,r) of the Woods-Saxon flux, on a grid in (log z_gamma, log r).
// f_WS(z_gamma,r) = alpha Z^2/(pi^2 z_gamma) * I(z_gamma,r)^2 (arXiv:2606.05469 Eq. 16).
class WSFluxTable {
private:
    gsl_spline2d* spline;
    gsl_interp_accel* z_gamma_acc;
    gsl_interp_accel* r_acc;
    double lz_gamma_min, lz_gamma_max, lr_min, lr_max;

    struct IntegrandParams { double z_gamma, r, mn2; const WSFormFactor* FA; };

    static double integrand(double kt, void* p)
    {
        IntegrandParams* par = (IntegrandParams*)p;
        double t = kt * kt + par->z_gamma * par->z_gamma * par->mn2;
        gsl_sf_result res;
        int status = gsl_sf_bessel_Jn_e(1, kt * par->r, &res);
        double J1 = (status == GSL_SUCCESS) ? res.val : 0.0;
        return kt * kt * (*par->FA)(std::sqrt(t)) / t * J1;
    }

    static double raw_I(double z_gamma, double r, double mn, const WSFormFactor& FA, double q_max)
    {
        IntegrandParams par{z_gamma, r, mn * mn, &FA};
        gsl_function F; F.function = &integrand; F.params = &par;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(2000);
        double result, err;
        gsl_integration_qag(&F, 0.0, q_max, 0, 1e-5, 2000, GSL_INTEG_GAUSS61, w, &result, &err);
        gsl_integration_workspace_free(w);
        return result;
    }

public:
    WSFluxTable(const WSFormFactor& FA, double mn, double q_max,
                double z_gamma_min, double z_gamma_max, int nz_gamma,
                double r_min, double r_max, int nr)
    {
        lz_gamma_min = std::log(z_gamma_min); lz_gamma_max = std::log(z_gamma_max);
        lr_min = std::log(r_min); lr_max = std::log(r_max);

        std::vector<double> lz_gamma(nz_gamma), lr(nr);
        for (int i = 0; i < nz_gamma; i++) lz_gamma[i] = lz_gamma_min + (lz_gamma_max - lz_gamma_min) * i / (nz_gamma - 1.0);
        for (int j = 0; j < nr; j++) lr[j] = lr_min + (lr_max - lr_min) * j / (nr - 1.0);

        spline = gsl_spline2d_alloc(gsl_interp2d_bicubic, nz_gamma, nr);
        std::vector<double> grid_z(nz_gamma * nr);
        for (int i = 0; i < nz_gamma; i++) {
            double z_gamma = std::exp(lz_gamma[i]);
            for (int j = 0; j < nr; j++) {
                double r = std::exp(lr[j]);
                double I = raw_I(z_gamma, r, mn, FA, q_max);
                double logI2 = std::log(std::max(I * I, 1e-300));
                gsl_spline2d_set(spline, grid_z.data(), i, j, logI2);
            }
        }
        gsl_spline2d_init(spline, lz_gamma.data(), lr.data(), grid_z.data(), nz_gamma, nr);
        z_gamma_acc = gsl_interp_accel_alloc();
        r_acc = gsl_interp_accel_alloc();
    }

    WSFluxTable(const WSFluxTable&) = delete;
    WSFluxTable& operator=(const WSFluxTable&) = delete;

    // I(z_gamma,r)^2, without the prefactor.
    double I2(double z_gamma, double r) const
    {
        double lz_gamma_ = std::log(z_gamma), lr_ = std::log(r);
        if (lz_gamma_ < lz_gamma_min) lz_gamma_ = lz_gamma_min; else if (lz_gamma_ > lz_gamma_max) lz_gamma_ = lz_gamma_max;
        if (lr_ < lr_min) lr_ = lr_min; else if (lr_ > lr_max) lr_ = lr_max;
        return std::exp(gsl_spline2d_eval(spline, lz_gamma_, lr_, z_gamma_acc, r_acc));
    }

    ~WSFluxTable()
    {
        gsl_spline2d_free(spline);
        gsl_interp_accel_free(z_gamma_acc);
        gsl_interp_accel_free(r_acc);
    }
};


class TASpline {
private:
    gsl_spline* spline;
    gsl_interp_accel* acc;
    double s_min, s_max;

public:
    TASpline(const std::vector<double>& s,
             const std::vector<double>& T)
    {
        spline = gsl_spline_alloc(gsl_interp_cspline, s.size());
        acc = gsl_interp_accel_alloc();
        gsl_spline_init(spline, s.data(), T.data(), s.size());

        s_min = s.front();
        s_max = s.back();
    }


    TASpline(const TASpline&) = delete;
    TASpline& operator=(const TASpline&) = delete;

    double operator()(double s)
    {
        if (s < s_min || s > s_max)
            return 1.0;

        return gsl_spline_eval(spline, s, acc);
    }

    ~TASpline()
    {
        gsl_spline_free(spline);
        gsl_interp_accel_free(acc);
    }
};


class GammaAA {
private:
    TASpline& TA;
    double b_min;
    double b_max;

public:
    GammaAA(TASpline& ta)
        : TA(ta), b_min(0.0), b_max(1e6)
    {}

    double operator()(double b)
    {
        if (b < b_min)
            throw std::out_of_range("b too small");

        if (b > b_max)
            return 1.0;

        return TA(b);
    }
};

GammaAA& gamma_aa();

// Thickness function T_B(s) of the Woods-Saxon density, with int d2s T_B(s) = B_mass.
class WSThickness {
private:
    gsl_spline* spline;
    gsl_interp_accel* acc;
    double s_max_;

public:
    WSThickness(double RA, double a, double B_mass, double s_max, int n = 150)
        : s_max_(s_max)
    {
        double zmax = RA + 12.0 * a;
        std::vector<double> sg(n), raw(n);
        for (int i = 0; i < n; i++) {
            double s = s_max * i / (n - 1.0);
            sg[i] = s;
            int M = 1000; double h = zmax / M, sum = 0.0;
            for (int k = 0; k <= M; k++) {
                double zp = k * h;
                double r  = std::sqrt(zp * zp + s * s);
                double rho = 1.0 / (1.0 + std::exp((r - RA) / a));
                double w = (k == 0 || k == M) ? 1.0 : (k % 2 ? 4.0 : 2.0);
                sum += w * rho;
            }
            raw[i] = 2.0 * sum * h / 3.0;   // x2: the integral is symmetric in z'
        }
        double h2 = sg[1] - sg[0], sum2 = 0.0;
        for (int i = 0; i < n; i++) {
            double w = (i == 0 || i == n - 1) ? 1.0 : (i % 2 ? 4.0 : 2.0);
            sum2 += w * raw[i] * sg[i];
        }
        double integral = 2.0 * M_PI * sum2 * h2 / 3.0;
        std::vector<double> Tg(n);
        for (int i = 0; i < n; i++) Tg[i] = raw[i] * (B_mass / integral);

        spline = gsl_spline_alloc(gsl_interp_cspline, n);
        acc = gsl_interp_accel_alloc();
        gsl_spline_init(spline, sg.data(), Tg.data(), n);
    }
    WSThickness(const WSThickness&) = delete;
    WSThickness& operator=(const WSThickness&) = delete;
    double s_max() const { return s_max_; }
    double operator()(double s) const { return (s >= s_max_) ? 0.0 : gsl_spline_eval(spline, s, acc); }
    ~WSThickness() { gsl_spline_free(spline); gsl_interp_accel_free(acc); }
};

// H(r) = int d2s T_B(s) Gamma_AA(|r-s|) P_EMD(|r-s|), the geometric factor of the
// effective flux (arXiv:2404.09731 Eq. 4). Spline in r up to r_switch; above it
// H(r) = B_mass Gamma_AA(r) P_EMD(r).
class EffFluxRadial {
private:
    gsl_spline* spline;
    gsl_interp_accel* acc;
    double r_min_, r_switch_;
    double B_mass, S;
    std::string channel;

    static double channel_factor(const std::string& channel, double b, double S)
    {
        double P_no_emd = std::exp(-S / (b * b));
        if (channel == "An0n") return P_no_emd;
        if (channel == "Xn0n") return P_no_emd * (1.0 - P_no_emd);
        return 1.0;   // AnAn
    }

    static double H_raw(double r, const WSThickness& TB, const std::string& channel, double S)
    {
        int Ns = 150, Nphi = 100;
        double sh = TB.s_max() / Ns, sum_s = 0.0;
        for (int i = 0; i <= Ns; i++) {
            double s = i * sh;
            double Tval = TB(s);
            double ph = M_PI / Nphi, sum_phi = 0.0;
            for (int j = 0; j <= Nphi; j++) {
                double phi = j * ph;
                double b = std::sqrt(r * r + s * s - 2.0 * r * s * std::cos(phi));
                double val = gamma_aa()(b) * channel_factor(channel, b, S);
                double w = (j == 0 || j == Nphi) ? 1.0 : (j % 2 ? 4.0 : 2.0);
                sum_phi += w * val;
            }
            double phi_integral = 2.0 * (sum_phi * ph / 3.0);   // x2: symmetric about phi = pi
            double w = (i == 0 || i == Ns) ? 1.0 : (i % 2 ? 4.0 : 2.0);
            sum_s += w * s * Tval * phi_integral;
        }
        return sum_s * sh / 3.0;
    }

public:
    // Log grid in r from r0 to r_switch.
    EffFluxRadial(const WSThickness& TB, const std::string& channel_, double S_,
                  double B_mass_, double r_switch, int nr = 400, double r0 = 0.5)
        : r_min_(r0), r_switch_(r_switch), B_mass(B_mass_), S(S_), channel(channel_)
    {
        std::vector<double> lrg(nr), Hg(nr);
        double lr0 = std::log(r0), lr1 = std::log(r_switch);
        for (int i = 0; i < nr; i++) {
            double lr = lr0 + (lr1 - lr0) * i / (nr - 1.0);
            lrg[i] = lr;
            Hg[i] = H_raw(std::exp(lr), TB, channel, S);
        }
        spline = gsl_spline_alloc(gsl_interp_cspline, nr);
        acc = gsl_interp_accel_alloc();
        gsl_spline_init(spline, lrg.data(), Hg.data(), nr);
    }
    EffFluxRadial(const EffFluxRadial&) = delete;
    EffFluxRadial& operator=(const EffFluxRadial&) = delete;

    double operator()(double r) const
    {
        if (r >= r_switch_) return B_mass * gamma_aa()(r) * channel_factor(channel, r, S);
        double rc = (r < r_min_) ? r_min_ : r;
        return gsl_spline_eval(spline, std::log(rc), acc);
    }

    ~EffFluxRadial() { gsl_spline_free(spline); gsl_interp_accel_free(acc); }
};


// All the fluxes are f = dN/dz_gamma, with z_gamma = omega/E_beam.
// The factor dz_gamma/dq+ is applied in integrand_inclusive().
//
// Transverse distances:
//   r   : from the centre of the emitting nucleus to the point where the photon is absorbed
//   s   : from the centre of the target nucleus to that point
//   b   : between the centres of the two nuclei, b = |r - s|
//   b_d : dipole impact parameter in the target (the glauber_mve_<b_d> files), not used here

// Sets sqrt(s), the nucleus and the flux settings (TARGET, CHANNEL, FLUX_MODEL, GAMMA_AA_ONE).
void set_collision_parameters(void* par);

// Prepares the flux of par->flux_model:
//   EFF   : effective flux, arXiv:2404.09731 Eq. 4 (default)
//   PL    : point-like flux, integrated over b in the cross section
//   WS    : Woods-Saxon flux, integrated over b in the cross section
//   TABLE : f(z_gamma) from the file FLUX_FILE, lines "z_gamma  f(z_gamma)"
// TARGET=pA (p+Pb): WS (default) or PL with Gamma_pA(b) and no EMD factor, or TABLE.
void init_photon_flux(void* par);

// Name of the flux in use, for the output header.
const std::string& photon_flux_info();

// true if the flux is already integrated over b (EFF, TABLE).
bool flux_is_integrated(void* par);

// f(z_gamma) for EFF and TABLE.
double integrated_photon_flux(double z_gamma);

// TABLE: error if the table starts above z_gamma_lowest.
void require_flux_covers(double z_gamma_lowest, void* par);

// PL and WS: integrand of the b integral, 2 pi b f(z_gamma,r=b) Gamma_AA(b) P_EMD(b).
double photon_flux(double b, double z_gamma, void* par);

// Point-like and Woods-Saxon f(z_gamma,r).
double flux_density(double z_gamma, double r, void* par);
double flux_density_WS(double z_gamma, double r, void* par);

// p+Pb: Gamma_pA(b) = exp(-sigma_NN T_A(b)) (arXiv:2606.05469 Eq. 26).
void init_pA_flux(double sigma_NN_mb, double RA_fm = 6.49, double a_fm = 0.54, double B_mass = 208.0);
double GammaPA(double b);

// Reads the table "b  Gamma_AA(b)".
void load_data_and_initialize(const std::string& filename);

// Builds the Woods-Saxon form factor and flux table.
void init_ws_form_factor(double RA_fm = 6.49, double a_fm = 0.54, double mn = 0.938);

// Builds f_eff(z_gamma), or reads it from input/flux_cache/.
void init_effective_flux(const std::string& channel, void* par,
                          double RA_fm = 6.49, double a_fm = 0.54, double B_mass = 208.0);

// Reads f(z_gamma) from a file with lines "z_gamma  f(z_gamma)".
void init_table_flux(const std::string& filename);

// The Gamma_AA(b) table read by load_data_and_initialize().
GammaAA& gamma_aa();

#endif
