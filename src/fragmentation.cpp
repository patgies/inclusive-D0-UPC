#include <cmath>
#include <memory>
#include <stdexcept>
#include "LHAPDF/LHAPDF.h"

namespace {
double g_lhapdf_scale = 1.5;
bool g_lhapdf_scale_set = false;
}

void SetLhapdfScale(double Q)
{
    g_lhapdf_scale = Q;
    g_lhapdf_scale_set = true;
}

double theta_step(double x) {
    return (x >= 0.0) ? 1.0 : 0.0;
}

// 1. Pseudoscalar state DP(z)

double BCFY_DP(double z, double r)
{ 
    if (z <= 0.0 || z >=1.0) return 0.0;
    
    double denom = 1.0 - (1.0 - r)*z;

    double factor = r * z * pow(1.-z,2) / pow(denom, 6);

   double poly = 6.0 - 18.0 * (1.0 - 2.0 * r) * z 
                  + (21.0 - 74.0 * r + 68.0 * r * r) * z * z 
                  - 2.0 * (1.0 - r) * (6.0 - 19.0 * r + 18.0 * r * r) * z * z * z 
                  + 3.0 * std::pow(1.0 - r, 2) * (1.0 - 2.0 * r + 2.0 * r * r) * z * z * z * z;
                  
    return factor * poly;
}
 

// 2. Vector state DV(z) 

double BCFY_DV(double z, double r) {
    if (z <= 0.0 || z >= 1.0) return 0.0;
    
    double denom = 1.0 - (1.0 - r) * z;

    double factor = 3.0 * (r * z * std::pow(1.0 - z, 2)) / std::pow(denom, 6); 
    
    double poly = 2.0 - 2.0 * (3.0 - 2.0 * r) * z 
                  + 3.0 * (3.0 - 2.0 * r + 4.0 * r * r) * z * z 
                  - 2.0 * (1.0 - r) * (4.0 - r + 2.0 * r * r) * z * z * z 
                  + std::pow(1.0 - r, 2) * (3.0 - 2.0 * r + 2.0 * r * r) * z * z * z * z;
                  
    return factor * poly;
}



double Dc_to_D0(double z, double r)
{
    const double mD = 1.8648;   // D0 mass in GeV
    const double mDstar = 2.0067; // D*0 mass in GeV
    const double mass_ratio = mDstar / mD;

    double DP_part = 0.168 * BCFY_DP(z, r);

    // Component 2: Vector state feed-down with kinematic scaling
    double z_scaled = mass_ratio * z;
    double DV_part = 0.0;


    double step = theta_step( (mD / mDstar) - z );
    DV_part = 0.39 * step * BCFY_DV(z_scaled, r) * mass_ratio;


    return DP_part + DV_part;
}

// Kniehl & Kramer fragmentation: D(z) = N * z * (1-z)^2 / ((1-z)^2 + eps*z)^2
double D_kniehl_kramer(double z, double N, double eps)
{
    if (z <= 0.0 || z >= 1.0) return 0.0;
    double one_minus_z = 1.0 - z;
    double denom = one_minus_z * one_minus_z + eps * z;
    return N * z * one_minus_z * one_minus_z / (denom * denom);
}

// LHAPDF fragmentation function D(z,Q) = xfxQ(pid, z, Q) / z, evaluated with the
// central member of a given LHAPDF FF set (same convention as
// partonCS/scripts/plot_fragmentation.py: D_lhapdf_cached).
// The PDF set is loaded once (lazily, on first call) and cached for reuse,
// since gsl_monte_vegas calls this O(1e5-1e7) times per point.
double D_lhapdf(double z, const std::string& setname, int pid)
{
    if (z <= 0.0 || z >= 1.0) return 0.0;

    static std::unique_ptr<LHAPDF::PDF> pdf;
    static std::string loaded_setname;

    if (!pdf || loaded_setname != setname) {
        pdf.reset(LHAPDF::mkPDF(setname, 0));
        loaded_setname = setname;
    }

    const double Q = g_lhapdf_scale_set ? g_lhapdf_scale : 1.5;
    return pdf->xfxQ(pid, z, Q) / z;
}