#include "amplitudelib.hpp"
#include "def.hpp"
#include "tools.hpp"
#include "photon_flux.hpp"
#include "interpolation.hpp"
#include "fourier.h"
#include "hymnd_grid.hpp"
#include "kk_grid.hpp"
#include "bcfy_grid.hpp"
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <gsl/gsl_errno.h>

using namespace std;

int main(int argc, char* argv[])
{
    if (argc != 3 && argc != 4) {
        cerr << "Error: expected 2 or 3 arguments, got " << argc - 1 << "." << endl;
        cerr << "Usage: " << argv[0] << " <pD0> [<dipole_file>] <y>" << endl;
        cerr << "  pD0  D0 meson transverse momentum [GeV], required" << endl;
        cerr << "  dipole_file path to a dipole amplitude data file;" << endl;
        cerr << "  if omitted, read from the DIPOLE_FILE environment variable" << endl;
        cerr << "  y  rapidity, required" << endl;
        return 1;
    }

    double pD0 = StrToReal(argv[1]);

    string datafile;
    if (argc == 4) {
        datafile = argv[2];
    } else if (getenv("DIPOLE_FILE")) {
        datafile = getenv("DIPOLE_FILE");
    }

    double y = StrToReal(argv[argc - 1]);

    if (datafile.empty()) {
        cerr << "Error: no dipole_file given and DIPOLE_FILE environment variable is not set." << endl;
        return 1;
    }

    AmplitudeLib inst(datafile);
    inst.SetOutOfRangeErrors(false);
    inst.SetInterpolationMethod(LINEAR_LINEAR);

    // Some dipole files have no x0 in the header.
    if (getenv("DIPOLE_X0")) {
        inst.SetX0(StrToReal(getenv("DIPOLE_X0")));
    }

    gsl_set_error_handler_off();

    parameters param;
    param.dipole = &inst;

    param.pD0 = pD0;
    param.m   = 1.5;
    param.m2  = param.m * param.m;
    set_collision_parameters(&param);

    param.r            = 0.1;
    param.N_kk         = 0.694;
    param.eps_kk       = 0.101;
    param.frag_type = FragmentationType::KniehlKramer;
    if (getenv("FRAG_TYPE")) {
        string frag_type_env = getenv("FRAG_TYPE");
        if (frag_type_env == "BCFY") {
            param.frag_type = FragmentationType::BCFY;
        } else if (frag_type_env == "KniehlKramer") {
            param.frag_type = FragmentationType::KniehlKramer;
        } else if (frag_type_env == "HymnD") {
            param.frag_type = FragmentationType::HymnD;
        } else {
            cerr << "Error: unknown FRAG_TYPE '" << frag_type_env << "'. Expected BCFY, KniehlKramer or HymnD." << endl;
            return 1;
        }
    }

    double mt0 = sqrt(param.pD0*param.pD0 + param.m2);

    double scale_factor = getenv("SCALE_FACTOR") ? StrToReal(getenv("SCALE_FACTOR")) : 1.0;
    double frag_scale = scale_factor * mt0;
    if (frag_scale < param.m) frag_scale = param.m;

    string hymnD_file = getenv("HYMND_FILE")
        ? getenv("HYMND_FILE")
        : "input/HymnD/prompt-D0-1-109_0000.dat";
    const int hymnD_charm_flavor = 4;
    if (param.frag_type == FragmentationType::HymnD) {
        param.D_frag_interp = MakeHymnDZInterpolator(hymnD_file, hymnD_charm_flavor, frag_scale);
    }
    if (param.frag_type == FragmentationType::KniehlKramer) {
        param.D_frag_interp = MakeKniehlKramerInterpolator(frag_scale);
    }
    if (param.frag_type == FragmentationType::BCFY) {
        param.D_frag_interp = MakeBCFYInterpolator(frag_scale);
    }
    param.z_h_min = 0.05;
    param.z_h_max = 1.0;

    try {
        init_photon_flux(&param);
    } catch (const std::exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    param.qpmax   = 800.0;
    param.lmax    = 50.0;
    param.calls   = getenv("CALLS") ? (size_t)StrToReal(getenv("CALLS")) : (size_t)2e5;

    // Upper limit of the b integral: the flux is zero beyond 60/(z_gamma*mn).
    // z_gamma_min is the smallest photon energy fraction possible (z_h = 1).
    double z_gamma_min = mt0 * exp(y) / param.ss;
    param.bmax   = 60.0 / (z_gamma_min * param.mn);
    try {
        require_flux_covers(z_gamma_min, &param);
    } catch (const std::exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    init_workspace_fourier(1000);
    set_fourier_precision(1.0e-6, 1.0e-6);

    cout << "# fragmentation : ";
    switch (param.frag_type) {
        case FragmentationType::KniehlKramer:
            cout << "Kniehl & Kramer, evolved (scale=" << scale_factor << "*mt0=" << frag_scale << ")";
            break;
        case FragmentationType::HymnD:
            cout << "HymnD (" << hymnD_file << ", member 0, Q=" << scale_factor << "*mt0=" << frag_scale << ")";
            break;
        case FragmentationType::BCFY:
        default:
            cout << "BCFY (r=" << param.r << "), evolved (scale=" << scale_factor << "*mt0=" << frag_scale << ")";
            break;
    }
    cout << endl;
    cout << "# target        : " << param.target << " (sqrt(s_NN) = " << param.ss << " GeV)" << endl;
    cout << "# channel       : " << param.channel << endl;
    cout << "# flux          : " << photon_flux_info() << endl;
    cout << "# y  dsigma_dyd^2pD0" << endl;

    param.y = y;

    double xbj = (mt0 / param.ss) * exp(-y);
    if (xbj > 0.01) xbj = 0.01;
    param.xbj  = xbj;

    param.Sk_interp = inst.MakeSkInterpolator(xbj, param.lmax);

    double res = D0CrossSection_inclusive(static_cast<void*>(&param));
    cout << y << "  " << res << endl;

    return 0;
}
