#include "amplitudelib.hpp"
#include "def.hpp"
#include "tools.hpp"
#include "gamma_aa.hpp"
#include "interpolation.hpp"
#include "fourier.h"
#include "fragmentation.hpp"
#include "hymnd_grid.hpp"
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
        cerr << "  pD0         D0 meson transverse momentum [GeV], required" << endl;
        cerr << "  dipole_file path to a dipole amplitude data file;" << endl;
        cerr << "              if omitted, read from the DIPOLE_FILE environment variable" << endl;
        cerr << "  y           rapidity, required" << endl;
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

    // Some dipole files miss x0 in the header, so fix it here if needed.
    if (getenv("DIPOLE_X0")) {
        inst.SetX0(StrToReal(getenv("DIPOLE_X0")));
    }

    gsl_set_error_handler_off();

    load_data_and_initialize("./inputs/Gamma_AA.dat");

    parameters param;
    param.dipole = &inst;

    param.pD0 = pD0;
    param.m   = 1.5;
    param.m2  = param.m * param.m;
    param.ss  = 5360.0;

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
        : "inputs/prompt-D0-1-109/prompt-D0-1-109_0000.dat";
    const int hymnD_charm_flavor = 4;
    if (param.frag_type == FragmentationType::HymnD) {
        param.D_frag_interp = MakeHymnDZInterpolator(hymnD_file, hymnD_charm_flavor, frag_scale);
    }
    param.zmin = 0.05;
    param.zmax = 1.0;

    param.alpha   = 1.0/137.0;
    param.Z       = 82.0;
    param.mn      = (208 * 0.931) / 208;
    param.S       = pow(17.4, 2) / pow(0.197327, 2);
    param.channel = getenv("CHANNEL") ? getenv("CHANNEL") : "An0n";

    param.bmin    = 14.2 / 0.197327;
    param.bmax    = 650.0;
    param.qpmax   = 800.0;
    param.lmax    = 50.0;
    param.calls   = 2e5;

    init_workspace_fourier(1000);
    set_fourier_precision(1.0e-6, 1.0e-6);

    cout << "# fragmentation : ";
    switch (param.frag_type) {
        case FragmentationType::KniehlKramer:
            cout << "Kniehl & Kramer (N=" << param.N_kk << ", eps=" << param.eps_kk << ")";
            break;
        case FragmentationType::HymnD:
            cout << "HymnD (" << hymnD_file << ", member 0, Q=" << scale_factor << "*mt0=" << frag_scale << ")";
            break;
        case FragmentationType::BCFY:
        default:
            cout << "BCFY (r=" << param.r << ")";
            break;
    }
    cout << endl;
    cout << "# channel       : " << param.channel << endl;
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
