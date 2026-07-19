#include "amplitudelib.hpp"
#include "def.hpp"
#include "tools.hpp"
#include "gamma_aa.hpp"
#include "interpolation.hpp"
#include "fourier.h"
#include "fragmentation.hpp"
#include "lhapdf_grid.hpp"
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

    gsl_set_error_handler_off();

    load_data_and_initialize("./data/Gamma_AA.dat");

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
        } else if (frag_type_env == "LHAPDF") {
            param.frag_type = FragmentationType::LHAPDF;
        } else {
            cerr << "Error: unknown FRAG_TYPE '" << frag_type_env << "'. Expected BCFY, KniehlKramer or LHAPDF." << endl;
            return 1;
        }
    }

    // Grid file for FragmentationType::LHAPDF: member 0 (central value) of an
    // LHAPDF lhagrid1-format fragmentation function set. Override with LHAPDF_FILE.
    string lhapdf_file = getenv("LHAPDF_FILE")
        ? getenv("LHAPDF_FILE")
        : "data/prompt-D0-1-109/prompt-D0-1-109_0000.dat";
    const int lhapdf_charm_flavor = 4;  // PDG id for the charm quark
    if (param.frag_type == FragmentationType::LHAPDF) {
        param.D_frag_interp = MakeLHAPDFZInterpolator(lhapdf_file, lhapdf_charm_flavor, param.m);
    }
    param.zmin = 0.05;
    param.zmax = 1.0;

    param.alpha   = 1.0/137.0;
    param.Z       = 82.0;
    param.mn      = (208 * 0.931) / 208;  // M/A, matching old-inclusiveD0's gamma_L = (ss/2)/(M/A)
    param.S       = pow(17.4, 2) / pow(0.197327, 2);
    param.channel = getenv("CHANNEL") ? getenv("CHANNEL") : "An0n";

    param.bmin    = 14.2 / 0.197327;
    param.bmax    = 650.0;
    param.qpmax   = 800.0;
    param.rmax    = 150.0;
    param.lmax    = 50.0;
    param.calls   = 2e5;

    double mt0 = sqrt(param.pD0*param.pD0 + param.m2);

    init_workspace_fourier(1000);
    set_fourier_precision(1.0e-6, 1.0e-6);

    cout << "# fragmentation : ";
    switch (param.frag_type) {
        case FragmentationType::KniehlKramer:
            cout << "Kniehl & Kramer (N=" << param.N_kk << ", eps=" << param.eps_kk << ")";
            break;
        case FragmentationType::LHAPDF:
            cout << "LHAPDF (" << lhapdf_file << ", member 0, Q=" << param.m << ")";
            break;
        case FragmentationType::BCFY:
        default:
            cout << "BCFY (r=" << param.r << ")";
            break;
    }
    cout << endl;
    cout << "# channel       : " << param.channel << endl;
    cout << "# y  dsigma_dy" << endl;

    param.y = y;

    double xbj = (mt0 / param.ss) * exp(-y);
    if (xbj > 0.01) xbj = 0.01;
    param.xbj  = xbj;
    param.qpmin = (mt0 / sqrt(2.0)) * exp(y);

    param.Sk_interp = inst.MakeSkInterpolator(xbj, param.lmax);

    double res = D0CrossSection_inclusive(static_cast<void*>(&param));
    cout << y << "  " << res << endl;

    return 0;
}
