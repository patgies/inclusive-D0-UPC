#include "amplitudelib.hpp"
#include "def.hpp"
#include "tools.hpp"
#include "gamma_aa.hpp"
#include "interpolation.hpp"
#include "fourier.h"
#include "fragmentation.hpp"
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
    if (argc != 4) {
        cerr << "Error: expected 3 arguments, got " << argc - 1 << "." << endl;
        cerr << "Usage: " << argv[0] << " <pD0> <dipole_file> <y>" << endl;
        cerr << "  pD0         D0 meson transverse momentum [GeV], required" << endl;
        cerr << "  dipole_file path to a dipole amplitude data file, required" << endl;
        cerr << "  y           rapidity, required" << endl;
        return 1;
    }

    double pD0      = StrToReal(argv[1]);
    string datafile = argv[2];
    double y        = StrToReal(argv[3]);

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
    param.frag_type       = FragmentationType::KniehlKramer;
    param.zmin = 0.05;
    param.zmax = 1.0;

    param.alpha   = 1.0/137.0;
    param.Z       = 82.0;
    param.mn      = 0.938;
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
