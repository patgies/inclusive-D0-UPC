#include "amplitudelib.hpp"
#include "def.hpp"
#include "tools.hpp"
#include "gamma_aa.hpp"
#include "interpolation.hpp"
#include "fourier.h"
#include "fragmentation.hpp"
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <gsl/gsl_errno.h>

using namespace std;

int main(int argc, char* argv[])
{
    string datafile = (argc > 2) ? argv[2] : "./data/proton/mve.dat";

    AmplitudeLib inst(datafile);
    inst.SetOutOfRangeErrors(false);
    inst.SetInterpolationMethod(LINEAR_LINEAR);

    gsl_set_error_handler_off();

    load_data_and_initialize("./data/Gamma_AA.dat");

    parameters param;
    param.dipole = &inst;

    param.pD0 = StrToReal(argv[1]);
    param.m   = 1.5;
    param.m2  = param.m * param.m;
    param.ss  = 5360.0;

    param.r            = 0.1;
    param.N_kk         = 0.694;
    param.eps_kk       = 0.101;
    param.frag_type       = FragmentationType::LHAPDF;
    param.lhapdf_setname  = "prompt-D0-1-109";
    param.lhapdf_pid      = 4;      // charm quark
    param.lhapdf_Q        = 1.5;    
    SetLhapdfScale(param.lhapdf_Q);
    param.zmin = 0.05;
    param.zmax = 1.0;

    param.alpha   = 1.0/137.0;
    param.Z       = 82.0;
    param.mn      = 0.938;
    param.S       = pow(17.4, 2) / pow(0.197327, 2);
    param.channel = "An0n";

    param.bmin    = 14.2 / 0.197327;
    param.bmax    = 650.0;
    param.qpmax   = 800.0;
    param.rmax    = 150.0;
    param.lmax    = 50.0;
    param.calls   = 2e5;

    double mt0 = sqrt(param.pD0*param.pD0 + param.m2);

    vector<double> y_vals = (argc > 3)
        ? vector<double>{ StrToReal(argv[3]) }
        : vector<double>{0.0, 0.5, 1.0, 1.5,2.0, 2.5, 3.0, 3.5,4.0};

    init_workspace_fourier(1000);
    set_fourier_precision(1.0e-6, 1.0e-6);

    cout << "# ============================================================" << endl;
    cout << "# Inclusive D0 photoproduction cross section: gamma + A -> D0 + X" << endl;
    cout << "# ------------------------------------------------------------" << endl;
    cout << "# Parameters used:" << endl;
    cout << "#   dipole amplitude file : " << datafile << endl;
    cout << "#   charm mass m          : " << param.m << " GeV" << endl;
    cout << "#   sqrt(s)               : " << param.ss << " GeV" << endl;
    cout << "#   fragmentation         : ";
    switch (param.frag_type) {
        case FragmentationType::KniehlKramer:
            cout << "Kniehl & Kramer (N=" << param.N_kk << ", eps=" << param.eps_kk << ")";
            break;
        case FragmentationType::LHAPDF:
            cout << "LHAPDF (" << param.lhapdf_setname << ", pid=" << param.lhapdf_pid
                 << ", Q=" << param.lhapdf_Q << " GeV)";
            break;
        case FragmentationType::BCFY:
        default:
            cout << "BCFY (r=" << param.r << ")";
            break;
    }
    cout << endl;
    cout << "#   z range (frag.)       : [" << param.zmin << ", " << param.zmax << "]" << endl;
    cout << "#   alpha_EM              : " << param.alpha << endl;
    cout << "#   Z (nucleus charge)    : " << param.Z << endl;
    cout << "#   nucleon mass mn       : " << param.mn << " GeV" << endl;
    cout << "#   EMD channel           : " << param.channel << endl;
    cout << "#   impact param b range  : [" << param.bmin << ", " << param.bmax << "] GeV^-1" << endl;
    cout << "#   qp max                : " << param.qpmax << " GeV" << endl;
    cout << "#   gluon l max           : " << param.lmax << " GeV" << endl;
    cout << "#   VEGAS calls           : " << param.calls << endl;
    cout << "#   pD0                   : " << param.pD0 << " GeV" << endl;
    cout << "# ============================================================" << endl;
    cout << "# y  dsigma_dy" << endl;

    for (double y : y_vals) {
        param.y = y;

        double xbj = (mt0 / param.ss) * exp(-y);
        if (xbj > 0.01) xbj = 0.01;
        param.xbj  = xbj;
        param.qpmin = (mt0 / sqrt(2.0)) * exp(y);

        // Precompute S_k at this y's xbj (sequential: Fourier workspace not thread-safe)
        param.Sk_interp = inst.MakeSkInterpolator(xbj, param.lmax);

        double res = D0CrossSection_inclusive(static_cast<void*>(&param));
        cout << y << "  " << res << endl;
    }

    return 0;
}
