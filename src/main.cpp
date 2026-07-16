#include "def.hpp"
#include "tools.hpp"  
#include "fourier.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <sstream>

using namespace std;


int main(int argc, char* argv[])  
{

    string datafile = "./data/Pb/mve/glauber_mve_2";

    AmplitudeLib dipole(datafile);
    dipole.SetOutOfRangeErrors(false);
    dipole.SetInterpolationMethod(LINEAR_LINEAR);
    //dipole.SetFTMethod(ACC_SERIES);
    //set_fpu_state();
    init_workspace_fourier(1000); // max number of besselJ zeros
    set_fourier_precision(1.0e-12,1.0e-12);

    
    parameters param;
    param.dipole = &dipole;
    param.calls   = 2e5;
    param.pD0 = StrToReal(argv[1]);
    param.m   = 1.5;
    param.m2  = param.m * param.m;
    param.ss  = 5360.0;
    param.alpha   = 1.0/137.0;
    param.Z       = 82.0;
    param.mn      = 0.938;
    
    // Fragmentation function 
    param.r            = 0.1;
    param.N_kk         = 0.694;
    param.eps_kk       = 0.101;
    param.use_kniehl_kramer = false;
    param.zmin = 0.05;
    param.zmax = 1.0;
    
    // Photon flux
    param.S       = pow(17.4, 2) / pow(0.197327, 2);
    param.channel = "An0n";
    param.bmin    = 14.2 / 0.197327;
    param.bmax    = 650.0;
    param.qpmax   = 800.0;
    
    param.rmax    = 150.0;
    param.lmax    = 50.0;
    
    
   
    param.Sk_interp = inst.MakeSkInterpolator(xbj, param.lmax);
    
    double result = integralZf(&par);

    cout << par.kD << " " << result << endl;

    return 0;
}
