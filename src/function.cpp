#include "def.hpp"
#include "tools.hpp"
#include "amplitudelib.hpp"
#include <cmath>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_bessel.h>
#include <fstream>
#include <iostream>

using namespace std;



// Hard amplitude in momentum space 
// pc  = charm quark transverse momentum (= pD0/zh)
// kp  = charm quark plus-momentum = (mT/sqrt(2))*exp(y)
// l   = gluon transverse momentum magnitude
// phi = angle between pc and l vectors
// qp  = photon plus-momentum
static double F_hard(double m, double kp, double pc, double l, double phi, double qp)
{
    double z   = kp / qp;
    double abr = pc*pc + l*l - 2.0*pc*l*cos(phi);  // |pc - l|^2

    double den1 = (m*m + abr) * (m*m + abr);
    double den2 = (m*m + pc*pc) * (m*m + pc*pc);
    double den3 = (m*m + abr) * (m*m + pc*pc);

    double par1 = 1.0/den1 + 1.0/den2 - 2.0/den3;
    double par2 = abr/den1 + pc*pc/den2 - 2.0*(pc*pc - l*pc*cos(phi))/den3;

    return l * (2.0*z*m*m*par1 + 2.0*z*(z*z + (1.0-z)*(1.0-z))*par2);
}


//  Collinear approximation expression

double col_ap(double m, double kp, double k, double l, double qp)  
{
   double z = kp / qp;
   
   double num = 4*SQR(m)*SQR(k) + 2*(SQR(z) + SQR(1-z))*(SQR(m)*SQR(m)+SQR(k)*SQR(k));
   double den = SQR(k*k + m*m)*SQR(k*k + m*m);

   return l*l*l*(z*num/den);

}

// Photon flux including 0n emission probability

double modflux(double qp, double b, int Z, double S, double ss, double M, double A)
{
     double alpha =1./137.;
     double w = qp /sqrt(2.);
     
     //0n neutron emission
     double p = exp(-S/(b*b));
     double gamma= (ss/2.)/(M/A);
     //flux
     double flux = (1./qp)*(SQR(Z)*alpha * SQR(w) / (SQR(M_PI) * SQR(gamma)) ) * SQR(besselFunction(1,w*b*5.068/gamma)); 
     
     return 2.*M_PI*b*p*flux;    //2pi and b from the 2D integral 

}
     

double integrand(double* vec, size_t dim, void* p)
{
    parameters *par = ((parameters*)p);   

    double x =  vec[0];
    double l =  vec[1];
    double qp =  vec[2];
    double b =  vec[3];    
    
    //Dipole amplitudes
    double D = par -> dipole -> S_k( l, par -> xbj, FUNDAMENTAL, 1.0);
    double photflux = modflux(qp, b, par -> Z, par -> S, par -> ss, par ->M, par->A);
    double function = F(par -> mc, par -> kp, par -> kD / par -> zf, l, x, qp);
    
    double fun2 = col_ap(par -> mc, par -> kp, par -> kD/ par -> zf, l, qp);

    //cout << "D: " << D << ", photflux: " << photflux << ", function: " << function << endl;


     return D*photflux*function;
}



