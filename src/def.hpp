#ifndef def_hpp 
#define def_hpp
#include "amplitudelib.hpp"


using namespace std;

struct parameters          
{
    AmplitudeLib *dipole;   
    string datafile;  
    double xbj;
    double Nc; 
    double mc;
    double kp; // k^+
    double kD;
    int Z;
    double S;
    double ss;
    double N;
    double e;
    double mt;
    double y;
    double A;
    size_t calls; 
    double zf;
    double M;
};


double integrand(double* vec, size_t dim, void* p);

double partonicCrossSection( void* p);

double integralZf(parameters* par); 

#endif
