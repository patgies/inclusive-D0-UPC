#include "def.hpp"
#include "amplitudelib.hpp"
#include "fragmentation.hpp"
#include "interpolation.hpp"
#include <cmath>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>

using namespace std;


// Hard amplitude in momentum space — same formula as inclusive-D0-UPC/src/function.cpp.
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



// par->Sk_interp must be set by the caller before launching VEGAS.

double integrand_inclusive(double* vec, size_t /*dim*/, void* p)
{
    parameters* par = (parameters*)p;

    double zh   = vec[0];
    double u_qp = vec[1];   // rescaled to [pp_zh, qpmax]
    double b    = vec[2];
    double phi  = vec[3];
    double l    = vec[4];

    //Kinematics from zh
    double pc    = par->pD0 / zh;
    double mt    = sqrt(pc*pc + par->m2);
    double pp_zh = (mt / sqrt(2.0)) * exp(par->y);

    //Rescale qp
    double jac_qp = par->qpmax - pp_zh;
    if (jac_qp <= 0.0) return 0.0;
    double qp = pp_zh + u_qp * jac_qp;

    //Photon flux 
    double flux = photon_flux(b, qp, p);

    //Gluon TMD: S_k from precomputed grid (set by caller)
    if (!par->Sk_interp) return 0.0;
    double Sk = par->Sk_interp->Evaluate(l);
    if (Sk < 0.0) Sk = 0.0;

    // Inclusive hard factor 
    double fhard = F_hard(par->m, pp_zh, pc, l, phi, qp);

    // Fragmentation function ( KK or BCFY)
    double D_frag = par->use_kniehl_kramer
                    ? D_kniehl_kramer(zh, par->N_kk, par->eps_kk)
                    : Dc_to_D0(zh, par->r);
    double frag_weight = D_frag / (zh * zh);

    return jac_qp * frag_weight * flux * Sk * fhard;
}


// Assumes par->Sk_interp is already set (precomputed in main before parallel launch).
double D0CrossSection_inclusive(void* p)
{
    parameters* par = (parameters*)p;

    const gsl_rng_type* T;
    gsl_rng* rng;
    gsl_rng_env_setup();
    T   = gsl_rng_default;
    rng = gsl_rng_alloc(T);

    gsl_monte_function F;
    F.f      = &integrand_inclusive;
    F.dim    = 5;
    F.params = par;

    // {zh, u_qp, b, phi, l}
    double low[] = {par->zmin, 0.0, par->bmin, 0.0,  0.0 };
    double up[]  = {par->zmax, 1.0, par->bmax, 2.0*M_PI, par->lmax };

    double res, err;
    gsl_monte_vegas_state* s = gsl_monte_vegas_alloc(F.dim);

    gsl_monte_vegas_integrate(&F, low, up, F.dim, par->calls/10, rng, s, &res, &err);

    int iter = 0;
    do {
        gsl_monte_vegas_integrate(&F, low, up, F.dim, par->calls, rng, s, &res, &err);
        iter++;
    } while (fabs(gsl_monte_vegas_chisq(s) - 1.0) > 0.25 && iter < 3);

    gsl_monte_vegas_free(s);
    gsl_rng_free(rng);

    return res;
}