#ifndef def_hpp 
#define def_hpp
#include <memory>
#include <string>
#include "amplitudelib.hpp"
#include "interpolation.hpp"

enum class FragmentationType { BCFY, KniehlKramer, HymnD };

struct parameters
{
    AmplitudeLib *dipole;

    // Kinematics
    double pD0;     // D0 transverse momentum
    double m, m2;   // charm mass and its square
    double y;       // D0 rapidity
    double ss;      // sqrt(s_NN)
    double xbj;     // x of the dipole amplitude

    // Fragmentation (c -> D0)
    double r;            // BCFY parameter
    double N_kk, eps_kk; // Kniehl-Kramer parameters
    FragmentationType frag_type = FragmentationType::BCFY;
    double z_h_min, z_h_max;   // range of the z_h integral (pc = pD0/z_h)

    // Photon flux
    double alpha, Z, mn, S;   // alpha_em, charge, nucleon mass, EMD area
    std::string target = "AA";      // AA (Pb+Pb) | pA (p+Pb: Gamma_pA, no EMD factor)
    std::string channel;
    bool gamma_aa_one = false;        // true: Gamma_AA(b) = 1
    std::string flux_model = "EFF";   // EFF | PL | WS | TABLE

    // VEGAS integration limits
    double bmin, bmax, qpmax, lmax;
    size_t calls;

    // Dipole in momentum space, S(l)
    std::unique_ptr<Interpolator> Sk_interp;

    // Fragmentation function D(z_h) at the fragmentation scale
    std::unique_ptr<Interpolator> D_frag_interp;
};

// dsigma/(dy d2pD0) without prefactors. Needs par->Sk_interp and par->D_frag_interp.
double D0CrossSection_inclusive(void* p);

#endif
