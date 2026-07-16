#ifndef def_hpp 
#define def_hpp
#include <memory>
#include <string>
#include "amplitudelib.hpp"
#include "interpolation.hpp"

enum class FragmentationType { BCFY, KniehlKramer, LHAPDF };

struct parameters
{
    AmplitudeLib *dipole;

    //Kinematics (charm/D0) 
    double pD0;     // target D0 transverse momentum (the produced hadron, not the charm quark)
    double m, m2;   // charm mass, mass^2
    double y;       // rapidity
    double ss;      // sqrt(s)
    double xbj;     // representative xbj (at z=1), used only to seed the interpolator cache
    double qpmin;   // tightest valid qp lower bound, reached at z=1 (pc=pD0)

    // Fragmentation (c -> D0)
    double r;            // BCFY non-perturbative parameter
    double N_kk, eps_kk; // Kniehl & Kramer parameters (N=0.694, eps=0.101)
    FragmentationType frag_type = FragmentationType::BCFY;
    std::string lhapdf_setname; // LHAPDF FF set name, e.g. "prompt-D0-1-109"
    int lhapdf_pid;             // parton id fragmenting into the D0 (4 = charm)
    double lhapdf_Q;            // fragmentation scale Q passed to xfxQ (GeV)
    double zmin, zmax;   // fragmentation z integration range

    //Photon flux
    double alpha, Z, mn, S;
    std::string channel;

    // VEGAS integration box 
    double bmin, bmax, qpmax, rmax, lmax;
    size_t calls;

    // Precomputed S_k grid for momentum-space evaluation
    std::unique_ptr<Interpolator> Sk_interp;

    // Fixed photon momentum for differential-in-qp calculation
    double qp_fixed = 0.0;
};

// Assumes par->Sk_interp is already set (precomputed in main before parallel launch).
double D0CrossSection_inclusive(void* p);

#endif
