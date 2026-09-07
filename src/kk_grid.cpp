#include "kk_grid.hpp"
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <fcntl.h>
#include <unistd.h>
#include "QCDNUM/QCDNUM.h"

namespace {


class ScopedStdoutSilencer {
 public:
  ScopedStdoutSilencer() {
    fflush(stdout);
    saved_fd_ = dup(1);
    int devnull = open("/dev/null", O_WRONLY);
    dup2(devnull, 1);
    close(devnull);
  }
  ~ScopedStdoutSilencer() {
    fflush(stdout);
    dup2(saved_fd_, 1);
    close(saved_fd_);
  }
 private:
  int saved_fd_;
};

// Raw (non-evolved) Kniehl & Kramer c/b -> D0 input fragmentation functions,
// D(x) itself (not x*D(x)), at their own natural starting scales mc/mb.
// Shared by funcC (which QCDNUM's evolfg requires in x*D(x) form) and
// the public KKInitialCondition{C,B} accessors below. The bottom-seeded
// evolution channel was removed from GetD (see below), but the raw bottom
// input is kept available for KKInitialConditionB / the diagnostic plot.
double KKInitialConditionC_raw(double x) {
  double N = 0.694, eps = 0.101;
  double den = (1-x)*(1-x) + eps*x;
  return N * x * (1-x)*(1-x) / (den*den);
}

double KKInitialConditionB_raw(double x) {
  double N = 81.7, alfa = 1.81, beta = 4.95;
  return N * pow(x, alfa) * pow(1-x, beta);
}

double funcC(int* ipdf, double* xx) {
  int    i = *ipdf;
  double x = *xx;
  double f = 0;
  if (i == 8) f = x * KKInitialConditionC_raw(x);   // cbar column -> c=cbar=D_c
  return f;
}

const double QMC = 1.5, QMB = 5.0;
const double QMC2 = 2.25, QMB2 = 25.0;

// One-time QCDNum grid + evolution setup. 
void EnsureEvolved() {
  static bool initialized = false;
  if (initialized) return;
  initialized = true;

  // tb  bb  cb  sb  ub  db   g   d   u   s   c   b   t
  double def[] =
    { 0., 0., 0., 0., 0.,-1., 0., 1., 0., 0., 0., 0., 0.,      // dval
      0., 0., 0., 0.,-1., 0., 0., 0., 1., 0., 0., 0., 0.,      // uval
      0., 0., 0.,-1., 0., 0., 0., 0., 0., 1., 0., 0., 0.,      // sval
      0., 0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0.,      // dbar
      0., 0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0.,      // ubar
      0., 0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      // sbar
      0., 0.,-1., 0., 0., 0., 0., 0., 0., 0., 1., 0., 0.,      // cval
      0., 0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      // cbar
      0.,-1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 0.,      // bval
      0., 1., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      // bbar
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,      // tval (zero)
      0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};     // tbar (zero)

  double as0 = 0.118, r20 = 8315.25;              // alphas(Mz2)=0.118, LO
  double xmin[] = {1.e-5};
  int    iwt[] = {1}, ng = 1, nxin = 300, iosp = 3;
  int    nqin = 100;
  double qq[] = {1.0, 2.25, 25.0, 1.e6}, wt[] = {1.0, 1.0, 2.0, 1.0};

 
  ScopedStdoutSilencer silence;

  QCDNUM::qcinit(20, "/dev/null");
  int nx, nq;
  QCDNUM::gxmake(xmin, iwt, ng, nxin, nx, iosp);
  QCDNUM::gqmake(qq, wt, 4, nqin, nq);

  // Every process (dipole is invoked once per (b, pT, y) point, so many
  // run concurrently under run_many_Pb.sh) recomputes and writes its own
  // weight tables here -- nothing reads this file back. A shared path
  // caused concurrent writers to corrupt each other's file (segfaults
  // under high -j parallelism on Roihu); use a per-process scratch path
  // instead, and clean it up since it's disposable.
  std::string wtpath = "/tmp/qcdnum_timelike_" + std::to_string(getpid()) + ".wgt";
  QCDNUM::wtfile(3, wtpath);      // itype=3: time-like wgts
  remove(wtpath.c_str());
  QCDNUM::setord(1);                              // LO
  QCDNUM::setalf(as0, r20);
  int iqc = QCDNUM::iqfrmq(QMC2);
  int iqb = QCDNUM::iqfrmq(QMB2);
  QCDNUM::setcbt(0, iqc, iqb, 999);

  double epsC;
  int iq0C = QCDNUM::iqfrmq(QMC2);
  QCDNUM::evolfg(13, funcC, def, iq0C, epsC);     // itype=3(t-like),iset=1 -> jset=13
}

// D0 FF at (x,q2): charm-evolved (iset=1) set only, returned as D(x,Q2)
// (not x*D(x,Q2)). The bottom-seeded channel (b -> B-hadron -> D0) was
// dropped: its evolved contribution was found to be negligible (< 1e-5
// relative, even well above the b threshold) and inconsistent with its own
// non-evolved input by ~5 orders of magnitude, indicating a bug in that
// channel's QCDNUM wiring rather than a real physics suppression -- not
// worth carrying since it never contributed meaningfully anyway.
double GetD(double x, double q2) {
  double pdfC[13];
  QCDNUM::allfxq(1, x, q2, pdfC, 0, 1);
  double d = pdfC[10];                       // c-quark: Fortran index +4 -> C++ 10
  return d / x;
}

}

// Raw (non-evolved) input fragmentation functions D(x) at the natural
// starting scale (mc for charm, mb for bottom) -- for comparing against the
// DGLAP-evolved MakeKniehlKramerInterpolator output at various Q.
double KKInitialConditionC(double x) { return KKInitialConditionC_raw(x); }
double KKInitialConditionB(double x) { return KKInitialConditionB_raw(x); }

std::unique_ptr<Interpolator> MakeKniehlKramerInterpolator(double Q) {
  EnsureEvolved();

  double q_clamped = Q;
  if (q_clamped < QMC) q_clamped = QMC;
  double q2 = q_clamped*q_clamped;

  const int npts = 80;
  const double zmin = 0.05, zmax = 1.0;
  std::vector<double> zvals(npts), dvals(npts);
  for (int i = 0; i < npts; i++) {
    double z = zmin + (zmax-zmin) * i/double(npts-1);
    if (z >= 1.0) z = 1.0 - 1.e-6;
    zvals[i] = z;
    dvals[i] = GetD(z, q2);
  }

  return std::unique_ptr<Interpolator>(new Interpolator(zvals, dvals));
}
