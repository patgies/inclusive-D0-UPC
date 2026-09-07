#include "bcfy_grid.hpp"
#include "fragmentation.hpp"
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

const double R = 0.1;  // light-quark/meson mass ratio (D mesons), same value
                        // used for the raw BCFY_DP/BCFY_DV in fragmentation.cpp

double funcP(int* ipdf, double* xx) {
  int    i = *ipdf;
  double x = *xx;
  double f = 0;
  if (i == 8) f = x * BCFY_DP(x, R);   // cbar column -> c=cbar seed
  return f;
}

double funcV(int* ipdf, double* xx) {
  int    i = *ipdf;
  double x = *xx;
  double f = 0;
  if (i == 8) f = x * BCFY_DV(x, R);   // same column: also charm-seeded
  return f;
}

double theta_step(double x) { return (x >= 0.0) ? 1.0 : 0.0; }

const double QMC = 1.5, QMC2 = 2.25;
const double mD = 1.8648, mDstar = 2.0067;
const double mass_ratio = mDstar / mD;

// One-time QCDNum grid + evolution setup for both the P and V channels.
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

  double as0 = 0.118, r20 = 8315.25;              // alphas(Mz2)=0.118
  double xmin[] = {1.e-5};
  int    iwt[] = {1}, ng = 1, nxin = 300, iosp = 3;
  int    nqin = 100;
  double qq[] = {1.0, 2.25, 25.0, 1.e6}, wt[] = {1.0, 1.0, 2.0, 1.0};

  ScopedStdoutSilencer silence;

  QCDNUM::qcinit(20, "/dev/null");
  int nx, nq;
  QCDNUM::gxmake(xmin, iwt, ng, nxin, nx, iosp);
  QCDNUM::gqmake(qq, wt, 4, nqin, nq);

  // Per-process scratch weight file -- a shared path caused concurrent
  // writers to corrupt each other's file when many (b, pT, y) points run
  // concurrently (see kk_grid.cpp's EnsureEvolved for the same fix).
  std::string wtpath = "/tmp/qcdnum_timelike_bcfy_" + std::to_string(getpid()) + ".wgt";
  QCDNUM::wtfile(3, wtpath);      // itype=3: time-like wgts
  remove(wtpath.c_str());
  QCDNUM::setord(1);                              // LO
  QCDNUM::setalf(as0, r20);
  int iqc = QCDNUM::iqfrmq(QMC2);
  QCDNUM::setcbt(0, iqc, 999, 999);               // no b, no t threshold

  double epsP;
  int iq0 = QCDNUM::iqfrmq(QMC2);
  QCDNUM::evolfg(13, funcP, def, iq0, epsP);      // itype=3(t-like),iset=1 -> jset=13
  double epsV;
  QCDNUM::evolfg(23, funcV, def, iq0, epsV);      // itype=3(t-like),iset=2 -> jset=23
}

// D0 FF at (x,q2): 0.168*P-evolved + 0.39*feeddown*V-evolved, combined AFTER
// evolution -- P and V are evolved as two independent DGLAP sets, and the
// hadron-level D*->D0 feed-down is applied on top (see testjobs/bcfyD0.cc in
// QCDnumFF for the derivation). Returned as D(x,Q2), not x*D(x,Q2).
double GetD(double x, double q2) {
  double pdfP[13];
  QCDNUM::allfxq(1, x, q2, pdfP, 0, 1);
  double xDtot = 0.168 * pdfP[10];              // c-quark: Fortran index +4 -> C++ 10

  double x_scaled = mass_ratio * x;
  if (theta_step(mD/mDstar - x) > 0.0 && x_scaled < 1.0) {
    double pdfV[13];
    QCDNUM::allfxq(2, x_scaled, q2, pdfV, 0, 1);
    xDtot += 0.39 * pdfV[10];
  }
  return xDtot / x;
}

}  // namespace

std::unique_ptr<Interpolator> MakeBCFYInterpolator(double Q) {
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
