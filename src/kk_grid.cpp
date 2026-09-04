#include "kk_grid.hpp"
#include <cmath>
#include <cstdio>
#include <vector>
#include <fcntl.h>
#include <unistd.h>
#include "QCDNUM/QCDNUM.h"

namespace {

// A few of QCDNum's Fortran WRITE statements are hardcoded to unit 6
// regardless of the lun passed to qcinit/setlun, so qcinit's own output
// redirection doesn't fully silence it. Redirect the OS-level stdout file
// descriptor instead -- this catches everything, since Fortran unit 6 is
// backed by fd 1 either way.
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

double funcC(int* ipdf, double* xx) {
  int    i = *ipdf;
  double x = *xx;
  double f = 0;
  if (i == 8) {                              // cbar column -> c=cbar=D_c
    double N = 0.694, eps = 0.101;
    double den = (1-x)*(1-x) + eps*x;
    f = N * x*x * (1-x)*(1-x) / (den*den);   // x * D_c(x,mc^2)
  }
  return f;
}

double funcB(int* ipdf, double* xx) {
  int    i = *ipdf;
  double x = *xx;
  double f = 0;
  if (i == 10) {                             // bbar column -> b=bbar=D_b
    double N = 81.7, alfa = 1.81, beta = 4.95;
    f = N * pow(x, alfa+1) * pow(1-x, beta); // x * D_b(x,mb^2)
  }
  return f;
}

const double QMC = 1.5, QMB = 5.0;
const double QMC2 = 2.25, QMB2 = 25.0;

// One-time QCDNum grid + evolution setup. Safe without extra locking: each
// ./dipole run is a separate OS process, never multiple threads sharing
// this state.
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

  // lun=-6 only suppresses the banner, not the citation box, READWT status
  // lines, or a handful of blank-line spacers that some QCDNum routines
  // write to a hardcoded unit 6 regardless of qcinit's lun. Those unprefixed
  // lines would corrupt run_many_Pb.sh's awk-based stdout parsing (it
  // expects the only non-'#' line to be the final result), so silence the
  // OS-level stdout fd for the whole setup instead of relying on qcinit.
  ScopedStdoutSilencer silence;

  QCDNUM::qcinit(20, "/dev/null");
  int nx, nq;
  QCDNUM::gxmake(xmin, iwt, ng, nxin, nx, iosp);
  QCDNUM::gqmake(qq, wt, 4, nqin, nq);
  QCDNUM::wtfile(3, "QCDnum/timelike.wgt");      // itype=3: time-like wgts
  QCDNUM::setord(1);                              // LO
  QCDNUM::setalf(as0, r20);
  int iqc = QCDNUM::iqfrmq(QMC2);
  int iqb = QCDNUM::iqfrmq(QMB2);
  QCDNUM::setcbt(0, iqc, iqb, 999);

  double epsC, epsB;
  int iq0C = QCDNUM::iqfrmq(QMC2);
  QCDNUM::evolfg(13, funcC, def, iq0C, epsC);     // itype=3(t-like),iset=1 -> jset=13

  int iq0B = QCDNUM::iqfrmq(QMB2);
  QCDNUM::evolfg(23, funcB, def, iq0B, epsB);     // itype=3(t-like),iset=2 -> jset=23
}

// Total D0 FF at (x,q2): sum of the charm-evolved (iset=1) and
// bottom-evolved (iset=2) sets, returned as D(x,Q2) (not x*D(x,Q2)).
double GetD(double x, double q2) {
  double pdfC[13], pdfB[13];
  QCDNUM::allfxq(1, x, q2, pdfC, 0, 1);
  double d = pdfC[10];                       // c-quark: Fortran index +4 -> C++ 10
  if (q2 >= QMB2) {
    QCDNUM::allfxq(2, x, q2, pdfB, 0, 1);
    d += pdfB[10];
  }
  return d / x;
}

}  // namespace

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
    // GetD/allfxq are undefined exactly at z=1 (1-z factors blow up);
    // nudge the last point in slightly.
    if (z >= 1.0) z = 1.0 - 1.e-6;
    zvals[i] = z;
    dvals[i] = GetD(z, q2);
  }

  return std::unique_ptr<Interpolator>(new Interpolator(zvals, dvals));
}
