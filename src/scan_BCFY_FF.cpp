// Standalone diagnostic: dump the DGLAP-evolved Braaten-Cheung-Fleming-Yuan
// (BCFY) c->D0 fragmentation function D(z, Q^2) at several mu values, for a
// couple of reference pT points, to visualize how the shape shifts with the
// evolution scale (see python/plot_kk_ff_shapes.py for the KK analogue).
//
// usage: ./build/bin/scan_BCFY_FF > out/bcfy_ff_shapes/bcfy_ff.csv
#include "bcfy_grid.hpp"
#include "fragmentation.hpp"
#include "interpolation.hpp"
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

namespace {
const double R = 0.1;

double BCFY_D0_IC(double z) { return Dc_to_D0(z, R); }

void DumpCurve(double pT, const std::string& label, double (*f)(double)) {
    for (int i = 0; i < 200; i++) {
        double z = 0.05 + (1.0 - 1.e-6 - 0.05) * i / 199.0;
        std::cout << pT << "," << label << "," << z << "," << f(z) << "\n";
    }
}
}

int main() {
    double mc = 1.5;
    std::vector<double> pTvals = {3.0, 9.0};
    std::vector<double> scale_factors = {0.25, 0.5, 1.0, 2.0, 4.0};

    std::cout << "pT,label,z,D\n";

    for (double pT : pTvals) {
        double mt = std::sqrt(mc*mc + pT*pT);

        // Raw (non-evolved) c->D0 input, i.e. 0.168*BCFY_DP + feed-down
        // BCFY_DV, evaluated directly at mu=mc -- same curve regardless of
        // pT, but repeated per-pT so each panel is self-contained.
        DumpCurve(pT, "IC (μ=mc)", BCFY_D0_IC);

        for (double scale_factor : scale_factors) {
            double Q = scale_factor * mt;
            auto interp = MakeBCFYInterpolator(Q);
            std::ostringstream label;
            label.precision(4);
            label << "μ=" << scale_factor << "*mt=" << Q;
            for (int i = 0; i < 200; i++) {
                double z = 0.05 + (1.0 - 1.e-6 - 0.05) * i / 199.0;
                double d = interp->Evaluate(z);
                std::cout << pT << "," << label.str() << "," << z << "," << d << "\n";
            }
        }
    }
    return 0;
}
