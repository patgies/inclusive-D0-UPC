#include "lhapdf_grid.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {

std::vector<double> ParseDoubles(const std::string& line)
{
    std::vector<double> out;
    std::istringstream iss(line);
    double v;
    while (iss >> v) out.push_back(v);
    return out;
}

std::vector<int> ParseInts(const std::string& line)
{
    std::vector<int> out;
    std::istringstream iss(line);
    int v;
    while (iss >> v) out.push_back(v);
    return out;
}

} // namespace

std::unique_ptr<Interpolator> MakeLHAPDFZInterpolator(
    const std::string& member_file, int pdg_flavor, double Q)
{
    std::ifstream in(member_file);
    if (!in) throw std::runtime_error("MakeLHAPDFZInterpolator: cannot open " + member_file);

    std::string line;
    // Skip header lines (PdfType:, Format:, ...) up to the first "---" separator.
    while (std::getline(in, line)) {
        if (line.compare(0, 3, "---") == 0) break;
    }

    // lhagrid1 format: repeated subgrid blocks, each xgrid/qgrid/flavors
    // header followed by nx*nq value rows (x outer loop, Q inner loop),
    // terminated by a "---" line. Subgrids are split at flavor thresholds
    // and overlap exactly at the boundary, so the first subgrid whose Q
    // range contains the requested Q is used.
    while (std::getline(in, line)) {
        std::vector<double> xgrid = ParseDoubles(line);
        if (xgrid.empty()) break;

        if (!std::getline(in, line))
            throw std::runtime_error("MakeLHAPDFZInterpolator: truncated Q-grid line in " + member_file);
        std::vector<double> qgrid = ParseDoubles(line);

        if (!std::getline(in, line))
            throw std::runtime_error("MakeLHAPDFZInterpolator: truncated flavor line in " + member_file);
        std::vector<int> flavors = ParseInts(line);

        int flavor_col = -1;
        for (size_t i = 0; i < flavors.size(); i++) {
            if (flavors[i] == pdg_flavor) { flavor_col = (int)i; break; }
        }
        if (flavor_col < 0)
            throw std::runtime_error("MakeLHAPDFZInterpolator: flavor id not found in grid " + member_file);

        size_t nx = xgrid.size();
        size_t nq = qgrid.size();
        std::vector<std::vector<double>> values(nx, std::vector<double>(nq));

        for (size_t ix = 0; ix < nx; ix++) {
            for (size_t iq = 0; iq < nq; iq++) {
                if (!std::getline(in, line))
                    throw std::runtime_error("MakeLHAPDFZInterpolator: truncated value grid in " + member_file);
                std::vector<double> row = ParseDoubles(line);
                if ((int)row.size() <= flavor_col)
                    throw std::runtime_error("MakeLHAPDFZInterpolator: short value row in " + member_file);
                values[ix][iq] = row[flavor_col];
            }
        }

        std::getline(in, line);  // consume the "---" line ending this subgrid block

        if (Q < qgrid.front() || Q > qgrid.back()) continue;  // not this subgrid

        // Bracket Q within this subgrid's Q points, then linearly interpolate
        // in Q for every x point to build a 1D f(x, Q) slice.
        size_t iq_hi = 1;
        while (iq_hi < qgrid.size() - 1 && qgrid[iq_hi] < Q) iq_hi++;
        size_t iq_lo = iq_hi - 1;
        double q_lo = qgrid[iq_lo], q_hi = qgrid[iq_hi];
        double t = (q_hi > q_lo) ? (Q - q_lo) / (q_hi - q_lo) : 0.0;

        std::vector<double> zvals(nx), dvals(nx);
        for (size_t ix = 0; ix < nx; ix++) {
            double xf = values[ix][iq_lo] * (1.0 - t) + values[ix][iq_hi] * t;
            zvals[ix] = xgrid[ix];
            dvals[ix] = xf / xgrid[ix];   // grid stores x*f(x,Q) -> convert to D(z,Q)
        }

        return std::unique_ptr<Interpolator>(new Interpolator(zvals, dvals));
    }

    throw std::runtime_error("MakeLHAPDFZInterpolator: Q outside all subgrids in " + member_file);
}
