"""
Inclusive D0 cross section in p+Pb UPCs (TARGET=pA), in bins of pT as a function of y,
as in Fig. 12 of arXiv:2606.05469. The lead nucleus moves towards positive y.
The bands are the fragmentation scale variation (Q between 0.5 and 2 times mT).

usage: python3 plotting_scripts/pPb_bins.py
Reads output/pPb/central_values/<FRAG>/ and output/pPb/scale_variation/<FRAG>/factor_<0.5|2.0>/.
"""
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cross_section import read_rapidity, factor_A, pi
import cms_comparison   # moves to the repo root and sets the plot style
from PbPb_bins import draw, FRAGS, SCALE_FACTORS, NARROW_PANELS

OUTPUT = "output/pPb"
sigma0 = 16.36   # mb, proton normalization of the MVe dipole. It replaces the b_d integral.


def load_pPb(pattern):
    # Reads the files "pD0  dsigma" of run_proton.sh. Returns {y: [(pT, cross section), ...]} in mb/GeV.
    results = {}
    for filename in sorted(glob.glob(pattern)):
        points = []
        with open(filename) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                pt, raw = (float(x) for x in line.split())
                points.append((pt, raw * factor_A * sigma0 * (2 * pi) * pt))
        results[read_rapidity(filename)] = points
    return results


def load_frag(frag):
    name = f"D0_pA_incl_{frag}_y*.dat"
    runs = [load_pPb(f"{OUTPUT}/central_values/{frag}/{name}")]
    if not runs[0]:
        return None
    for factor in SCALE_FACTORS:
        varied = load_pPb(f"{OUTPUT}/scale_variation/{frag}/factor_{factor}/{name}")
        if varied:
            runs.append(varied)
    return runs


def main():
    frag_runs = {frag: runs for frag in FRAGS for runs in [load_frag(frag)] if runs is not None}
    if not frag_runs:
        sys.exit(f"No results in {OUTPUT}/central_values/ -- run TARGET=pA run_scripts/run_proton.sh first.")
    draw(NARROW_PANELS, "plots/D0_incl_bins_pPb.pdf", frag_runs=frag_runs, sublabel="p+Pb")


if __name__ == "__main__":
    main()
