"""
Nuclear modification factor of Eq. (40) of arXiv:2503.16108 (Fig. 5 there):

    R_pA = dsigma_AA / (A * dsigma_pA),

both differential in y and pT. The proton baseline dsigma_pA uses the same photon flux as
Pb+Pb, so R_pA is not measurable but shows the size of the saturation effects (R_pA -> 1 at
large pT).

The Glauber files only go up to b_d = 32 GeV^-1, so the edge of the nucleus is missing. As in
the paper, it is added to the numerator as A_edge * dsigma_pA, with A_edge the number of
nucleons of the Woods-Saxon nucleus beyond the last b_d (EDGE=0 leaves it out).

Numerator:   output/<CHANNEL>/central_values/<FRAG>/            (run_scripts/run_nucleus.sh)
Denominator: output/<CHANNEL>/proton_baseline/<FRAG>/           (run_scripts/run_proton.sh, TARGET=AA)

usage: FRAG_TYPE=KniehlKramer python3 plotting_scripts/RpA.py
"""
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cross_section import load_results, load_proton_results, add_edge
import cms_comparison   # moves to the repo root and sets the plot style

CHANNEL = os.environ.get("CHANNEL", "An0n")
FRAG = os.environ.get("FRAG_TYPE", "KniehlKramer")   # the paper uses Kniehl-Kramer
PT_MAX = float(os.environ.get("PT_MAX", 5.0))        
A = 208
Y_TO_PLOT = [0.0, 0.5, 1.0, 1.5, 2.0]
# one hue, dark (y = 0) to light (y = 2)
COLORS = ["#08306b", "#08519c", "#2171b5", "#4292c6", "#6baed6"]
LINESTYLES = ["-", (0, (8, 2)), "--", "-.", ":"]


def r_pa():
    # {y: (pT values, R_pA)} for the rapidities that are in both runs
    nucleus_files = f"output/{CHANNEL}/central_values/{FRAG}/D0_incl_{FRAG}_{CHANNEL}_Pb_y*.dat"
    proton_files = f"output/{CHANNEL}/proton_baseline/{FRAG}/D0_proton_baseline_incl_{FRAG}_{CHANNEL}_y*.dat"
    proton = load_proton_results(proton_files)
    # the nucleus with its edge, on the rapidities of the proton run
    nucleus = {y: v for y, v in load_results(nucleus_files).items() if y in proton}
    nucleus = add_edge(nucleus, nucleus_files, proton_files)
    out = {}
    for y in Y_TO_PLOT:
        if y not in nucleus or y not in proton:
            continue
        aa, pa = dict(nucleus[y]), dict(proton[y])
        pts = sorted(pt for pt in set(aa) & set(pa) if pt <= PT_MAX + 1e-9)
        out[y] = (np.array(pts), np.array([aa[pt] / (A * pa[pt]) for pt in pts]))
    return out


def main():
    results = r_pa()
    if not results:
        sys.exit(f"No overlapping rapidities in output/{CHANNEL}/central_values/{FRAG}/ and "
                 f"output/{CHANNEL}/proton_baseline/{FRAG}/ -- run run_scripts/run_nucleus.sh and "
                 f"run_scripts/run_proton.sh (TARGET=AA) first.")

    os.makedirs("plots", exist_ok=True)
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    handles = []
    for i, (y, (pts, values)) in enumerate(results.items()):
        ax.plot(pts, values, color=COLORS[i], linestyle=LINESTYLES[i], lw=2.5)
        handles.append(Line2D([0], [0], color=COLORS[i], linestyle=LINESTYLES[i], lw=2.5, label=f"$y = {y:g}$"))
    ax.axhline(1.0, color="gray", linestyle=":", linewidth=1)
    ax.set_xlim(0, PT_MAX)
    ax.set_xlabel(r"$p_{D^0\perp}$ [GeV]", labelpad=12)
    ax.set_ylabel(r"$R_{pA}$", labelpad=10)
    ax.text(0.05, 0.93, rf"Pb + Pb, {CHANNEL}, {FRAG}", transform=ax.transAxes, ha="left", va="top", fontsize=16)
    ax.legend(handles=handles, loc="lower right", fontsize=16, frameon=False)
    fig.tight_layout()
    outname = f"plots/RpA_{CHANNEL}_{FRAG}.pdf"
    fig.savefig(outname, bbox_inches="tight")
    print(f"Saved: {outname}")


if __name__ == "__main__":
    main()
