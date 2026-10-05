"""
Inclusive D0 cross section (An0n) minus the diffractive one (0n0n), in the pT and y
bins of Fig. 9 of arXiv:2606.05469. With f(An0n) - f(Xn0n) = f(0n0n), this is Eq. (24)
of that paper: the cross section of the Xn0n selection.

The diffractive results (exclusive + diffractive parts) are the output of
diffractive-D0-UPC, copied to input/diffractive/0n0n/.

usage: python3 plotting_scripts/PbPb_bins_subtracted.py
"""
import glob
import math
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cross_section import (read_rapidity, read_data_file, group_by_pt, integrate_over_b,
                           load_results, GEVSQR_TO_MB)
from cms_comparison import load_cms_data   # also moves to the repo root and sets the plot style
from matplotlib.lines import Line2D
from PbPb_bins import bin_averages, FRAGS, FIG9_PANELS, BAND_HANDLE, style_panel, draw_cms, finish

INCLUSIVE = "output/An0n"
DIFFRACTIVE = "input/diffractive/0n0n"

alphae = 1 / 137
mc = 1.5
e_c = 2 / 3
Nc = 3

# (mu_F, mu_R) factors of the scale band: 0.5 <= mu_F/mu_R <= 2
SCALE_FACTORS = [0.5, 1.0, 2.0]
SCALE_COMBOS = [(f, r) for f in SCALE_FACTORS for r in SCALE_FACTORS if 0.5 <= f / r <= 2.0]


def alphas_run(mu):
    # One-loop running coupling, as in diffractive-D0-UPC/src/alphas_running.py
    alphas_mZ, mZ, Nf = 0.118, 91.2, 4
    b0 = (33 - 2 * Nf) / (12 * math.pi)
    return 1.0 / (1.0 / alphas_mZ + 2 * b0 * math.log(mu / mZ))


def diffractive_prefactor(process, pt, mu_r):
    # Prefactors of diffractive-D0-UPC/plotting_scripts/D0.py
    if process == "exclusive":
        return alphae * Nc * e_c**2 / (2 * math.pi**2)
    alphas = alphas_run(mu_r * math.sqrt(pt**2 + mc**2))
    return alphas * alphae * e_c**2 * (Nc**2 - 1) / (8 * math.pi**4)


def scale_dir(base, kind, frag, mu_f):
    if mu_f == 1.0:
        return f"{base}/central_values/{frag}"
    return f"{base}/scale_variation/{frag}/factor_{mu_f}"


def load_diffractive(frag, mu_f, mu_r):
    # exclusive + diffractive parts. Returns {y: [(pT, cross section), ...]} in mb/GeV.
    total = {}
    for process in ("exclusive", "diffractive"):
        pattern = f"{scale_dir(DIFFRACTIVE, process, frag, mu_f)}/D0_{process}_{frag}_0n0n_Pb_y*.dat"
        for filename in sorted(glob.glob(pattern)):
            y = read_rapidity(filename)
            b_list, pt_list, dsigma_list = read_data_file(filename)
            for pt, pairs in group_by_pt(b_list, pt_list, dsigma_list).items():
                value = (2 * math.pi) * integrate_over_b(pairs) * diffractive_prefactor(process, pt, mu_r) \
                        * (2 * math.pi) * pt * GEVSQR_TO_MB
                total.setdefault(y, {}).setdefault(pt, 0.0)
                total[y][pt] += value
    return {y: sorted(points.items()) for y, points in total.items()}


def load_inclusive(frag, mu_f):
    return load_results(f"{scale_dir(INCLUSIVE, None, frag, mu_f)}/D0_incl_{frag}_An0n_Pb_y*.dat")


def bands(frag, pt_lo, pt_hi, y_edges):
    # (central, min, max) of the inclusive, the diffractive and their difference
    incl, diff, sub = [], [], []
    for mu_f, mu_r in [(1.0, 1.0)] + [c for c in SCALE_COMBOS if c != (1.0, 1.0)]:
        i = bin_averages(load_inclusive(frag, mu_f), pt_lo, pt_hi, y_edges)
        d = bin_averages(load_diffractive(frag, mu_f, mu_r), pt_lo, pt_hi, y_edges)
        incl.append(i)
        diff.append(d)
        sub.append(i - d)
    return [(np.array(v)[0], np.min(v, axis=0), np.max(v, axis=0)) for v in (incl, diff, sub)]


# (label, color, line style) of the two curves
CURVES = [
    ("An0n", "#2166ac", "-"),    # inclusive
    ("An0n $-$ 0n0n (diff)", "#b2182b", "--"),   # An0n inclusive minus 0n0n diffractive
]


def draw(frag, filename, cms):
    label = FRAGS[frag][0]
    fig, axes = plt.subplots(2, 2, figsize=(16, 11.34), sharex=True)
    print(f"{label}: diffractive (0n0n) / inclusive (An0n), central scale")

    for ax, (pt_lo, pt_hi, y_edges, show_cms) in zip(axes.flat, FIG9_PANELS):
        incl, diff, sub = bands(frag, pt_lo, pt_hi, y_edges)
        for (central, low, high), (_, color, linestyle) in zip((incl, sub), CURVES):
            # one horizontal line per y bin, without the vertical lines joining the bins
            ax.hlines(central, y_edges[:-1], y_edges[1:], color=color, linestyle=linestyle, lw=3)
            for y_lo, y_hi, lo, hi in zip(y_edges[:-1], y_edges[1:], low, high):
                ax.fill_between([y_lo, y_hi], lo, hi, color=color, alpha=0.25, linewidth=0)
        top = np.nanmax(incl[2])
        ratios = "  ".join(f"{r:.3f}" for r in diff[0] / incl[0])
        print(f"   {pt_lo:g} < pT < {pt_hi:g} GeV, y bins {y_edges}: {ratios}")
        if show_cms and (pt_lo, pt_hi) in cms:
            top = max(top, draw_cms(ax, cms[(pt_lo, pt_hi)]))
        first = ax is axes.flat[0]
        style_panel(ax, pt_lo, pt_hi, -2.0, 2.0, top, legend_entries=5 if first else 0,
                    sublabel=label if first else None)

    handles = [Line2D([0], [0], color=color, linestyle=linestyle, lw=3, label=name)
               for name, color, linestyle in CURVES]
    handles.append(BAND_HANDLE)
    handles.append(Line2D([0], [0], color="black", marker="o", linestyle="none", markersize=9, label="CMS"))
    finish(fig, axes, handles, filename, legend_title="Inclusive", labelspacing=0.4)


def main():
    cms = load_cms_data()
    for frag in FRAGS:
        if not glob.glob(f"{DIFFRACTIVE}/central_values/{frag}/*.dat"):
            continue
        draw(frag, f"plots/D0_incl_corrected_{frag}.pdf", cms)


if __name__ == "__main__":
    main()
