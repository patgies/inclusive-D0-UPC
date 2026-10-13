"""
Inclusive D0 cross section (An0n) minus the diffractive one (0n0n), in the pT and y
bins of Fig. 9 of arXiv:2606.05469. With f(An0n) - f(Xn0n) = f(0n0n), this is Eq. (24)
of that paper: the cross section of the Xn0n selection.

The diffractive results (exclusive + diffractive parts) are the output of
diffractive-D0-UPC, copied to input/diffractive/0n0n/.

usage: python3 plotting_scripts/PbPb_bins_subtracted.py   (FLUX=Starlight for the STARlight flux)
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
from matplotlib.patches import Patch
from PbPb_bins import bin_averages, FRAGS, FIG9_PANELS, BAND_HANDLE, style_panel, draw_cms, finish

# FLUX=Starlight: inclusive runs of output/An0n_Starlight/ (FLUX_MODEL=TABLE) and diffractive runs with
# the STARlight flux, copied to input/diffractive/0n0n_Starlight/. Default: the EFF flux.
FLUX = os.environ.get("FLUX", "")
FLUX_NAME = {"Starlight": "STARlight"}.get(FLUX, FLUX)   # name shown in the figures
SUFFIX = f"_{FLUX}" if FLUX else ""
INCLUSIVE = f"output/An0n{SUFFIX}"
DIFFRACTIVE = f"input/diffractive/0n0n{SUFFIX}"

alphae = 1 / 137
mc = 1.5
e_c = 2 / 3
Nc = 3

# Dilute edge of the nucleus: the lead dipole files stop at b = 32 GeV^-1, where A T_A(b) sigma0/2 ~ 1.
# Beyond it the nuclear amplitude is A T_A(b) sigma0/2 times the proton one (arXiv:1309.6963), so that
# sigma_A^edge = EDGE_FACTOR * sigma_p, with sigma_p the proton-target result with the same flux,
# fragmentation function and scale (diffractive-D0-UPC output/proton_edge/, copied to
# input/diffractive/0n0n<SUFFIX>_proton_edge/). EDGE_FACTOR = (sigma0/2) int_{b>b_max} [A T_A]^2 d^2b,
# as nuclear_edge_factor(32) in diffractive-D0-UPC/plotting_scripts/D0.py. NUCLEAR_EDGE=0 switches it off.
NUCLEAR_EDGE = os.environ.get("NUCLEAR_EDGE", "1") == "1"
PROTON_EDGE = f"input/diffractive/0n0n{SUFFIX}_proton_edge"
EDGE_FACTOR = 10.21
SIGMA0_MB = 16.36   # sigma0/2 of the proton dipole fit


def proton_edge(process, frag, mu_f, mu_r):
    # {y: {pT: EDGE_FACTOR * proton cross section in mb/GeV}}, or {} if switched off or missing
    if not NUCLEAR_EDGE:
        return {}
    out = {}
    pattern = f"{PROTON_EDGE}/factor_{mu_f}/{frag}/D0_proton_baseline_{process}_{frag}_0n0n_y*.dat"
    for filename in sorted(glob.glob(pattern)):
        points = {}
        with open(filename) as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    pt, raw = (float(v) for v in line.split())
                    # diffractive_prefactor has no area: the proton gets sigma0/2 in mb (no GEVSQR_TO_MB)
                    points[pt] = EDGE_FACTOR * raw * diffractive_prefactor(process, pt, mu_r) * SIGMA0_MB \
                                 * (2 * math.pi) * pt
        out[read_rapidity(filename)] = points
    return out


# (mu_F, mu_R) factors of the scale band: 0.5 <= mu_F/mu_R <= 2
SCALE_FACTORS = [0.5, 1.0, 2.0]
SCALE_COMBOS = [(f, r) for f in SCALE_FACTORS for r in SCALE_FACTORS if 0.5 <= f / r <= 2.0]


def alphas_run(mu):
    # One-loop running coupling
    # Nf = 5 from mZ down to mb, Nf = 4 below, continuous at mb
    alphas_mZ, mZ, mb = 0.118, 91.2, 4.8
    b0 = lambda nf: (33 - 2 * nf) / (12 * math.pi)
    if mu >= mb:
        return 1.0 / (1.0 / alphas_mZ + 2 * b0(5) * math.log(mu / mZ))
    inv_mb = 1.0 / alphas_mZ + 2 * b0(5) * math.log(mb / mZ)
    return 1.0 / (inv_mb + 2 * b0(4) * math.log(mu / mb))


def flux_convention_factor(filename):
    # diffractive files written before the photon-flux
    # fix of Oct 2026 used dN/domega in the q+ integral and are sqrt(2) too large. New files have a
    # "# photon flux    : dN/dq+" header line.
    with open(filename) as f:
        for line in f:
            if not line.startswith("#"):
                break
            if line.startswith("# photon flux") and "dN/dq+" in line:
                return 1.0
    return 1.0 / math.sqrt(2.0)


def diffractive_prefactor(process, pt, mu_r):
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
        pattern = f"{scale_dir(DIFFRACTIVE, process, frag, mu_f)}/D0_{process}_{frag}_0n0n_Pb*_y*.dat"   # "_Pb_STARLIGHT_y" for the STARlight flux
        edge = proton_edge(process, frag, mu_f, mu_r)
        for filename in sorted(glob.glob(pattern)):
            y = read_rapidity(filename)
            flux_factor = flux_convention_factor(filename)
            b_list, pt_list, dsigma_list = read_data_file(filename)
            for pt, pairs in group_by_pt(b_list, pt_list, dsigma_list).items():
                value = (2 * math.pi) * integrate_over_b(pairs) * diffractive_prefactor(process, pt, mu_r) \
                        * (2 * math.pi) * pt * GEVSQR_TO_MB * flux_factor
                total.setdefault(y, {}).setdefault(pt, 0.0)
                total[y][pt] += value + edge.get(y, {}).get(pt, 0.0)
    return {y: sorted(points.items()) for y, points in total.items()}


def load_inclusive(frag, mu_f):
    return load_results(f"{scale_dir(INCLUSIVE, None, frag, mu_f)}/D0_incl_{frag}_An0n_Pb*_y*.dat")   # "_Pb_TABLE_y" for a table flux


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
    ("Xn0n", "#b2182b", "--"),   # An0n inclusive minus 0n0n diffractive
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
        style_panel(ax, pt_lo, pt_hi, -2.0, 2.0, top, legend_entries=4 if first else 0,
                    sublabel=(f"{label}\n{FLUX_NAME}" if FLUX else label) if first else None)

    handles = [Line2D([0], [0], color=color, linestyle=linestyle, lw=3, label=name)
               for name, color, linestyle in CURVES]
    # the band of the difference also has the mu_R variation of the diffractive part (alpha_s)
    handles.append(Patch(facecolor="0.3", alpha=0.25, edgecolor="none", label=r"$\mu_F, \mu_R \in [0.5, 2]\, m_T$"))
    handles.append(Line2D([0], [0], color="black", marker="o", linestyle="none", markersize=9, label="CMS"))
    finish(fig, axes, handles, filename, labelspacing=0.5)


def main():
    cms = load_cms_data()
    only = os.environ.get("FRAG_TYPE", "")
    if only and only not in FRAGS:
        sys.exit(f"Unknown FRAG_TYPE '{only}'. Expected one of: {', '.join(FRAGS)}.")
    for frag in ([only] if only else FRAGS):
        if not glob.glob(f"{DIFFRACTIVE}/central_values/{frag}/*.dat"):
            continue
        draw(frag, f"plots/D0_incl_corrected_{frag}{SUFFIX}.pdf", cms)


if __name__ == "__main__":
    main()
