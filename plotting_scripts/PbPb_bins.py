"""
Inclusive D0 cross section in Pb+Pb UPCs, in bins of pT as a function of y,
as in Fig. 9 of arXiv:2606.05469. The bands are the fragmentation scale
variation (Q between 0.5 and 2 times mT).

usage: python3 plotting_scripts/PbPb_bins.py
Reads output/<CHANNEL>/central_values/<FRAG>/ and
output/<CHANNEL>/scale_variation/<FRAG>/factor_<0.5|2.0>/.
"""
import glob
import os
import re
import sys

import numpy as np
from scipy.integrate import simpson
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator, ScalarFormatter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cross_section import load_results
from cms_comparison import load_cms_data   # also moves to the repo root and sets the plot style

CHANNEL = os.environ.get("CHANNEL", "An0n")
# FLUX=Starlight reads output/<CHANNEL>_Starlight/ (runs with FLUX_MODEL=TABLE). Default: the EFF flux.
FLUX = os.environ.get("FLUX", "")
OUTPUT = f"output/{CHANNEL}_{FLUX}" if FLUX else f"output/{CHANNEL}"
SUFFIX = f"_{FLUX}" if FLUX else ""
FLUX_LABEL = f"{FLUX} flux" if FLUX else None

# fragmentation function: (label, color, line style)
FRAGS = {
    "BCFY": ("BCFY", "#2166ac", "-"),
    "KniehlKramer": ("Kniehl-Kramer", "#1b7837", "--"),
    "HymnD": ("HymnD", "#b2182b", ":"),
}
SCALE_FACTORS = ["0.5", "2.0"]
# legend entry of the scale-variation bands (as in diffractive-D0-UPC)
BAND_HANDLE = Patch(facecolor="0.3", alpha=0.25, edgecolor="none", label=r"$\mu_F \in [0.5, 2]\, m_t$")

# Panels (pT bin, y bin edges, draw the CMS data) of the two figures.
# Fig. 9 of arXiv:2606.05469: the CMS bins. CMS has one y bin for 2 < pT < 5 GeV.
Y_EDGES = [-2.0, -1.0, 0.0, 1.0, 2.0]
FIG9_PANELS = [
    (2.0, 5.0, [-1.0, 1.0], True),
    (2.0, 5.0, Y_EDGES, False),
    (5.0, 8.0, Y_EDGES, True),
    (8.0, 12.0, Y_EDGES, True),
]
# Fixed upper limit of the vertical axis of a pT bin. The others are automatic.
Y_AXIS_MAX = {}   # e.g. {(2.0, 5.0): 8.0}
# Narrow pT bins, y from -3 to 3.
WIDE_Y_EDGES = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
NARROW_PANELS = [(lo, hi, WIDE_Y_EDGES, False)
               for lo, hi in [(0.0, 1.0), (1.0, 2.0), (2.0, 3.0), (3.0, 4.0), (6.0, 7.0), (10.0, 11.0)]]


def pt_average(points, pt_lo, pt_hi):
    # Average of d2sigma/dy dpT over the pT bin. The curve goes to zero at pT = 0.
    pt = np.array([0.0] + [p[0] for p in sorted(points)])
    cs = np.array([0.0] + [p[1] for p in sorted(points)])
    inside = (pt > pt_lo) & (pt < pt_hi)
    x = np.concatenate(([pt_lo], pt[inside], [pt_hi]))
    y = np.concatenate(([np.interp(pt_lo, pt, cs)], cs[inside], [np.interp(pt_hi, pt, cs)]))
    return simpson(y, x=x) / (pt_hi - pt_lo)


def bin_averages(results, pt_lo, pt_hi, y_edges):
    # Average over the pT bin and over each y bin. A y bin with missing y values gives nan.
    per_y = {y: pt_average(results[y], pt_lo, pt_hi) for y in results}
    values = []
    for y_lo, y_hi in zip(y_edges[:-1], y_edges[1:]):
        ys = sorted(y for y in per_y if y_lo - 1e-9 <= y <= y_hi + 1e-9)
        if len(ys) < 3 or abs(ys[0] - y_lo) > 1e-9 or abs(ys[-1] - y_hi) > 1e-9:
            values.append(np.nan)
            continue
        values.append(simpson([per_y[y] for y in ys], x=ys) / (y_hi - y_lo))
    return np.array(values)


def load_frag(frag):
    # results of the central scale and of the scale variation, for one fragmentation function
    name = f"D0_incl_{frag}_{CHANNEL}_Pb*_y*.dat"   # "_Pb_TABLE_y" for a table flux
    runs = [load_results(f"{OUTPUT}/central_values/{frag}/{name}")]
    if not runs[0]:
        return None
    for factor in SCALE_FACTORS:
        varied = load_results(f"{OUTPUT}/scale_variation/{frag}/factor_{factor}/{name}")
        if varied:
            runs.append(varied)
    return runs


def band(runs, pt_lo, pt_hi, y_edges):
    # central value, minimum and maximum of the runs in each y bin
    stack = np.array([bin_averages(run, pt_lo, pt_hi, y_edges) for run in runs])
    return stack[0], stack.min(axis=0), stack.max(axis=0)


class TrimmedFormatter(ScalarFormatter):
    # Tick labels without zeros at the end of a decimal: 0.01 instead of 0.010.
    def _set_format(self):
        super()._set_format()
        self.format = re.sub(r"%1\.\d+f", "%g", self.format)


def style_panel(ax, pt_lo, pt_hi, y_min, y_max, top, legend_entries=0, sublabel=None):
    # Axes, ticks and pT label of one panel (style of diffractive-D0-UPC).
    ax.set_xlim(y_min, y_max)
    # more room above the data in the panel with the legend
    # the legend takes about 0.12 of the panel height per row, from 0.91 down
    automatic = (1.05 / (0.89 - 0.12 * legend_entries) if legend_entries else 1.25) * top
    ax.set_ylim(0, Y_AXIS_MAX.get((pt_lo, pt_hi), automatic))
    ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.tick_params(labelsize=30, pad=10)
    # grey tick marks, black tick labels
    ax.tick_params(which="major", length=9, width=1.2, color="0.4", labelcolor="black")
    ax.tick_params(which="minor", length=4.5, width=1.0, color="0.4")
    formatter = TrimmedFormatter(useMathText=True)
    formatter.set_powerlimits((-2, 2))
    # the x10^n label in black, a little above the frame
    ax.yaxis.get_offset_text().set_fontsize(30)
    ax.yaxis.get_offset_text().set_color("black")
    ax.yaxis.OFFSETTEXTPAD = 10
    ax.yaxis.set_major_formatter(formatter)
    ax.text(0.94, 0.91, rf"$p_{{D^0\perp}} \in ({pt_lo:g}, {pt_hi:g})$ GeV", transform=ax.transAxes,
            ha="right", va="top", fontsize=28, bbox=dict(facecolor="white", edgecolor="none", pad=3),
            zorder=10)
    if sublabel:
        ax.text(0.94, 0.78, sublabel, transform=ax.transAxes, ha="right", va="top", fontsize=28,
                linespacing=1.5, bbox=dict(facecolor="white", edgecolor="none", pad=3), zorder=10)


def draw_cms(ax, rows):
    # CMS points with the statistical error, and the systematic error as a box. Returns the highest value.
    top = 0.0
    for y_lo, y_hi, value, stat, syst in rows:
        ax.errorbar(0.5 * (y_lo + y_hi), value, yerr=stat, fmt="o", color="black", markersize=9,
                    capsize=0, zorder=5)
        ax.fill_between([y_lo, y_hi], value - syst, value + syst, facecolor="none",
                        edgecolor="black", linewidth=1.2, zorder=4)
        top = max(top, value + syst)
    return top


def finish(fig, axes, handles, filename, legend_title=None, labelspacing=0.55):
    # Legend in the first panel, axis titles, and save.
    # legend in a white box, a little below the top; the pT and FF labels are drawn on top of the box
    legend = axes.flat[0].legend(handles=handles, loc="upper left", bbox_to_anchor=(0.04, 0.91), fontsize=26,
                                 frameon=True, facecolor="white", edgecolor="none", framealpha=1.0,
                                 title=legend_title, title_fontsize=26, alignment="left",
                                 labelspacing=labelspacing, borderpad=0.15, borderaxespad=0.25, handlelength=1.3,
                                 handletextpad=0.4)
    legend.set_zorder(5)
    for ax in axes[-1]:
        ax.set_xlabel(r"$y$", labelpad=12, fontsize=34)
    fig.supylabel(r"$d\sigma/dy\,dp_{D^0\perp}$ [mb/GeV]", fontsize=34, x=0.024)
    fig.tight_layout(h_pad=1.5, w_pad=3.0, rect=(0.035, 0, 1, 1))
    os.makedirs("plots", exist_ok=True)
    fig.savefig(filename, bbox_inches="tight")
    print(f"Saved: {filename}")


def draw(panels, filename, cms=None, frag_runs=None, sublabel=FLUX_LABEL):
    # frag_runs: {frag: [central results, scale variations...]}. Default: the Pb+Pb results.
    if frag_runs is None:
        frag_runs = {}
        for frag in FRAGS:
            runs = load_frag(frag)
            if runs is not None:
                frag_runs[frag] = runs
    if not frag_runs:
        sys.exit("No results found -- run the scripts in run_scripts/ first.")

    nrows = int(np.ceil(len(panels) / 2))
    fig, axes = plt.subplots(nrows, 2, figsize=(16, 5.67 * nrows), sharex=True, squeeze=False)
    y_min = min(panel[2][0] for panel in panels)
    y_max = max(panel[2][-1] for panel in panels)
    has_cms = cms is not None and any(panel[3] for panel in panels)
    n_legend = len(frag_runs) + 1 + (1 if has_cms else 0)   # + the band entry

    for ax, (pt_lo, pt_hi, y_edges, show_cms) in zip(axes.flat, panels):
        top = 0.0
        for frag, runs in frag_runs.items():
            _, color, linestyle = FRAGS[frag]
            central, low, high = band(runs, pt_lo, pt_hi, y_edges)
            # one horizontal line per y bin, without the vertical lines joining the bins
            ax.hlines(central, y_edges[:-1], y_edges[1:], color=color, linestyle=linestyle, lw=3)
            for y_lo, y_hi, lo, hi in zip(y_edges[:-1], y_edges[1:], low, high):
                ax.fill_between([y_lo, y_hi], lo, hi, color=color, alpha=0.25, linewidth=0)
            top = max(top, np.nanmax(high))
        if show_cms and cms is not None and (pt_lo, pt_hi) in cms:
            top = max(top, draw_cms(ax, cms[(pt_lo, pt_hi)]))
        first = ax is axes.flat[0]
        style_panel(ax, pt_lo, pt_hi, y_min, y_max, top, legend_entries=n_legend if first else 0,
                    sublabel=sublabel if first else None)

    handles = [Line2D([0], [0], color=color, linestyle=linestyle, lw=3, label=label)
               for frag, (label, color, linestyle) in FRAGS.items() if frag in frag_runs]
    handles.append(BAND_HANDLE)
    if has_cms:
        handles.append(Line2D([0], [0], color="black", marker="o", linestyle="none", markersize=9, label="CMS"))
    finish(fig, axes, handles, filename)


def main():
    draw(FIG9_PANELS, f"plots/D0_incl_CMS_{CHANNEL}{SUFFIX}.pdf", cms=load_cms_data())
    draw(NARROW_PANELS, f"plots/D0_incl_bins_{CHANNEL}{SUFFIX}.pdf")


if __name__ == "__main__":
    main()
