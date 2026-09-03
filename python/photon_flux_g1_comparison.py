import os

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator

from cross_section import load_results


os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 15,
    "axes.labelsize": 22,
    "axes.titlesize": 20,
    "xtick.labelsize": 22,
    "ytick.labelsize": 22,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 8,
    "ytick.major.size": 10,
    "xtick.minor.size": 4,
    "ytick.minor.size": 4,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
})

Y_VALUES_TO_PLOT = [ 0, 1, 2, 3]

PT_VALUES_TO_REPORT = [2.1, 5.1, 9.9]

LABEL_OFFSETS = {
    (0, 2.1): (0, 10),
    (1.0, 2.1): (0, 0),
    (2, 2.1): (0, -10),
    (3, 2.1): (0, -35),
    (0, 5.1): (0, 15),
    (1.0, 5.1): (0, -5),
    (2, 5.1): (0, -25),
    (3, 5.1): (0, -35),
    (0, 9.9): (0, 12),
    (1.0, 9.9): (0, -20),
    (2, 9.9): (0, -25),
    (3, 9.9): (0, -35),
}
DEFAULT_LABEL_OFFSETS = [(0, 20), (0, -34), (0, 20), (0, -34)]


def compute_percent_difference(results_g1, results_no_g1):
    # % diff of Gamma_AA=1 relative to the physical Gamma_AA(b) case:
    # (sigma_G1 - sigma_no_g1) / sigma_no_g1 * 100
    percent_diff_by_y = {}
    for y in Y_VALUES_TO_PLOT:
        if y not in results_g1 or y not in results_no_g1:
            continue
        g1_points = dict(results_g1[y])
        no_g1_points = dict(results_no_g1[y])

        row = []
        for pt in PT_VALUES_TO_REPORT:
            if pt not in g1_points or pt not in no_g1_points:
                continue
            percent_diff = (g1_points[pt] - no_g1_points[pt]) / no_g1_points[pt] * 100
            row.append((pt, percent_diff))
        percent_diff_by_y[y] = row
    return percent_diff_by_y


def format_percent(percent_diff):
    return "%d" % round(percent_diff)


def print_percent_difference(percent_diff_by_y):
    for y, row in percent_diff_by_y.items():
        print("y=" + str(y))
        for pt, percent_diff in row:
            print("  pT=%5.1f GeV:  %6.2f%%" % (pt, percent_diff))


def main():
    results_g1 = load_results("files/central/D0_incl_KniehlKramer_An0n_G1_Pb_y*.dat")
    results_no_g1 = load_results("files/central/D0_incl_KniehlKramer_An0n_Pb_y*.dat")

    percent_diff_by_y = compute_percent_difference(results_g1, results_no_g1)
    print_percent_difference(percent_diff_by_y)

    plt.figure(figsize=(8, 6.5))


    # Blues for y=0/y=1, browns for y=2/y=3.
    shades = ["#3f76c0", "#0d3f8f", "#b35a2e", "#4a1c08"]

    for i, y in enumerate(Y_VALUES_TO_PLOT):
        color = shades[i % len(shades)]

        if y in results_g1:
            points = sorted(results_g1[y])
            pt_values = [point[0] for point in points]
            cross_section_values = [point[1] for point in points]
            plt.plot(pt_values, cross_section_values, color=color, linestyle="-", label="y=" + str(y))

            g1_points = dict(results_g1[y])
            default_offset = DEFAULT_LABEL_OFFSETS[i % len(DEFAULT_LABEL_OFFSETS)]
            for pt, percent_diff in percent_diff_by_y.get(y, []):
                offset = LABEL_OFFSETS.get((y, pt), default_offset)
                plt.annotate(format_percent(percent_diff) + "\\%", xy=(pt, g1_points[pt]),
                             xytext=offset, textcoords="offset points",
                             color="black", fontsize=16, ha="center",
                             bbox=dict(facecolor="white", edgecolor="none",
                                       boxstyle="round,pad=0.2", alpha=1.0))

        if y in results_no_g1:
            points = sorted(results_no_g1[y])
            pt_values = [point[0] for point in points]
            cross_section_values = [point[1] for point in points]
            plt.plot(pt_values, cross_section_values, color=color, linestyle="--")

    plt.yscale("log")
    ymin, ymax = plt.ylim()
    plt.ylim(ymin, ymax * 3)
    plt.gca().yaxis.set_minor_locator(
        LogLocator(base=10, subs=(2, 3, 4, 5, 6, 7, 8, 9), numticks=100)
    )
    plt.xlabel(r"$p_{D\perp}$ [GeV]", labelpad=10)
    plt.ylabel(r"$\mathrm{d}\sigma/\mathrm{d}y\,\mathrm{d}p_{D\perp}$ [mb/GeV]", labelpad=10)
    plt.title(r"Pb + Pb $\to$ D$^0$ + X", pad=10, fontsize=24)

    y_legend = plt.legend(loc="lower left", fontsize=18, frameon=False)
    plt.gca().add_artist(y_legend)

    style_handles = [
        Line2D([0], [0], color="black", linestyle="-", label=r"$\Gamma_{AA}=1$"),
        Line2D([0], [0], color="black", linestyle="--", label=r"$\Gamma_{AA}(\mathbf{b})$"),
    ]
    plt.legend(handles=style_handles, loc="upper right", fontsize=20, frameon=False)

    plt.tight_layout()
    plt.savefig("plots/D0_incl_G1_flux_comparison_PbPb.pdf", dpi=600, bbox_inches="tight")


if __name__ == "__main__":
    main()
