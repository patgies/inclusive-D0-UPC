"""Quick plot of the HymnD pT spectrum."""

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

from cross_section import load_results

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 16,
    "axes.titlesize": 18,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 8,
    "ytick.major.size": 8,
    "xtick.minor.size": 4,
    "ytick.minor.size": 4,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
})

Y_VALUES = [0, 1, 2, 3]


def main():
    # load the data files
    results = load_results("files/D0_incl_HymnD_An0n_Pb_y*.dat")

    # choose the rapidity values you want to plot
    y_list = []
    for y in Y_VALUES:
        if float(y) in results:
            y_list.append(y)

    cmap = cm.coolwarm
    norm = mcolors.Normalize(vmin=min(Y_VALUES), vmax=max(Y_VALUES))

    fig, ax = plt.subplots(figsize=(6.5, 5.5))

    legend_handles = []
    for y in y_list:
        data = sorted(results[float(y)])
        pt = [p[0] for p in data]
        sigma = [p[1] for p in data]
        line, = ax.plot(pt, sigma, color=cmap(norm(y)), lw=1.2, alpha=0.75)
        legend_handles.append((line, r"$y_D = %d$" % y))

    lines = [h for h, _ in legend_handles]
    labels = [lab for _, lab in legend_handles]
    ax.legend(lines, labels, loc="upper right", fontsize=10, frameon=False)

    ax.set_yscale("log")
    ax.set_xlabel(r"$p_{D^0}$ [GeV]")
    ax.set_ylabel(r"$d\sigma/dy\,dp_T$ [mb/GeV]")
    ax.set_title(r"Pb + Pb $\to$ D$^0$ + X  (An0n, $\sqrt{s}=5.36$ TeV)")
    ax.set_xlim(0, 12)
    ax.set_ylim(1e-7, 1e2)

    plt.tight_layout()
    plt.savefig("plots/D0_incl_HymnD_pt_spectrum.pdf", dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
