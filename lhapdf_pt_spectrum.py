"""
Plots the D0 pT spectrum (dsigma/dy dpT vs pT) for the LHAPDF (HymnD)
fragmentation scheme, one curve per rapidity, over the y in [-3, 3] range
backfilled into files/D0_incl_LHAPDF_An0n_Pb_y*.dat.
"""
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

Y_MIN, Y_MAX = -3.0, 3.0


def main():
    results = load_results("files/D0_incl_LHAPDF_An0n_Pb_y*.dat")

    y_values = sorted(y for y in results if Y_MIN - 1e-9 <= y <= Y_MAX + 1e-9)

    norm = mcolors.Normalize(vmin=Y_MIN, vmax=Y_MAX)
    cmap = cm.coolwarm

    fig, ax = plt.subplots(figsize=(6.5, 5.5))

    for y in y_values:
        points = sorted(results[y])
        pt_values = [p[0] for p in points]
        cross_section_values = [p[1] for p in points]
        ax.plot(pt_values, cross_section_values, color=cmap(norm(y)), lw=1.5)

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label(r"$y_D$")

    ax.set_yscale("log")
    ax.set_xlabel(r"$p_{D^0}$ [GeV]")
    ax.set_ylabel(r"$d\sigma/dy\,dp_T$ [mb/GeV]")
    ax.set_title(r"HymnD: Pb + Pb $\to$ D$^0$ + X (An0n, $\sqrt{s}=5.36$ TeV)")
    ax.set_xlim(0, 12)

    plt.tight_layout()
    plt.savefig("plots/D0_incl_LHAPDF_pt_spectrum_by_y.pdf", dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
