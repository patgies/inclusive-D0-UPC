"""
Compares the theory predictions (BCFY, Kniehl-Kramer, HymnD) with the CMS D0 data
in Pb+Pb UPCs (HEPData, doi:10.17182/hepdata.156822).
"""
import csv
import glob
import os
import re

import numpy as np
from scipy.integrate import simpson
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator

from cross_section import (
    load_results,
    factor_A,
    GEVSQR_TO_MB,
    pi,
)

# Paths are relative to the repo root.
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

# Plot style
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 11,
    "axes.labelsize": 18,
    "axes.titlesize": 20,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
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


PT_BINS_Y_BINS = [
    (2.0, 5.0, [(-1.0, 1.0)]),
    (5.0, 8.0, [(-2.0, -1.0), (-1.0, 0.0), (0.0, 1.0), (1.0, 2.0)]),
    (8.0, 12.0, [(-2.0, -1.0), (-1.0, 0.0), (0.0, 1.0), (1.0, 2.0)]),
]


def load_cms_data():
    # Reads the 3 CMS csv files into one dictionary.
    data = {}
    pattern = "input/CMS_data/HEPData-ins2968597-v1-D^0_cross_section_for_*_GeV_in_PbPb_UPCs.csv"
    file_list = glob.glob(pattern)

    for filename in file_list:
        # The pT bin is in the file name.
        match = re.search(r"for_(\d+)_+p_T_+(\d+)_GeV", filename)
        pt_lo = float(match.group(1))
        pt_hi = float(match.group(2))

        f = open(filename)
        all_lines = f.readlines()
        f.close()

        # Skip comments and empty lines.
        good_lines = []
        for line in all_lines:
            if line.startswith("#"):
                continue
            if line.strip() == "":
                continue
            good_lines.append(line)

        reader = csv.reader(good_lines)
        rows = []
        first_line = True
        for row in reader:
            if first_line:
                # header line
                first_line = False
                continue
            y = float(row[0])
            y_lo = float(row[1])
            y_hi = float(row[2])
            value = float(row[3])
            stat_p = float(row[4])
            stat_m = float(row[5])
            syst_p = float(row[6])
            syst_m = float(row[7])
            stat = 0.5 * (stat_p - stat_m)
            syst = 0.5 * (syst_p - syst_m)
            rows.append((y_lo, y_hi, value, stat, syst))

        rows.sort()
        data[(pt_lo, pt_hi)] = rows

    return data


def interpolate(x, x_list, y_list):
    return float(np.interp(x, x_list, y_list))


def pt_integrate(pt_values, cs_values, pt_lo, pt_hi):
    # Sort by pT.
    pairs = list(zip(pt_values, cs_values))
    pairs.sort()

    sorted_pt = []
    sorted_cs = []
    for pt, cs in pairs:
        sorted_pt.append(pt)
        sorted_cs.append(cs)

    # Value of the curve at the edges of the bin
    edge_lo = interpolate(pt_lo, sorted_pt, sorted_cs)
    edge_hi = interpolate(pt_hi, sorted_pt, sorted_cs)

    # Points to integrate: bin edge, grid points inside the bin, other edge.
    x_points = [pt_lo]
    y_points = [edge_lo]
    for i in range(len(sorted_pt)):
        pt = sorted_pt[i]
        cs = sorted_cs[i]
        if pt > pt_lo and pt < pt_hi:
            x_points.append(pt)
            y_points.append(cs)
    x_points.append(pt_hi)
    y_points.append(edge_hi)

    return simpson(y_points, x=x_points)


def compute_theory_points(pattern):
    # results is a dictionary: y value -> list of (pT, cross section)
    results = load_results(pattern)

    out = []
    for pt_lo, pt_hi, y_bins in PT_BINS_Y_BINS:

        # Integrate over this pT bin for every y.
        per_y = {}
        for y in results:
            points = results[y]
            pt_list = []
            cs_list = []
            for point in points:
                pt_list.append(point[0])
                cs_list.append(point[1])
            per_y[y] = pt_integrate(pt_list, cs_list, pt_lo, pt_hi)

        y_points = []
        for y_lo, y_hi in y_bins:
            # y values inside this y bin
            y_list = []
            for y in per_y:
                if y_lo - 1e-9 <= y <= y_hi + 1e-9:
                    y_list.append(y)
            y_list.sort()

            cs_list = []
            for y in y_list:
                cs_list.append(per_y[y])

            y_integral = simpson(cs_list, x=y_list)
            avg = y_integral / (y_hi - y_lo) / (pt_hi - pt_lo)
            y_points.append((y_lo, y_hi, avg))

        out.append((pt_lo, pt_hi, y_points))

    return out


def compute_hymnd_scale_theory_points():
    # HymnD scale band: Q = mt0 (central), mt0/2 and 2*mt0.
    # No error bars if the scale runs are missing.
    # Not 0.25 and 4 as for Kniehl-Kramer: 0.25*mt0 is below the D0 mass for pD0 < 7.3 GeV.
    central_pattern = "output/An0n/HymnD_band/member_0000/D0_incl_HymnD_An0n_Pb_y*.dat"
    low_pattern = "output/An0n/scale_variation/HymnD/factor_0.5/D0_incl_HymnD_An0n_Pb_y*.dat"
    high_pattern = "output/An0n/scale_variation/HymnD/factor_2.0/D0_incl_HymnD_An0n_Pb_y*.dat"

    if len(glob.glob(central_pattern)) == 0:
        return None
    if len(glob.glob(low_pattern)) == 0:
        return None
    if len(glob.glob(high_pattern)) == 0:
        return None

    central_theory = compute_theory_points(central_pattern)
    low_theory = compute_theory_points(low_pattern)
    high_theory = compute_theory_points(high_pattern)

    out = []
    for i in range(len(PT_BINS_Y_BINS)):
        pt_lo, pt_hi, y_bins = PT_BINS_Y_BINS[i]
        y_points = []
        for j in range(len(y_bins)):
            y_lo, y_hi = y_bins[j]

            central = central_theory[i][2][j][2]
            low_value = low_theory[i][2][j][2]
            high_value = high_theory[i][2][j][2]

            # low_value can be larger than high_value
            lowest = min(low_value, high_value, central)
            highest = max(low_value, high_value, central)

            y_points.append((y_lo, y_hi, central, lowest, highest))
        out.append((pt_lo, pt_hi, y_points))

    return out


def compute_hymnd_replica_theory_points():
    # HymnD replica band: member_0000 is the central fit, members 0001-0100 give the 68% band.
    # Skipped if the members are missing (run_scripts/run_members.sh).
    member_dirs = sorted(glob.glob("output/An0n/HymnD_band/member_*"))
    if len(member_dirs) < 2:
        return None

    central_theory = compute_theory_points(
        os.path.join(member_dirs[0], "D0_incl_HymnD_An0n_Pb_y*.dat")
    )
    replica_theory = []
    for member_dir in member_dirs[1:]:
        pattern = os.path.join(member_dir, "D0_incl_HymnD_An0n_Pb_y*.dat")
        replica_theory.append(compute_theory_points(pattern))

    out = []
    for i in range(len(PT_BINS_Y_BINS)):
        pt_lo, pt_hi, y_bins = PT_BINS_Y_BINS[i]
        y_points = []
        for j in range(len(y_bins)):
            y_lo, y_hi = y_bins[j]

            central = central_theory[i][2][j][2]
            values = []
            for rep in replica_theory:
                values.append(rep[i][2][j][2])
            lo, hi = np.percentile(values, [16, 84])

            y_points.append((y_lo, y_hi, central, lo, hi))
        out.append((pt_lo, pt_hi, y_points))

    return out


def compute_kniehlkramer_scale_theory_points():
    # Kniehl-Kramer scale band: Q = mt0 (central), mt0/4 and 4*mt0 (run_scripts/run_scale_variation.sh).
    # No error bars if the scale runs are missing.
    central_pattern = "output/An0n/scale_variation/KniehlKramer/factor_1.0/D0_incl_KniehlKramer_An0n_Pb_y*.dat"
    low_pattern = "output/An0n/scale_variation/KniehlKramer/factor_0.25/D0_incl_KniehlKramer_An0n_Pb_y*.dat"
    high_pattern = "output/An0n/scale_variation/KniehlKramer/factor_4.0/D0_incl_KniehlKramer_An0n_Pb_y*.dat"

    if len(glob.glob(central_pattern)) == 0:
        return None
    if len(glob.glob(low_pattern)) == 0:
        return None
    if len(glob.glob(high_pattern)) == 0:
        return None

    central_theory = compute_theory_points(central_pattern)
    low_theory = compute_theory_points(low_pattern)
    high_theory = compute_theory_points(high_pattern)

    out = []
    for i in range(len(PT_BINS_Y_BINS)):
        pt_lo, pt_hi, y_bins = PT_BINS_Y_BINS[i]
        y_points = []
        for j in range(len(y_bins)):
            y_lo, y_hi = y_bins[j]

            central = central_theory[i][2][j][2]
            low_value = low_theory[i][2][j][2]
            high_value = high_theory[i][2][j][2]

            # low_value can be larger than high_value
            lowest = min(low_value, high_value, central)
            highest = max(low_value, high_value, central)

            y_points.append((y_lo, y_hi, central, lowest, highest))
        out.append((pt_lo, pt_hi, y_points))

    return out


def compute_bk_posterior_theory_points(frag_type="HymnD"):
    # BK band: 100 posterior samples of the BK initial condition, fragmentation function fixed.
    # No sample is the central one, so the band is the 16th-84th percentile
    # spread around the median, as a fraction.
    member_dirs = sorted(glob.glob("output/An0n/bk_band/member_*"))
    if len(member_dirs) < 2:
        return None

    member_theory = []
    for member_dir in member_dirs:
        pattern = os.path.join(member_dir, f"D0_incl_{frag_type}_An0n_Pb_y*.dat")
        if len(glob.glob(pattern)) == 0:
            return None
        member_theory.append(compute_theory_points(pattern))

    out = []
    for i in range(len(PT_BINS_Y_BINS)):
        pt_lo, pt_hi, y_bins = PT_BINS_Y_BINS[i]
        y_points = []
        for j in range(len(y_bins)):
            y_lo, y_hi = y_bins[j]

            values = []
            for member in member_theory:
                values.append(member[i][2][j][2])
            median = np.median(values)
            lo, hi = np.percentile(values, [16, 84])

            y_points.append((y_lo, y_hi, median, lo / median, hi / median))
        out.append((pt_lo, pt_hi, y_points))

    return out


def compute_combined_theory_points():
    # HymnD band: scale variation, replicas and BK, added in quadrature.
    scale = compute_hymnd_scale_theory_points()
    # replicas = compute_hymnd_replica_theory_points()
    replicas = None
    # bk = compute_bk_posterior_theory_points("HymnD")
    bk = None
    if scale is None:
        return None, False

    out = []
    for i in range(len(PT_BINS_Y_BINS)):
        pt_lo, pt_hi, y_bins = PT_BINS_Y_BINS[i]
        y_points = []
        for j in range(len(y_bins)):
            y_lo, y_hi = y_bins[j]

            _, _, central, lo_scale, hi_scale = scale[i][2][j]

            lo_sq = (central - lo_scale) ** 2
            hi_sq = (hi_scale - central) ** 2

            if replicas is not None:
                _, _, _, lo_repl, hi_repl = replicas[i][2][j]
                lo_sq += (central - lo_repl) ** 2
                hi_sq += (hi_repl - central) ** 2

            if bk is not None:
                _, _, _, frac_lo_bk, frac_hi_bk = bk[i][2][j]
                lo_bk = central * frac_lo_bk
                hi_bk = central * frac_hi_bk
                lo_sq += (central - lo_bk) ** 2
                hi_sq += (hi_bk - central) ** 2

            y_points.append((y_lo, y_hi, central, central - lo_sq ** 0.5, central + hi_sq ** 0.5))
        out.append((pt_lo, pt_hi, y_points))

    return out, bk is not None


def compute_bk_only_theory_points(pattern, frag_type):
    # BCFY and Kniehl-Kramer: only the BK band, around their central prediction.
    bk = compute_bk_posterior_theory_points(frag_type)
    if bk is None:
        return None

    central_theory = compute_theory_points(pattern)

    out = []
    for i in range(len(PT_BINS_Y_BINS)):
        pt_lo, pt_hi, y_bins = PT_BINS_Y_BINS[i]
        y_points = []
        for j in range(len(y_bins)):
            y_lo, y_hi, central = central_theory[i][2][j]
            _, _, _, frac_lo, frac_hi = bk[i][2][j]
            y_points.append((y_lo, y_hi, central, central * frac_lo, central * frac_hi))
        out.append((pt_lo, pt_hi, y_points))

    return out



FRAG_SCHEMES = [
    ("output/An0n/central_values/BCFY/D0_incl_BCFY_An0n_Pb_y*.dat", "BCFY", "+", 0.06, "BCFY", "empty"),
    ("output/An0n/scale_variation/KniehlKramer/factor_1.0/D0_incl_KniehlKramer_An0n_Pb_y*.dat", "Kniehl-Kramer", "o", -0.06, "KniehlKramer", "empty"),
    (
        "output/An0n/HymnD_band/member_0000/D0_incl_HymnD_An0n_Pb_y*.dat",
        r"HymnD",
        ".", -0.12, "HymnD", "solid",
    ),
]


def main():
    cms_data = load_cms_data()

    plt.figure(figsize=(6.5, 5.5))

    colors = {2.0: "#0a8a1e", 5.0: "#e00000", 8.0: "#0b3c8c"}
    labels = {
        2.0: r"$2 < p_{D\perp} < 5$",
        5.0: r"$5 < p_{D\perp} < 8$",
        8.0: r"$8 < p_{D\perp} < 12$",
    }


    for pt_lo, pt_hi, y_bins in PT_BINS_Y_BINS:
        color = colors[pt_lo]
        rows = cms_data[(pt_lo, pt_hi)]
        for row in rows:
            y_lo = row[0]
            y_hi = row[1]
            value = row[2]
            stat = row[3]
            syst = row[4]
            total_error = (stat ** 2 + syst ** 2) ** 0.5
            box = plt.Rectangle(
                (y_lo, value - total_error), y_hi - y_lo, 2 * total_error,
                fill=True, facecolor=color, alpha=0.15, edgecolor=color, linewidth=1.2,
            )
            plt.gca().add_patch(box)


    for pattern, frag_label, marker, dx, frag_type, style in FRAG_SCHEMES:
        theory = compute_theory_points(pattern)
        for pt_lo, pt_hi, y_points in theory:
            color = colors[pt_lo]
            y_centers = []
            values = []
            point_dx = 0.08 if frag_type == "HymnD" else (-0.08 if frag_type == "BCFY" else 0.0)
            for y_lo, y_hi, avg in y_points:
                y_centers.append(0.5 * (y_lo + y_hi) + point_dx)
                values.append(avg)
            if style == "hatched":

                pc = plt.scatter(
                    y_centers, values, marker=marker, s=36,
                    facecolor=color, edgecolor="black", linewidth=0.8, zorder=3,
                )
                pc.set_hatch("///")
            else:
                facecolor = color if style == "solid" else "none"
                markersize = 5 if marker == "o" else (8 if marker == "." else None)
                if style == "empty":
                    # White circle first, so the band does not show through the open marker.
                    plt.plot(
                        y_centers, values, marker, color="white", markerfacecolor="white",
                        markersize=markersize, zorder=2,
                    )
                plt.plot(
                    y_centers, values, marker, color=color, markerfacecolor=facecolor,
                    markersize=markersize, zorder=3,
                )


    bands_drawn = {}
    for pattern, frag_label, marker, dx, frag_type, style in FRAG_SCHEMES:
        if frag_type == "HymnD":
            band, bk_included = compute_combined_theory_points()
        elif frag_type == "KniehlKramer":
            band = compute_kniehlkramer_scale_theory_points()
            bk_included = False
        else:
            # band = compute_bk_only_theory_points(pattern, frag_type)
            # bk_included = band is not None
            band = None
        if band is None:
            continue
        bands_drawn[frag_type] = bk_included
        # Uncertainty as a shaded box.
        box_width = 0.02
        box_dx = 0.08 if frag_type == "HymnD" else 0.0
        for pt_lo, pt_hi, y_points in band:
            color = colors[pt_lo]
            for y_lo, y_hi, central, low, high in y_points:
                y_center = 0.5 * (y_lo + y_hi) + box_dx
                box = plt.Rectangle(
                    (y_center - box_width / 2, low), box_width, high - low,
                    facecolor=color, edgecolor="none", alpha=0.85,
                    fill=True, zorder=1,
                )
                plt.gca().add_patch(box)
                cap_half_width = box_width * 0.75
                for y_edge in (low, high):
                    plt.plot(
                        [y_center - cap_half_width, y_center + cap_half_width],
                        [y_edge, y_edge], color=color, linewidth=1.0, zorder=1,
                    )

    # Legend: fragmentation functions
    handles = []
    for pattern, frag_label, marker, dx, frag_type, style in FRAG_SCHEMES:
        if style == "hatched":

            one_handle = plt.scatter(
                [], [], marker=marker, s=36,
                facecolor="black", edgecolor="black", linewidth=0.8, label=frag_label,
            )
            one_handle.set_hatch("///")
        else:
            one_handle = plt.Line2D(
                [0], [0], color="black", marker=marker, linestyle="",
                markerfacecolor="black" if style == "solid" else "none", label=frag_label,
                markersize=5 if marker == "o" else (8 if marker == "." else None),
            )
        handles.append(one_handle)

    scheme_legend = plt.legend(
        handles=handles, title="Fragmentation function", frameon=False,
        loc="lower left", bbox_to_anchor=(0.31, 0.0), handletextpad=0.3,
        fontsize=12, title_fontsize=12,
    )
    plt.gca().add_artist(scheme_legend)

    # Legend: CMS pT bins
    style_handles = []
    for pt_lo in colors:
        patch = Patch(
            fill=True, facecolor=colors[pt_lo], alpha=0.15,
            edgecolor=colors[pt_lo], linewidth=1, label=labels[pt_lo],
        )
        style_handles.append(patch)
    plt.legend(
        handles=style_handles, title="CMS data", frameon=False,
        loc="lower left", bbox_to_anchor=(0.0, 0.0),
        fontsize=12, title_fontsize=12,
    )

    #separator = plt.Line2D(
    #    [0.32, 0.32], [0.0, 0.13], transform=plt.gca().transAxes,
    #    color="gray", linewidth=0.8,
    #)
    #plt.gca().add_line(separator)

    plt.yscale("log")
    plt.xlabel(r"$y$")
    plt.ylabel(r"$\mathrm{d}\sigma/\mathrm{d}y\,\mathrm{d}p_{D\perp}$ [mb/GeV]", labelpad=10)
    plt.title(r"Pb + Pb $\to$ D$^0$ + X ")

    # Figure caption
    caption_lines = []
    if "HymnD" in bands_drawn:
        caption_lines.append(
            r"\textbf{HymnD band.} Combines in quadrature: factorization-scale "
            r"variation ($Q=0.5$-$2\times \sqrt{m_c^2+k_{D\perp}^2}$),"
        )
        caption_lines.append(r"fit uncertainty from the 100 HymnD fit replicas")
        if bands_drawn["HymnD"]:
            caption_lines[-1] += r","
            caption_lines.append(
                r"and BK initial-condition uncertainty from a 100-sample dipole-amplitude posterior."
            )
        else:
            caption_lines[-1] += r"."

    if "KniehlKramer" in bands_drawn:
        caption_lines.append(
            r"\textbf{Kniehl-Kramer band.} Factorization-scale variation "
            r"($Q=0.5$-$2\times \sqrt{m_c^2+k_{D\perp}^2}$) only;"
        )
        caption_lines.append(r"no fit-replica study for this fixed analytic form.")

    other_bk_only = [ft for ft in ("BCFY",) if bands_drawn.get(ft)]
    if other_bk_only:
        caption_lines.append(
            r"\textbf{" + "/".join(other_bk_only) + r" bands.} BK initial-condition "
            r"uncertainty (100-sample dipole-amplitude posterior) only,"
        )
        caption_lines.append(r"these schemes have no scale/replica study.")

    # for i, line in enumerate(caption_lines):
    #     plt.figtext(0.01, -0.05 - 0.03 * i, line, fontsize=7, color="dimgray", ha="left")

    plt.xlim(-2.4, 2.4)
    plt.ylim(1e-4, 1e1)
    plt.gca().xaxis.set_major_locator(MultipleLocator(1))
    plt.savefig("plots/D0_incl_CMS_comparison.pdf", dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
