"""
This script makes a plot that compares our theory predictions
(three fragmentation schemes: BCFY, Kniehl-Kramer, and LHAPDF)
to the CMS data for D0 mesons in Pb+Pb collisions.

CMS data comes from HEPData (doi:10.17182/hepdata.156822).
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

# make the plot look nicer 
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

# these are the pT bins and, for each pT bin, the y bins we need
PT_BINS_Y_BINS = [
    (2.0, 5.0, [(-1.0, 1.0)]),
    (5.0, 8.0, [(-2.0, -1.0), (-1.0, 0.0), (0.0, 1.0), (1.0, 2.0)]),
    (8.0, 12.0, [(-2.0, -1.0), (-1.0, 0.0), (0.0, 1.0), (1.0, 2.0)]),
]


def load_cms_data():
    # this function reads the 3 CMS csv files and puts everything
    # into one big dictionary so we can look it up later
    data = {}
    pattern = "data/HEPData-ins2968597-v1-D^0_cross_section_for_*_GeV_in_PbPb_UPCs.csv"
    file_list = glob.glob(pattern)

    for filename in file_list:
        # figure out which pT bin this file is for by looking at the name
        match = re.search(r"for_(\d+)_+p_T_+(\d+)_GeV", filename)
        pt_lo = float(match.group(1))
        pt_hi = float(match.group(2))

        f = open(filename)
        all_lines = f.readlines()
        f.close()

        # keep only the lines that are not comments and not empty
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
                # this is the header line, skip it
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
    # simple linear interpolation, same thing np.interp does,
    # written out by hand so it's clear what's going on
    return float(np.interp(x, x_list, y_list))


def pt_integrate(pt_values, cs_values, pt_lo, pt_hi):
    # put the points in order of increasing pT first
    pairs = list(zip(pt_values, cs_values))
    pairs.sort()

    sorted_pt = []
    sorted_cs = []
    for pt, cs in pairs:
        sorted_pt.append(pt)
        sorted_cs.append(cs)

    # get the value of the curve exactly at the edges of the bin
    edge_lo = interpolate(pt_lo, sorted_pt, sorted_cs)
    edge_hi = interpolate(pt_hi, sorted_pt, sorted_cs)

    # now build the list of x and y points to integrate:
    # bin edge, then all the grid points inside the bin, then the other edge
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

        # for every y value we have, integrate the cross section
        # over this pT bin
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
            # collect the y values that fall inside this y bin
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


def compute_lhapdf_scale_theory_points():
    # LHAPDF's error bars come from changing the factorization scale Q
    # (instead of using the 101 replicas). We need 3 runs:
    # the normal one (Q = mt0) and the two "out/lhapdf_scale" ones
    # (Q = mt0/2 and Q = mt0*2). If those two haven't been produced
    # yet, we just don't draw the error bars.
    central_pattern = "files/D0_incl_LHAPDF_An0n_Pb_y*.dat"
    low_pattern = "out/lhapdf_scale/factor_0.5/files/D0_incl_LHAPDF_An0n_Pb_y*.dat"
    high_pattern = "out/lhapdf_scale/factor_2.0/files/D0_incl_LHAPDF_An0n_Pb_y*.dat"

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

            # low_value and high_value are not necessarily in the
            # right order, so use min/max to be safe
            lowest = min(low_value, high_value, central)
            highest = max(low_value, high_value, central)

            y_points.append((y_lo, y_hi, central, lowest, highest))
        out.append((pt_lo, pt_hi, y_points))

    return out


def compute_lhapdf_replica_theory_points():
    # LHAPDF's other source of uncertainty: the 100 MC fit replicas (as
    # opposed to the scale variation above, which always uses the central
    # member). member_0000 is the central fit; members 0001-0100 are the
    # replicas the 68% band is built from. These come from
    # run_lhapdf_members_oberon.sh; if it hasn't been run yet, we skip this
    # band.
    member_dirs = sorted(glob.glob("files/lhapdf/member_*"))
    if len(member_dirs) < 2:
        return None

    central_theory = compute_theory_points(
        os.path.join(member_dirs[0], "files/D0_incl_LHAPDF_An0n_Pb_y*.dat")
    )
    replica_theory = []
    for member_dir in member_dirs[1:]:
        pattern = os.path.join(member_dir, "files/D0_incl_LHAPDF_An0n_Pb_y*.dat")
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


def compute_bk_posterior_theory_points(frag_type="LHAPDF"):
    # Uncertainty source shared by every fragmentation scheme: 100 posterior
    # samples of the BK initial-condition fit ("LOmvefit" -- Q_{s,0}^2, e_c,
    # C^2, sigma0/2; see bk/posteriorsamples_100_LOmvefit.dat), each
    # independently BK-evolved to its own dipole amplitude
    # (data/Pb/bk_posterior/member_<NNNN>), with the fragmentation function
    # held fixed at frag_type's central member/scale -- the opposite of
    # compute_lhapdf_replica_theory_points, which varies the fragmentation
    # function at fixed dipole. Comes from run_bk_posterior_members_oberon.sh
    # (FRAG_TYPE=<frag_type> when it was run).
    #
    # Unlike the LHAPDF replica set, there is no reason posterior sample
    # "member 0" should coincide with the real central (data/Pb/mve)
    # dipole -- the 100 samples are independent draws from the fit
    # posterior, not variations around a labeled central member. So instead
    # of splitting off one member as "the" central value (as
    # compute_lhapdf_replica_theory_points does), we take the *fractional*
    # 16th/84th percentile spread of all 100 samples around their own
    # median, to be re-applied to the real central value by whichever
    # caller combines it in -- same reasoning as the fractional workaround
    # used for the once-stale LHAPDF replica set.
    member_dirs = sorted(glob.glob("data/Pb/bk_posterior/member_*"))
    if len(member_dirs) < 2:
        return None

    member_theory = []
    for member_dir in member_dirs:
        pattern = os.path.join(member_dir, f"files/D0_incl_{frag_type}_An0n_Pb_y*.dat")
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
    # Three independent LHAPDF-curve uncertainty sources: scale variation,
    # fit (replica) uncertainty, and now the BK initial-condition
    # (posterior-sample dipole) uncertainty. Combined in quadrature since
    # each varies one thing (Q, the fragmentation-function fit, or the
    # dipole amplitude) while holding the other two fixed.
    scale = compute_lhapdf_scale_theory_points()
    replicas = compute_lhapdf_replica_theory_points()
    bk = compute_bk_posterior_theory_points("LHAPDF")
    if scale is None or replicas is None:
        return None, False

    out = []
    for i in range(len(PT_BINS_Y_BINS)):
        pt_lo, pt_hi, y_bins = PT_BINS_Y_BINS[i]
        y_points = []
        for j in range(len(y_bins)):
            y_lo, y_hi = y_bins[j]

            _, _, central, lo_scale, hi_scale = scale[i][2][j]
            _, _, _, lo_repl, hi_repl = replicas[i][2][j]

            lo_sq = (central - lo_scale) ** 2 + (central - lo_repl) ** 2
            hi_sq = (hi_scale - central) ** 2 + (hi_repl - central) ** 2

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
    # BCFY and Kniehl-Kramer have no scale-variation or fit-replica study
    # (they're fixed analytic forms, not a PDF/FF-fit with its own members)
    # -- so the only uncertainty source available for them is the BK
    # initial-condition (posterior-sample dipole) one, applied on top of
    # that scheme's own central prediction.
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


# list of (data file pattern, name to show in legend, marker shape, offset,
# FRAG_TYPE tag used to look up that scheme's BK-posterior sweep output)
FRAG_SCHEMES = [
    ("files/D0_incl_BCFY_An0n_Pb_y*.dat", "BCFY", "o", 0.06, "BCFY"),
    ("files/D0_incl_KniehlKramer_An0n_Pb_y*.dat", "Kniehl-Kramer", "x", -0.06, "KniehlKramer"),
    ("files/D0_incl_LHAPDF_An0n_Pb_y*.dat", "HymnD", "+", -0.12, "LHAPDF"),
]


def main():
    cms_data = load_cms_data()

    plt.figure(figsize=(6.5, 5.5))

    colors = {2.0: "tab:green", 5.0: "tab:red", 8.0: "tab:blue"}
    labels = {
        2.0: r"$2 < k_{D\perp} < 5$",
        5.0: r"$5 < k_{D\perp} < 8$",
        8.0: r"$8 < k_{D\perp} < 12$",
    }

    # draw a box for every CMS data point (value +/- stat and syst error)
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
                fill=False, edgecolor=color, linewidth=1.2,
            )
            plt.gca().add_patch(box)

    # draw a marker for each theory scheme
    for pattern, frag_label, marker, dx, frag_type in FRAG_SCHEMES:
        theory = compute_theory_points(pattern)
        for pt_lo, pt_hi, y_points in theory:
            color = colors[pt_lo]
            y_centers = []
            values = []
            for y_lo, y_hi, avg in y_points:
                y_centers.append(0.5 * (y_lo + y_hi))
                values.append(avg)
            plt.plot(y_centers, values, marker, color=color, markerfacecolor="none")

    # draw each scheme's uncertainty band, if we have one -- same faded-box
    # style for all of them (same color as the CMS data box in that pT
    # bin). HymnD (LHAPDF) gets scale variation + fit/replica + BK
    # initial-condition combined in quadrature; BCFY and Kniehl-Kramer have
    # no scale/replica study of their own, so they only get the
    # BK-initial-condition band, on top of their own central prediction.
    bands_drawn = {}
    for pattern, frag_label, marker, dx, frag_type in FRAG_SCHEMES:
        if frag_type == "LHAPDF":
            band, bk_included = compute_combined_theory_points()
        else:
            band = compute_bk_only_theory_points(pattern, frag_type)
            bk_included = band is not None
        if band is None:
            continue
        bands_drawn[frag_type] = bk_included
        for pt_lo, pt_hi, y_points in band:
            color = colors[pt_lo]
            for y_lo, y_hi, central, low, high in y_points:
                box = plt.Rectangle(
                    (y_lo, low), y_hi - y_lo, high - low,
                    facecolor=color, edgecolor="none", alpha=0.3,
                )
                plt.gca().add_patch(box)

    # legend for the theory schemes
    handles = []
    for pattern, frag_label, marker, dx, frag_type in FRAG_SCHEMES:
        one_handle = plt.Line2D(
            [0], [0], color="black", marker=marker, linestyle="",
            markerfacecolor="none", label=frag_label,
        )
        handles.append(one_handle)
    scheme_legend = plt.legend(
        handles=handles, loc="lower right", title="Fragmentation scheme", framealpha=1,
    )
    plt.gca().add_artist(scheme_legend)

    # legend for the CMS pT bins 
    style_handles = []
    for pt_lo in colors:
        patch = Patch(fill=False, edgecolor=colors[pt_lo], linewidth=1, label=labels[pt_lo])
        style_handles.append(patch)
    plt.legend(handles=style_handles, loc="lower left", title="CMS data", framealpha=1)

    plt.yscale("log")
    plt.xlabel(r"$y_D$")
    plt.ylabel(r"$d\sigma/dk_{D\perp}dy_D$ [mb/GeV]")
    plt.title(r"Pb + Pb $\to$ D$^0$ + X  (An0n, $\sqrt{s}=5.36$ TeV)")

    # Figure caption, split across lines rather than one long figtext -- a
    # single very long line forces bbox_inches="tight" to widen the whole
    # saved canvas to fit it, shrinking the plot itself.
    caption_lines = []
    if "LHAPDF" in bands_drawn:
        caption_lines.append(
            r"\textbf{HymnD band.} Combines in quadrature: factorization-scale "
            r"variation ($Q=0.5$-$2\times m_T$, $m_T^2=m_c^2+k_{D\perp}^2$),"
        )
        caption_lines.append(r"fit (replica) uncertainty from the 100 HymnD fit replicas")
        if bands_drawn["LHAPDF"]:
            caption_lines[-1] += r","
            caption_lines.append(
                r"and BK initial-condition uncertainty from a 100-sample dipole-amplitude posterior."
            )
        else:
            caption_lines[-1] += r"."

    other_bk_only = [ft for ft in ("BCFY", "KniehlKramer") if bands_drawn.get(ft)]
    if other_bk_only:
        caption_lines.append(
            r"\textbf{" + "/".join(other_bk_only) + r" bands.} BK initial-condition "
            r"uncertainty (100-sample dipole-amplitude posterior) only,"
        )
        caption_lines.append(r"these schemes have no scale/replica study.")

    for i, line in enumerate(caption_lines):
        plt.figtext(0.01, -0.05 - 0.03 * i, line, fontsize=7, color="dimgray", ha="left")

    plt.xlim(-2.4, 2.4)
    plt.ylim(1e-4, 1e1)
    plt.gca().xaxis.set_major_locator(MultipleLocator(1))
    plt.savefig("plots/D0_incl_dsigma_dkDperp_dyD_PbPb.pdf", dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
