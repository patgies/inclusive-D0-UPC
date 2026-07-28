import glob
import math
import os
import statistics
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.integrate import simpson


pi = math.pi
alpha_em = 1 / 137
e_charm_squared = 4 / 9  # (2/3)^2 electric charge of the charm quark
Nc = 3

# unit conversion (GeV^-2 -> mb)
FMGEV = 5.068
GEVSQR_TO_NB = 1.0e7 / (FMGEV * FMGEV)
GEVSQR_TO_MB = GEVSQR_TO_NB * 1e-6

# Overall normalization for the nucleus (Pb) target.
# b (glauber_mve_<b>) is in GeV^-1.
# The (2*pi)^4 comes from the prefactor and the momentum-space dipole integral.
factor_A = alpha_em * e_charm_squared * Nc / (2 * pi) ** 4

A = 208
ATA = 30.756  # dimensionless


def read_rapidity(filename):
    # look for the line that says "fixed rapidity y : <value>"
    # and give back that value as a number
    f = open(filename)
    for line in f:
        if "fixed rapidity y" in line:
            parts = line.split(":")
            text_after_colon = parts[-1]
            f.close()
            return float(text_after_colon)
    f.close()
    raise ValueError("no rapidity header found in " + filename)


def read_data_file(filename):
    # read a "b  pT  dsigma_dy" file and put the 3 columns into 3 lists
    b_list = []
    pt_list = []
    dsigma_list = []

    f = open(filename)
    for line in f:
        line = line.strip()
        if line == "":
            continue
        if line.startswith("#"):
            continue
        columns = line.split()
        b_list.append(float(columns[0]))
        pt_list.append(float(columns[1]))
        dsigma_list.append(float(columns[2]))
    f.close()

    return b_list, pt_list, dsigma_list


def group_by_pt(b_list, pt_list, dsigma_list):
    # put the (b, dsigma) pairs into a dictionary, one list per pT value
    groups = {}
    for i in range(len(pt_list)):
        b = b_list[i]
        pt = pt_list[i]
        dsigma = dsigma_list[i]
        if pt not in groups:
            groups[pt] = []
        groups[pt].append((b, dsigma))
    return groups


def integrate_over_b(pairs):
    # integrate b * dsigma over b using Simpson's rule
    pairs = list(pairs)
    pairs.sort()  # (b, dsigma) tuples sort by b first, which is what we want

    b_values = []
    weighted_values = []
    for pair in pairs:
        b = pair[0]
        dsigma = pair[1]
        b_values.append(b)
        weighted_values.append(b * dsigma)

    return simpson(weighted_values, x=b_values)


def load_results(pattern):
    # read every file that matches "pattern" and return
    # a dictionary: y value -> list of (pT, cross section)
    results = {}

    file_list = sorted(glob.glob(pattern))
    for filename in file_list:
        y = read_rapidity(filename)
        b_list, pt_list, dsigma_list = read_data_file(filename)
        pt_groups = group_by_pt(b_list, pt_list, dsigma_list)

        results[y] = []
        for pt in pt_groups:
            pairs = pt_groups[pt]
            b_integral = integrate_over_b(pairs)

            # (2*pi)^2 comes from d^2p -> dp*p*(2*pi) (the pT integral) and
            # the angular part of the 2D impact-parameter integral.
            cross_section = (2 * pi) * b_integral * factor_A * (2 * pi) * pt * GEVSQR_TO_MB

            results[y].append((pt, cross_section))

    return results


def load_lhapdf_band(member_dir_pattern):
    # Combine one run_many_Pb.sh output set per LHAPDF replica member into
    # a mean +/- standard-deviation band per rapidity. Returns an empty
    # dictionary if there are no member directories yet (that is, if
    # run_lhapdf_members.sh hasn't been run).
    member_dirs = sorted(glob.glob(member_dir_pattern))
    if len(member_dirs) == 0:
        return {}

    per_member_results = []
    for member_dir in member_dirs:
        one_pattern = os.path.join(member_dir, "files/D0_incl_LHAPDF_An0n_Pb_y*.dat")
        per_member_results.append(load_results(one_pattern))

    band = {}
    for y in per_member_results[0]:
        pt_values = []
        for pair in per_member_results[0][y]:
            pt_values.append(pair[0])
        pt_values.sort()

        cross_sections_by_pt = []
        for pt in pt_values:
            values = []
            for results in per_member_results:
                points = dict(results[y])
                values.append(points[pt])
            cross_sections_by_pt.append(values)

        means = []
        stds = []
        for values in cross_sections_by_pt:
            means.append(statistics.mean(values))
            stds.append(statistics.pstdev(values))

        band[y] = (pt_values, means, stds)

    return band


def main():
    results_g1 = load_results("files/D0_incl_KniehlKramer_An0n_G1_Pb_y*.dat")
    results_no_g1 = load_results("files/D0_incl_KniehlKramer_An0n_Pb_y*.dat")
    results_bcfy = load_results("files/D0_incl_BCFY_An0n_Pb_y*.dat")
    lhapdf_band = load_lhapdf_band("files/lhapdf/member_*")
    if lhapdf_band:
        results_lhapdf = None
    else:
        results_lhapdf = load_results("files/D0_incl_LHAPDF_An0n_Pb_y*.dat")

    # Plot the results for each rapidity y

    plt.figure(figsize=(7, 5))

    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    # pick one color per y value (only for y <= 2.0, that's all we plot)
    y_values_to_plot = []
    for y in results_g1:
        if y <= 2.0:
            y_values_to_plot.append(y)
    y_values_to_plot.sort()

    colors = {}
    for i in range(len(y_values_to_plot)):
        y = y_values_to_plot[i]
        colors[y] = color_cycle[i % len(color_cycle)]

    for y in sorted(results_g1):
        if y > 2.0:
            continue
        points = sorted(results_g1[y])
        pt_values = []
        cross_section_values = []
        for point in points:
            pt_values.append(point[0])
            cross_section_values.append(point[1])
        plt.plot(pt_values, cross_section_values, color=colors[y], linestyle="-", label="y=" + str(y))

    for y in sorted(results_no_g1):
        if y > 2.0:
            continue
        points = sorted(results_no_g1[y])
        pt_values = []
        cross_section_values = []
        for point in points:
            pt_values.append(point[0])
            cross_section_values.append(point[1])
        plt.plot(pt_values, cross_section_values, color=colors.get(y), linestyle="--")

    for y in sorted(results_bcfy):
        if y > 2.0:
            continue
        points = sorted(results_bcfy[y])
        pt_values = []
        cross_section_values = []
        for point in points:
            pt_values.append(point[0])
            cross_section_values.append(point[1])
        plt.plot(pt_values, cross_section_values, color=colors.get(y), linestyle=":")

    if lhapdf_band:
        for y in sorted(lhapdf_band):
            if y > 2.0:
                continue
            pt_values, means, stds = lhapdf_band[y]
            color = colors.get(y)
            lower = []
            upper = []
            for i in range(len(means)):
                mean = means[i]
                std = stds[i]
                lower_value = mean - std
                if lower_value < 1e-12:
                    lower_value = 1e-12
                lower.append(lower_value)
                upper.append(mean + std)
            plt.fill_between(pt_values, lower, upper, color=color, alpha=0.2, linewidth=0)
            plt.plot(pt_values, means, color=color, linestyle="-.")
    else:
        for y in sorted(results_lhapdf):
            if y > 2.0:
                continue
            points = sorted(results_lhapdf[y])
            pt_values = []
            cross_section_values = []
            for point in points:
                pt_values.append(point[0])
                cross_section_values.append(point[1])
            plt.plot(pt_values, cross_section_values, color=colors.get(y), linestyle="-.")

    plt.yscale("log")
    plt.xlabel(r"$p_{D^0}$ [GeV]")
    plt.ylabel(r"$d\sigma/dy\,dp_T$ [mb/GeV]")
    plt.title("Inclusive $D^0$ photoproduction, Pb+Pb UPC")

    y_legend = plt.legend(loc="upper right")
    plt.gca().add_artist(y_legend)

    if lhapdf_band:
        lhapdf_legend_text = "LHAPDF (mean +/- std over replicas)"
    else:
        lhapdf_legend_text = "LHAPDF (member 0, no errors)"

    style_handles = [
        Line2D([0], [0], color="black", linestyle="-", label="KniehlKramer, G1"),
        Line2D([0], [0], color="black", linestyle="--", label="KniehlKramer, no G1"),
        Line2D([0], [0], color="black", linestyle=":", label="BCFY"),
        Line2D([0], [0], color="black", linestyle="-.", label=lhapdf_legend_text),
    ]
    plt.legend(handles=style_handles, loc="lower left")

    plt.tight_layout()
    plt.savefig("plots/D0_incl_dsigma_dy_dpt_PbPb_g1.png", dpi=150)


if __name__ == "__main__":
    main()
