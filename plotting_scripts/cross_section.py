import glob
import math
import os
import statistics
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.integrate import simpson


pi = math.pi
alpha_em = 1 / 137
e_charm_squared = 4 / 9  # charm charge squared
Nc = 3

# GeV^-2 -> mb
FMGEV = 5.068
GEVSQR_TO_NB = 1.0e7 / (FMGEV * FMGEV)
GEVSQR_TO_MB = GEVSQR_TO_NB * 1e-6

# Prefactor of the cross section. The dipole impact parameter b_d (called b below) is in GeV^-1.
factor_A = alpha_em * e_charm_squared * Nc / (2 * pi) ** 4

A = 208
ATA = 30.756  # dimensionless


def read_rapidity(filename):
    # Reads y from the header line "fixed rapidity y : <value>".
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
    # Reads the columns b, pT, dsigma of one file.
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
    # (b, dsigma) pairs for each pT
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
    # Integral of b * dsigma over b (Simpson's rule)
    pairs = list(pairs)
    pairs.sort()  # sort by b

    b_values = []
    weighted_values = []
    for pair in pairs:
        b = pair[0]
        dsigma = pair[1]
        b_values.append(b)
        weighted_values.append(b * dsigma)

    return simpson(weighted_values, x=b_values)


def load_results(pattern):
    # Reads the files that match pattern. Returns {y: [(pT, cross section), ...]}.
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

            # 2*pi from the b angle, 2*pi*pT from d2p -> dpT
            cross_section = (2 * pi) * b_integral * factor_A * (2 * pi) * pt * GEVSQR_TO_MB

            results[y].append((pt, cross_section))

    return results


def load_hymnd_band(member_dir_pattern):
    # Mean and standard deviation over the HymnD members, for each y.
    member_dirs = sorted(glob.glob(member_dir_pattern))
    if len(member_dirs) == 0:
        return {}

    per_member_results = []
    for member_dir in member_dirs:
        one_pattern = os.path.join(member_dir, "D0_incl_HymnD_An0n_Pb_y*.dat")
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
    results_g1 = load_results("output/An0n/central_values/KniehlKramer/D0_incl_KniehlKramer_An0n_G1_Pb_y*.dat")
    results_no_g1 = load_results("output/An0n/central_values/KniehlKramer/D0_incl_KniehlKramer_An0n_Pb_y*.dat")
    results_bcfy = load_results("output/An0n/central_values/BCFY/D0_incl_BCFY_An0n_Pb_y*.dat")
    results_hymnd = load_results("output/An0n/scale_variation/HymnD/factor_1.0/D0_incl_HymnD_An0n_Pb_y*.dat")
    hymnd_band = load_hymnd_band("output/An0n/HymnD_band/member_*")
    if hymnd_band:
        results_hymnd = None
    else:
        results_hymnd = load_results("output/An0n/central_values/HymnD/D0_incl_HymnD_An0n_Pb_y*.dat")

    # One curve per y

    plt.figure(figsize=(7, 5))

    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    # one color per y (only y <= 2 is plotted)
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

    for y in sorted(results_hymnd):
        if y > 2.0:
            continue
        points = sorted(results_hymnd[y])
        pt_values = []
        cross_section_values = []
        for point in points:
            pt_values.append(point[0])
            cross_section_values.append(point[1])
        plt.plot(pt_values, cross_section_values, color=colors.get(y), linestyle=(0, (3, 1, 1, 1)))

    if hymnd_band:
        for y in sorted(hymnd_band):
            if y > 2.0:
                continue
            pt_values, means, stds = hymnd_band[y]
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
        for y in sorted(results_hymnd):
            if y > 2.0:
                continue
            points = sorted(results_hymnd[y])
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

    if hymnd_band:
        hymnd_legend_text = "HymnD (mean +/- std over replicas)"
    else:
        hymnd_legend_text = "HymnD (member 0, no errors)"

    style_handles = [
        Line2D([0], [0], color="black", linestyle="-", label="KniehlKramer, G1"),
        Line2D([0], [0], color="black", linestyle="--", label="KniehlKramer, no G1"),
        Line2D([0], [0], color="black", linestyle=":", label="BCFY"),
        Line2D([0], [0], color="black", linestyle="-.", label=hymnd_legend_text),
        Line2D([0], [0], color="black", linestyle=(0, (3, 1, 1, 1)), label="HymnD"),
    ]
    plt.legend(handles=style_handles, loc="lower left")

    plt.tight_layout()
    plt.savefig("plots/D0_incl_dsigma_dy_dpt_PbPb_g1.png", dpi=150)


if __name__ == "__main__":
    main()
