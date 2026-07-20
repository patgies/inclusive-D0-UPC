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

#  unit conversion (GeV^-2 -> mb)
FMGEV = 5.068
GEVSQR_TO_NB = 1.0e7 / (FMGEV * FMGEV)
GEVSQR_TO_MB = GEVSQR_TO_NB * 1e-6

# Overall normalization for the nucleus (Pb) target.
# b (glauber_mve_<b>) is in GeV^-1.
# The (2*pi)^4 comes from the prefactor and the momentum-space dipole integral.
factor_A = alpha_em * e_charm_squared * Nc / (2 * pi) ** 4

A = 208
ATA = 30.756  # dimensionless

#Given (b, dsigma) pairs at a fixed pT, integrate b*dsigma over b using Simpson's rule.

def read_rapidity(filename):
    """Find the line '# fixed rapidity y : <value>' in the file header and return it."""
    with open(filename) as f:
        for line in f:
            if "fixed rapidity y" in line:
                text_after_colon = line.split(":")[-1]
                return float(text_after_colon)
    raise ValueError(f"no rapidity header found in {filename}")


def read_data_file(filename):
    """
    Read one 'b  pT  dsigma_dy' data file, skipping comment lines (starting
    with '#'). Returns three plain lists of the same length: b, pt, dsigma.
    """
    b_list = []
    pt_list = []
    dsigma_list = []

    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line == "" or line.startswith("#"):
                continue
            columns = line.split()
            b_list.append(float(columns[0]))
            pt_list.append(float(columns[1]))
            dsigma_list.append(float(columns[2]))

    return b_list, pt_list, dsigma_list


def group_by_pt(b_list, pt_list, dsigma_list):
    """Group the (b, dsigma) pairs by their pT value: {pt: [(b, dsigma), ...]}."""
    groups = {}
    for b, pt, dsigma in zip(b_list, pt_list, dsigma_list):
        if pt not in groups:
            groups[pt] = []
        groups[pt].append((b, dsigma))
    return groups


def integrate_over_b(pairs):
    """
    Given (b, dsigma) pairs at a fixed pT, integrate b*dsigma over b using
    Simpson's rule. This is the radial measure of the 2D impact-parameter
    (b) integral, int b * dsigma_dy(b) db.
    """
    pairs = sorted(pairs, key=lambda pair: pair[0])
    b_values = [pair[0] for pair in pairs]
    weighted_values = [b * dsigma for b, dsigma in pairs]
    return simpson(weighted_values, x=b_values)


def load_results(pattern):
    """Read all files matching pattern and return {y: [(pt, cross_section), ...]}."""
    results = {}

    for filename in sorted(glob.glob(pattern)):
        y = read_rapidity(filename)
        b_list, pt_list, dsigma_list = read_data_file(filename)
        pt_groups = group_by_pt(b_list, pt_list, dsigma_list)

        results[y] = []
        for pt, pairs in pt_groups.items():
            b_integral = integrate_over_b(pairs)

            # (2*pi)^2 comes from d^2p -> dp*p*(2*pi) (the pT integral) and
            # the angular part of the 2D impact-parameter integral.
            cross_section = (2 * pi) * b_integral * factor_A * (2 * pi) * pt * GEVSQR_TO_MB

            results[y].append((pt, cross_section))

    return results


def load_lhapdf_band(member_dir_pattern):
    """
    Read one run_many_Pb.sh output set per LHAPDF replica member (each
    matching member_dir_pattern, e.g. 'out/lhapdf/member_*') and combine them
    into a mean +/- standard-deviation band per rapidity. The LHAPDF set's
    ErrorType is "replicas" (see data/prompt-D0-1-109/prompt-D0-1-109.info),
    so the standard deviation across members is the correct uncertainty
    measure -- not a Hessian eigenvector-pair formula.

    Returns {y: (pt_values, mean_values, std_values)}, or {} if no member
    directories are found yet (run_lhapdf_members.sh hasn't been run).
    """
    member_dirs = sorted(glob.glob(member_dir_pattern))
    if not member_dirs:
        return {}

    per_member_results = [
        load_results(os.path.join(member_dir, "files/D0_incl_LHAPDF_An0n_Pb_y*.dat"))
        for member_dir in member_dirs
    ]

    band = {}
    for y in per_member_results[0]:
        pt_values = sorted(pt for pt, _ in per_member_results[0][y])

        cross_sections_by_pt = []
        for pt in pt_values:
            values = []
            for results in per_member_results:
                points = dict(results[y])
                values.append(points[pt])
            cross_sections_by_pt.append(values)

        means = [statistics.mean(values) for values in cross_sections_by_pt]
        stds = [statistics.pstdev(values) for values in cross_sections_by_pt]
        band[y] = (pt_values, means, stds)

    return band


def main():
    results_g1 = load_results("files/D0_incl_KniehlKramer_An0n_G1_Pb_y*.dat")
    results_no_g1 = load_results("files/D0_incl_KniehlKramer_An0n_Pb_y*.dat")
    results_bcfy = load_results("files/D0_incl_BCFY_An0n_Pb_y*.dat")
    lhapdf_band = load_lhapdf_band("files/lhapdf/member_*")
    results_lhapdf = None if lhapdf_band else load_results("files/D0_incl_LHAPDF_An0n_Pb_y*.dat")

    #Plot the results for each rapidity y

    plt.figure(figsize=(7, 5))

    color_cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    colors = {
        y: color_cycle[i % len(color_cycle)]
        for i, y in enumerate(sorted(y for y in results_g1 if y <= 2.0))
    }

    for y in sorted(results_g1):
        if y > 2.0:
            continue
        points = sorted(results_g1[y], key=lambda pair: pair[0])
        pt_values = [pair[0] for pair in points]
        cross_section_values = [pair[1] for pair in points]
        plt.plot(pt_values, cross_section_values, color=colors[y], linestyle="-", label=f"y={y}")

    for y in sorted(results_no_g1):
        if y > 2.0:
            continue
        points = sorted(results_no_g1[y], key=lambda pair: pair[0])
        pt_values = [pair[0] for pair in points]
        cross_section_values = [pair[1] for pair in points]
        plt.plot(pt_values, cross_section_values, color=colors.get(y), linestyle="--")

    for y in sorted(results_bcfy):
        if y > 2.0:
            continue
        points = sorted(results_bcfy[y], key=lambda pair: pair[0])
        pt_values = [pair[0] for pair in points]
        cross_section_values = [pair[1] for pair in points]
        plt.plot(pt_values, cross_section_values, color=colors.get(y), linestyle=":")

    if lhapdf_band:
        for y in sorted(lhapdf_band):
            if y > 2.0:
                continue
            pt_values, means, stds = lhapdf_band[y]
            color = colors.get(y)
            lower = [max(mean - std, 1e-12) for mean, std in zip(means, stds)]
            upper = [mean + std for mean, std in zip(means, stds)]
            plt.fill_between(pt_values, lower, upper, color=color, alpha=0.2, linewidth=0)
            plt.plot(pt_values, means, color=color, linestyle="-.")
    else:
        for y in sorted(results_lhapdf):
            if y > 2.0:
                continue
            points = sorted(results_lhapdf[y], key=lambda pair: pair[0])
            pt_values = [pair[0] for pair in points]
            cross_section_values = [pair[1] for pair in points]
            plt.plot(pt_values, cross_section_values, color=colors.get(y), linestyle="-.")

    plt.yscale("log")
    plt.xlabel(r"$p_{D^0}$ [GeV]")
    plt.ylabel(r"$d\sigma/dy\,dp_T$ [mb/GeV]")
    plt.title("Inclusive $D^0$ photoproduction, Pb+Pb UPC")

    y_legend = plt.legend(loc="upper right")
    plt.gca().add_artist(y_legend)
    style_handles = [
        Line2D([0], [0], color="black", linestyle="-", label="KniehlKramer, G1"),
        Line2D([0], [0], color="black", linestyle="--", label="KniehlKramer, no G1"),
        Line2D([0], [0], color="black", linestyle=":", label="BCFY"),
        Line2D(
            [0], [0], color="black", linestyle="-.",
            label="LHAPDF (mean +/- std over replicas)" if lhapdf_band else "LHAPDF (member 0, no errors)",
        ),
    ]
    plt.legend(handles=style_handles, loc="lower left")

    plt.tight_layout()
    plt.savefig("plots/D0_incl_dsigma_dy_dpt_PbPb_g1.png", dpi=150)


if __name__ == "__main__":
    main()

