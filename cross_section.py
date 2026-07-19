import glob
import math
import matplotlib.pyplot as plt
from scipy.integrate import simpson


pi = math.pi
alpha_em = 1 / 137
e_charm_squared = 4 / 9  # (2/3)^2, electric charge of the charm quark
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


def main():
    filenames = sorted(glob.glob("files/D0_incl_KniehlKramer_An0n_Gamma1_Pb_y*.dat"))

    # results[y] = list of (pt, cross_section) pairs
    results = {}

    for filename in filenames:
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


    #Plot the results for each rapidity y
    
    plt.figure(figsize=(7, 5))

    for y in sorted(results):
        if y > 2.0:
            continue
        points = sorted(results[y], key=lambda pair: pair[0])
        pt_values = [pair[0] for pair in points]
        cross_section_values = [pair[1] for pair in points]
        plt.plot(pt_values, cross_section_values, label=f"y={y}")

    plt.yscale("log")
    plt.xlabel(r"$p_{D^0}$ [GeV]")
    plt.ylabel(r"$d\sigma/dy\,dp_T$ [mb/GeV]")
    plt.title("Inclusive $D^0$ photoproduction, Pb+Pb UPC")
    plt.legend()
    plt.tight_layout()
    plt.savefig("plots/D0_incl_dsigma_dy_dpt_PbPb_gamma1.png", dpi=150)


