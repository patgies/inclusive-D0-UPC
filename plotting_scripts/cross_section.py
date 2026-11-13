import glob
import math
import os
from scipy.integrate import simpson, trapezoid


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
SIGMA0_HALF = 16.36   # mb, proton normalization of the MVe dipole. It replaces the b_d integral.

# The Glauber dipole files stop at b_d = 32 GeV^-1, which leaves out the edge of the nucleus
# (about 20 nucleons). As in arXiv:2503.16108, the edge is dilute and is added as
# A_edge * (proton cross section with the same flux). EDGE=0 leaves it out.
EDGE = os.environ.get("EDGE", "1") != "0"
R_WS, A_WS = 6.49, 0.54   # Woods-Saxon radius and diffuseness [fm]
HBARC = 0.197327


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


def load_proton_results(pattern):
    # Reads the files "pD0  dsigma" of run_proton.sh. Returns {y: [(pT, cross section), ...]} in mb/GeV.
    results = {}
    for filename in sorted(glob.glob(pattern)):
        points = []
        with open(filename) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                pt, raw = (float(x) for x in line.split())
                points.append((pt, raw * factor_A * SIGMA0_HALF * (2 * pi) * pt))
        results[read_rapidity(filename)] = points
    return results


def nucleons_beyond(b_max):
    # Number of nucleons of the Woods-Saxon nucleus at transverse distance > b_max [GeV^-1]
    import numpy as np
    s = np.linspace(0.0, 25.0, 2501)
    z = np.linspace(-40.0, 40.0, 4001)
    TA = np.array([trapezoid(1.0 / (1.0 + np.exp((np.sqrt(z**2 + x**2) - R_WS) / A_WS)), z) for x in s])
    weight = 2 * pi * s * TA
    inside = s <= b_max * HBARC
    return A * (1.0 - trapezoid(weight[inside], s[inside]) / trapezoid(weight, s))


def add_edge(nucleus, nucleus_pattern, proton_pattern):
    # nucleus + A_edge * proton, point by point. Stops if a proton point is missing.
    if not EDGE or not nucleus:
        return nucleus
    b_max = max(float(line.split()[0]) for f in glob.glob(nucleus_pattern)
                for line in open(f) if line.strip() and line[0] != "#")
    a_edge = nucleons_beyond(b_max)
    proton = load_proton_results(proton_pattern)
    out = {}
    for y, points in nucleus.items():
        p = dict(proton.get(y, []))
        missing = [pt for pt, _ in points if pt not in p]
        if missing:
            raise SystemExit(f"Edge of the nucleus: no proton result for y = {y:g}, pT = {missing[0]:g} in "
                             f"{proton_pattern} -- run run_scripts/run_proton.sh (TARGET=AA), or EDGE=0.")
        out[y] = [(pt, value + a_edge * p[pt]) for pt, value in points]
    return out


def proton_dir(base, frag, scale="1.0"):
    # proton reference of output/<CHANNEL>/central_values/<frag>/ (scale 1) or of scale_variation/
    return f"{base}/proton_baseline/{frag}" if scale == "1.0" else f"{base}/proton_baseline/factor_{scale}/{frag}"
