#!/usr/bin/env python3
# Plots the three c -> D0 fragmentation functions used in src/main.cpp
# (BCFY, Kniehl & Kramer, LHAPDF) as a function of z.

import numpy as np
import matplotlib.pyplot as plt

# Parameters, matching src/main.cpp
r = 0.1
N_kk, eps_kk = 0.694, 0.101
Q = 1.5                 # charm mass, param.m in main.cpp
lhapdf_flavor = 4       # PDG id for the charm quark
lhapdf_dir = "../data/prompt-D0-1-109"
lhapdf_set = "prompt-D0-1-109"
lhapdf_num_replicas = 100   # members 0001-0100; member 0000 is the central value


def bcfy_dp(z, r):
    z = np.asarray(z, dtype=float)
    denom = 1.0 - (1.0 - r) * z
    factor = r * z * (1.0 - z)**2 / denom**6
    poly = (6.0 - 18.0 * (1.0 - 2.0 * r) * z
            + (21.0 - 74.0 * r + 68.0 * r**2) * z**2
            - 2.0 * (1.0 - r) * (6.0 - 19.0 * r + 18.0 * r**2) * z**3
            + 3.0 * (1.0 - r)**2 * (1.0 - 2.0 * r + 2.0 * r**2) * z**4)
    val = factor * poly
    return np.where((z <= 0.0) | (z >= 1.0), 0.0, val)


def bcfy_dv(z, r):
    z = np.asarray(z, dtype=float)
    denom = 1.0 - (1.0 - r) * z
    factor = 3.0 * (r * z * (1.0 - z)**2) / denom**6
    poly = (2.0 - 2.0 * (3.0 - 2.0 * r) * z
            + 3.0 * (3.0 - 2.0 * r + 4.0 * r**2) * z**2
            - 2.0 * (1.0 - r) * (4.0 - r + 2.0 * r**2) * z**3
            + (1.0 - r)**2 * (3.0 - 2.0 * r + 2.0 * r**2) * z**4)
    val = factor * poly
    return np.where((z <= 0.0) | (z >= 1.0), 0.0, val)


def dc_to_d0(z, r):
    mD, mDstar = 1.8648, 2.0067
    mass_ratio = mDstar / mD

    DP_part = 0.168 * bcfy_dp(z, r)

    z_scaled = mass_ratio * z
    step = np.where((mD / mDstar) - np.asarray(z, dtype=float) >= 0.0, 1.0, 0.0)
    DV_part = 0.39 * step * bcfy_dv(z_scaled, r) * mass_ratio

    return DP_part + DV_part


def d_kniehl_kramer(z, N, eps):
    z = np.asarray(z, dtype=float)
    one_minus_z = 1.0 - z
    denom = one_minus_z**2 + eps * z
    val = N * z * one_minus_z**2 / denom**2
    return np.where((z <= 0.0) | (z >= 1.0), 0.0, val)


def parse_lhapdf_grid(path, pdg_flavor, Q):
    with open(path) as f:
        lines = f.read().splitlines()

    i = 0
    while lines[i].strip() != "---":
        i += 1
    i += 1

    while i < len(lines) and lines[i].strip():
        xgrid = np.array([float(v) for v in lines[i].split()])
        qgrid = np.array([float(v) for v in lines[i + 1].split()])
        flavors = [int(v) for v in lines[i + 2].split()]
        col = flavors.index(pdg_flavor)

        nx, nq = len(xgrid), len(qgrid)
        values = np.empty((nx, nq))
        row = i + 3
        for ix in range(nx):
            for iq in range(nq):
                values[ix, iq] = float(lines[row].split()[col])
                row += 1
        i = row + 1  # skip the "---" line ending this subgrid

        if qgrid[0] <= Q <= qgrid[-1]:
            iq_hi = np.searchsorted(qgrid, Q)
            iq_hi = min(max(iq_hi, 1), nq - 1)
            iq_lo = iq_hi - 1
            t = (Q - qgrid[iq_lo]) / (qgrid[iq_hi] - qgrid[iq_lo])
            xf = values[:, iq_lo] * (1.0 - t) + values[:, iq_hi] * t
            return xgrid, xf / xgrid   # x*f(x,Q) -> D(z,Q)

    raise RuntimeError(f"Q={Q} outside all subgrids in {path}")


def lhapdf_member_path(member):
    return f"{lhapdf_dir}/{lhapdf_set}_{member:04d}.dat"


z = np.linspace(0.001, 0.999, 500)

plt.plot(z, dc_to_d0(z, r), label="BCFY")
plt.plot(z, d_kniehl_kramer(z, N_kk, eps_kk), label="Kniehl & Kramer")

z_lhapdf, D_central = parse_lhapdf_grid(lhapdf_member_path(0), lhapdf_flavor, Q)
mask = (z_lhapdf >= z.min()) & (z_lhapdf <= z.max())

# 68% CL band from the 16th/84th percentiles of the replica members
# (ErrorType: replicas in prompt-D0-1-109.info), not assumed Gaussian.
replicas = np.array([
    parse_lhapdf_grid(lhapdf_member_path(m), lhapdf_flavor, Q)[1]
    for m in range(1, lhapdf_num_replicas + 1)
])
D_lo, D_hi = np.percentile(replicas, [16, 84], axis=0)

line, = plt.plot(z_lhapdf[mask], D_central[mask], label="LHAPDF")
plt.fill_between(z_lhapdf[mask], D_lo[mask], D_hi[mask],
                  color=line.get_color(), alpha=0.3, label="LHAPDF 68% CL")

plt.xlabel("z")
plt.ylabel("D(z)")
plt.title("c -> D0 fragmentation functions")
plt.legend()
plt.tight_layout()
plt.savefig("fragmentation_functions.png", dpi=150)
plt.show()
