import glob
import math

import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import simpson

pi = math.pi
alfae = 1/137
eq = 4/9  # e_charm=(2/3)^2
Nc = 3

# Conversion factors
FMGEV = 5.068
GEVSQRTONB = 1.0e7/(FMGEV*FMGEV)
GEVSQRTOMB = GEVSQRTONB*1e-6

facA = alfae*eq*Nc/(2*pi)**4  # nucleus target -- b (glauber_mve_<b>) is already GeV^-1, no FMGEV needed here
# In the factor above, the (2*pi)^4 comes from the prefactor and the integration over the dipole in momentum space

A = 208
ATA = 30.756  # dimensionless

# Proton-target normalization (facp) needs `radius` defined and a proton result
# file, neither of which exist yet -- add back once proton runs are available.

## AA collisions: b is in fm (the glauber_mve_<b> file tag), dsigma_dy is the
## raw per-(b, pT, y) dipole output from run_many_Pb.sh. Collapse the b
## dimension with Simpson's rule, weighted by b (2D polar radial measure):
##   integral(pT, y) = int b * dsigma_dy(b, pT, y) db
## the remaining angular 2*pi factor is applied below alongside the pT one.

def read_y(fname):
    with open(fname) as f:
        for line in f:
            if "fixed rapidity y" in line:
                return float(line.split(":")[-1])
    raise ValueError(f"no rapidity header found in {fname}")


def main():
    files = sorted(glob.glob("files/D0_incl_KniehlKramer_An0n_Pb_y*.dat"))

    rows = []
    for fname in files:
        y = read_y(fname)
        df = pd.read_csv(fname, comment="#", sep=r"\s+", names=["b", "pt", "dsigma_dy"])
        for pt, group in df.groupby("pt"):
            group = group.sort_values("b")
            integral = simpson(group["b"]*group["dsigma_dy"], x=group["b"])
            rows.append({"y": y, "pt": pt, "res": integral})

    result_A = pd.DataFrame(rows)
    result_A["res"] = (2*pi)*result_A["res"]*facA*(2*pi)*result_A["pt"]*GEVSQRTOMB
    # The factor of (2*pi)^2 comes from d^2p -> dp*p*(2pi) and the 2D impact parameter b-integral (angular part)

    plt.figure(figsize=(7, 5))
    for y, group in result_A.groupby("y"):
        group = group.sort_values("pt")
        plt.plot(group["pt"], group["res"], label=f"y={y}")

    plt.yscale("log")
    plt.xlabel(r"$p_{D^0}$ [GeV]")
    plt.ylabel(r"$d\sigma/dy\,dp_T$ [mb/GeV]")
    plt.title("Inclusive $D^0$ photoproduction, Pb+Pb UPC")
    plt.legend()
    plt.tight_layout()
    plt.savefig("plots/D0_incl_dsigma_dy_dpt_PbPb.png", dpi=150)


if __name__ == "__main__":
    main()
