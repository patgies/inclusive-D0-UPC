import argparse
import numpy as np

# Converts the flux tables of P. Paakkinen (arXiv:2404.09731) to the format of FLUX_MODEL=TABLE,
# lines "z_gamma  f  z_gamma*f" with f = dN/dz_gamma.
# usage: python3 src/make_flux_table.py

TABLE_DIR = "input/Starlight_photon_flux"
COLUMN = {"AnAn": 3, "An0n": 4}   # columns: j  y  y^(1/4)  log_AnAn  log_An0n  (y = z_gamma)


def flux(z_gamma, x_nodes, logf_nodes):
    """f(z_gamma) from the table: interpolation in z_gamma^(1/4) of Eqs. 34-37 of arXiv:2404.09731.
    Below the first table point z_gamma*f(z_gamma) = A ln(c/z_gamma), fixed by the first two points (as in diffractive-D0-UPC)."""
    z0, z1 = x_nodes[0] ** 4, x_nodes[1] ** 4
    if z_gamma < z0:
        g0, g1 = z0 * np.exp(logf_nodes[0]), z1 * np.exp(logf_nodes[1])
        A = (g0 - g1) / np.log(z1 / z0)
        logc = g0 / A + np.log(z0)
        return A * (logc - np.log(z_gamma)) / z_gamma
    d = z_gamma ** 0.25 - x_nodes
    if np.any(np.abs(d) < 1e-14):
        return np.exp(logf_nodes[np.argmin(np.abs(d))])
    n = len(x_nodes) - 1
    w = np.array([(-1.0 if j % 2 else 1.0) * (0.5 if j in (0, n) else 1.0) for j in range(n + 1)]) / d
    return np.exp(np.sum(w * logf_nodes) / np.sum(w))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zmin", type=float, default=1e-6)
    ap.add_argument("--n", type=int, default=200)
    args = ap.parse_args()

    for kind in ("PL", "WS"):
        source = f"{TABLE_DIR}/log-flux-tbl-{kind}.dta"
        table = np.loadtxt(source)
        with open(source) as fh:
            first_line = fh.readline().strip("# \n")
        for channel, col in COLUMN.items():
            x_nodes, logf_nodes = table[:, 2], table[:, col]
            z_gamma = np.logspace(np.log10(args.zmin), np.log10(table[-1, 1]), args.n)
            z_gamma[-1] = table[-1, 1]
            f = np.array([flux(z_gamma_i, x_nodes, logf_nodes) for z_gamma_i in z_gamma])
            out = f"{TABLE_DIR}/flux_{kind}_{channel}.dat"
            header = "\n".join([
                "Photon flux f(z_gamma) = dN/dz_gamma, z_gamma = omega/E_beam = photon energy / beam energy per nucleon",
                f"{channel} flux of P. Paakkinen, from {source} with src/make_flux_table.py. Please cite arXiv:2404.09731.",
                first_line,
                f"below z_gamma = {table[0, 1]:g} (first point of the table): z_gamma*f = A ln(c/z_gamma), fixed by the first two points",
                "z_gamma  f(z_gamma)  z_gamma*f(z_gamma)",
            ])
            np.savetxt(out, np.c_[z_gamma, f, z_gamma * f], fmt="%.10e", header=header, comments="# ")
            print(f"Saved: {out}")


if __name__ == "__main__":
    main()
