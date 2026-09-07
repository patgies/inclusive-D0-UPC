#!/usr/bin/env python3
# Plots the DGLAP-evolved Braaten-Cheung-Fleming-Yuan (BCFY) c->D0
# fragmentation function D(z, Q^2) at a few different scales mu, one panel
# per reference pT, to visualize how the shape shifts (see
# plot_kk_ff_shapes.py for the Kniehl-Kramer analogue).
#
# usage: python3 python/plot_bcfy_ff_shapes.py

import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("out/bcfy_ff_shapes/bcfy_ff.csv")
pTvals = sorted(df["pT"].unique())

fig, axes = plt.subplots(1, len(pTvals), figsize=(7 * len(pTvals), 5))

for ax, pT in zip(axes, pTvals):
    sub = df[df["pT"] == pT]
    for label, group in sub.groupby("label", sort=False):
        style = "--" if label.startswith("IC") else "-"
        ax.plot(group["z"], group["D"], style, label=label)
    ax.set_xlabel("z")
    ax.set_ylabel(r"$D(z, \mu^2)$")
    ax.set_title(f"pT = {pT:g} GeV")
    ax.legend(frameon=False)
fig.suptitle("LO BCFY FF with DGLAP, $mt=\sqrt{p^2+1.5^2}$")
fig.savefig("plots/bcfy_FragmFunc.pdf", dpi=150, bbox_inches="tight")
