"""Supplementary Figure 2 — choosing the number of clonal clusters.

A  the selection curve per mouse: number of clusters against the correlation
   cutoff, with the chosen cutoff marked. The cutoff is where the curve
   flattens, i.e. where adding stringency stops splitting real structure.
B  the resulting clusters, raw barcode trajectories behind the LOESS consensus
   that every downstream analysis actually uses.

Run:  python figures/supplementary/figS2_cluster_selection.py
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _common import CLONE_COLORS, COLONISED, io, style, time_axis, to_days
from succession.config import DATA

style.apply()


def main():
    """Returns {'A': fig, 'B': fig}."""
    figs = {}
    fig, axes = plt.subplots(2, 4, figsize=(15, 6.4))
    for ax, m in zip(axes.flat, COLONISED):
        sel = pd.read_csv(DATA / "barcodes" / "clusters" / f"{m}_threshold_selection.csv")
        ax.plot(sel["cutoff"], sel["n_clusters"], "-o", ms=2.5, color="#2C3E50", lw=1.1)
        n_final = io.load_clone_loess(m)["cluster"].nunique()
        hit = sel[sel["n_clusters"] == n_final]
        if len(hit):
            ax.axvline(hit["cutoff"].iloc[0], color="#C0392B", ls="--", lw=1.0)
        ax.axhline(n_final, color="#C0392B", ls=":", lw=0.8)
        ax.set_title(f"{m}  ({n_final} clusters)", fontweight="bold")
        ax.set_xlabel("Correlation cutoff")
    for ax in axes[:, 0]:
        ax.set_ylabel("Number of clusters")
    fig.suptitle("Supplementary Figure 2A — cluster-number selection", y=0.99)
    fig.tight_layout()
    figs["A"] = style.save(fig, "figS2A_cluster_selection")

    fig, axes = plt.subplots(2, 4, figsize=(15, 6.4), sharey=True)
    for ax, m in zip(axes.flat, COLONISED):
        raw = pd.read_csv(DATA / "barcodes" / "clusters" / f"{m}_clustered_series.csv")
        lo = io.load_clone_loess(m)
        for cl, grp in raw.groupby("cluster"):
            colour = CLONE_COLORS.get(f"C{cl}", "#999999")
            for _, b in grp.groupby("ID"):
                b = b.sort_values("time")
                # rasterized: one line per barcode, far too many for vector EPS
                ax.plot(to_days(b["time"].to_numpy(), group="colonised"),
                        b["frequency"], color=colour, lw=0.25, alpha=0.25,
                        rasterized=True)
        for clone, grp in lo.groupby("clone"):
            ax.plot(to_days(grp["index"].to_numpy(), group="colonised"),
                    10 ** grp["log10_freq"],
                    color=CLONE_COLORS.get(clone, "#999999"), lw=1.6)
        ax.set_yscale("log"); ax.set_ylim(1e-7, 1e-1)
        ax.set_title(m, fontweight="bold")
        time_axis(ax, m)
    for ax in axes[:, 0]:
        ax.set_ylabel("Frequency")
    fig.suptitle("Supplementary Figure 2B — cluster trajectories, raw behind LOESS", y=0.99)
    fig.tight_layout()
    figs["B"] = style.save(fig, "figS2B_cluster_trajectories")
    return figs


if __name__ == "__main__":
    print("Supplementary Figure 2")
    main()
