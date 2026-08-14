"""Supplementary Figure 1 — barcode dynamics on a linear frequency scale.

The same trajectories as Figure 1C-D. A log axis shows the rare majority; a
linear axis shows which lineages actually carry the population, and makes the
handful of clones reaching percent-level frequency obvious.

Run:  python figures/supplementary/figS1_barcode_linear.py
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from _common import COHORT_1, COLONISED, io, style, time_axis, to_days

style.apply()


def main(top_n: int = 200):
    bc = io.load_barcode_trajectories()
    bc = bc[bc["freq"] > 0]
    peak = bc.groupby(["mouse", "ID"])["freq"].max().reset_index()
    keep = (peak.sort_values("freq", ascending=False)
            .groupby("mouse").head(top_n)[["mouse", "ID"]])
    keep["keep"] = True
    bc = bc.merge(keep, on=["mouse", "ID"], how="left")
    bc = bc[bc["keep"].notna()]

    fig, axes = plt.subplots(2, 4, figsize=(15, 6.4), sharey=True)
    for ax, m in zip(axes.flat, COLONISED):
        d = bc[bc["mouse"] == m].sort_values(["ID", "index"])
        ids = d["ID"].to_numpy()
        x = to_days(d["index"].to_numpy(), group="colonised")
        y = d["freq"].to_numpy()
        hexes = d["hex"].to_numpy()
        cut = np.flatnonzero(np.r_[True, ids[1:] != ids[:-1], True])
        segs = [np.column_stack((x[a:b], y[a:b])) for a, b in zip(cut[:-1], cut[1:]) if b - a > 1]
        cols = [hexes[a] for a, b in zip(cut[:-1], cut[1:]) if b - a > 1]
        # rasterized: thousands of segments would make the EPS enormous
        ax.add_collection(LineCollection(segs, colors=cols, linewidths=0.6,
                                         alpha=0.85, rasterized=True))
        ax.set_xlim(0, 16.5 if m in COHORT_1 else 15.5)
        ax.set_ylim(0, max(0.05, float(y.max()) * 1.05))
        ax.set_title(m, fontweight="bold")
        time_axis(ax, m)
    for ax in axes[:, 0]:
        ax.set_ylabel("Barcode frequency (linear)")
    fig.suptitle(f"Supplementary Figure 1 — barcode dynamics, linear scale "
                 f"(top {top_n} per mouse)", y=0.99)
    fig.tight_layout()
    return style.save(fig, "figS1_barcode_linear")


if __name__ == "__main__":
    print("Supplementary Figure 1")
    main()
