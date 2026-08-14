"""Supplementary Figure 6 — co-clustering pooled across all mice.

A  the global shape-based distance matrix over every clonal cluster from every
   mouse, ordered by hierarchical clustering. Blocks that mix mice are lineages
   behaving the same way in different animals.
B  cluster-overlap coefficients between mice.

Run:  python figures/supplementary/figS6_global_coclustering.py
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import squareform

from _common import MOUSE_COLORS, style
from succession.config import DATA

style.apply()


def _matrix(path, index_col):
    df = pd.read_csv(path)
    labels = df[index_col].tolist()
    M = df[labels].to_numpy(float)
    M = (M + M.T) / 2.0
    np.fill_diagonal(M, 0.0)
    return labels, M


def main():
    labels, M = _matrix(DATA / "coclustering" / "all_sbd_dist.csv", "series")
    order = leaves_list(linkage(squareform(M, checks=False), method="average"))
    lab = [labels[i] for i in order]
    Mo = M[np.ix_(order, order)]

    fig, (ax, axb) = plt.subplots(1, 2, figsize=(14, 6),
                                  gridspec_kw={"width_ratios": [1.25, 1]})
    im = ax.imshow(Mo, cmap="viridis_r", aspect="auto")
    fig.colorbar(im, ax=ax, shrink=0.8, label="shape-based distance")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title("A   Global SBD distance, all clonal clusters")
    # a colour strip showing which mouse each row came from
    mice = [l.split(".")[0] for l in lab]
    strip = np.array([[list(MOUSE_COLORS).index(mm) if mm in MOUSE_COLORS else -1
                       for mm in mice]])
    ax.imshow(strip, extent=(0, len(lab), len(lab), len(lab) + len(lab) * 0.02),
              cmap="tab10", aspect="auto")

    ov = pd.read_csv(DATA / "coclustering" / "similarity_all_overlap.csv")
    cols = [c for c in ov.columns if c != "cluster"]
    O = ov[cols].to_numpy(float)
    im2 = axb.imshow(O, cmap="magma", aspect="auto", vmin=0, vmax=1)
    fig.colorbar(im2, ax=axb, shrink=0.8, label="overlap coefficient")
    axb.set_xticks([]); axb.set_yticks([])
    axb.set_title("B   Cluster overlap between mice")
    fig.suptitle("Supplementary Figure 6 — co-clustering pooled across mice", y=0.99)
    fig.tight_layout()
    return style.save(fig, "figS6_global_coclustering")


if __name__ == "__main__":
    print("Supplementary Figure 6")
    main()
