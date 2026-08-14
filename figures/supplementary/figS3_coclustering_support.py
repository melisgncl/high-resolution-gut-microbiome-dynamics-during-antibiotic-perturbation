"""Supplementary Figure 3 — per-mouse co-clustering, with cluster membership.

Each dendrogram is UPGMA on shape-based distances between clonal-cluster and
family trajectories, with leaves coloured by the co-cluster they were assigned
to. Whether a clone sits with Paenibacillaceae or with Enterobacteriaceae is the
per-animal version of the association Figure 3C-D quantifies.

The published panel also carries pvclust bootstrap support. Those resampling
values are not among the shipped tables, so support is not annotated here; the
topology and membership are.

Run:  python figures/supplementary/figS3_coclustering_support.py
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

from _common import COLONISED, io, style
from succession.config import DATA

style.apply()


def main():
    palette = plt.get_cmap("tab10")
    fig, axes = plt.subplots(2, 4, figsize=(16, 8.4))
    for ax, m in zip(axes.flat, COLONISED):
        d = io.load_sbd_distance(m)
        labels = sorted(set(d["series1"]) | set(d["series2"]))
        pos = {l: i for i, l in enumerate(labels)}
        M = np.zeros((len(labels), len(labels)))
        for _, r in d.iterrows():
            i, j = pos[r["series1"]], pos[r["series2"]]
            M[i, j] = M[j, i] = r["dist"]
        Z = linkage(squareform(M, checks=False), method="average")

        memb = {}
        path = DATA / "coclustering" / f"{m}_sbd_clusters.csv"
        if path.exists():
            cl = pd.read_csv(path)
            memb = dict(zip(cl["series"], cl["sbd_cluster"]))

        dendrogram(Z, labels=labels, orientation="right", ax=ax,
                   color_threshold=0, above_threshold_color="#555555",
                   leaf_font_size=6)
        for tick in ax.get_yticklabels():
            k = memb.get(tick.get_text())
            if k is not None:
                tick.set_color(palette((int(k) - 1) % 10))
            if not tick.get_text().startswith("C"):
                tick.set_fontstyle("italic")
        ax.set_title(m, fontweight="bold")
        ax.set_xticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    fig.suptitle("Supplementary Figure 3 — co-clustering of clonal lineages with "
                 "bacterial families (leaf colour = co-cluster)", y=0.99)
    fig.tight_layout()
    return style.save(fig, "figS3_coclustering")


if __name__ == "__main__":
    print("Supplementary Figure 3")
    main()
