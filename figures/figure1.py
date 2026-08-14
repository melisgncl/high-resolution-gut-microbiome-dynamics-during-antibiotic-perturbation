"""Figure 1 — colonisation dynamics of barcoded E. coli.

B    CFU per gram
C-D  barcode frequency trajectories, top 1000 coloured
E-F  dominant clonal clusters
G-H  16S community composition at family level
I    co-clustering of clonal lineages with bacterial families

Panel A is an experimental schematic and has no code.

Run:  python figures/figure1.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from succession import io, style                                   # noqa: E402
from succession.config import (COHORT_1, COHORT_2, COLONISED,       # noqa: E402
                               CLONE_COLORS, FAMILY_COLORS, FAMILY_GATE,
                               LEGEND_FAMILIES, MOUSE_COLORS,
                               N_COLOUR_BARCODES, PSEUDOCOUNT_16S_PLOT)
from succession.timeaxis import to_days                            # noqa: E402

style.apply()
CFU_COLORS = {"m1": "#1F5C8B", "m2": "#2E86C1", "m3": "#5DADE2", "m4": "#AED6F1",
              "m5": "#C2185B", "m6": "#D81B60", "m7": "#8E24AA", "m8": "#5E35B1"}


#: Sampling days for a colonised mouse: 3 h, 6 h, 12 h, then daily.
def _sample_days(upto):
    return [3 / 24, 6 / 24, 12 / 24] + [float(d) for d in range(1, upto + 1)]


def _time_axis(ax, upto=16):
    """A tick at every sample, labelled only where the published panel labels."""
    labelled = {3 / 24: "3h", 2.0: "2d", 6.0: "6d", 10.0: "10d",
                float(upto): f"{upto}d"}
    days = _sample_days(upto)
    ax.set_xticks(days)
    ax.set_xticklabels([labelled.get(d, "") for d in days],
                       rotation=45, ha="right")


# ── B ────────────────────────────────────────────────────────────────────────
def panel_b():
    cfu = io.load_cfu()
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    for m in COLONISED:
        d = cfu[cfu["mouse"] == m].sort_values("hours")
        ax.plot(d["hours"] / 24, d["cfu"], color=CFU_COLORS[m], lw=1.2, label=m)
    ax.set_yscale("log")
    ax.set_xlabel("Time post-colonisation (days)")
    ax.set_ylabel("CFU per gram")
    ax.set_title("B   E. coli load")
    ax.legend(ncol=2, loc="lower right")
    return style.save(fig, "fig1B_cfu")


# ── C-D ──────────────────────────────────────────────────────────────────────
def _polylines(sub):
    """(x, y, colour) per barcode, sorted within each trajectory.

    Rows must be grouped by ID and ordered in time before drawing, otherwise
    matplotlib joins points in table order and the panel fills with spurious
    diagonals.
    """
    sub = sub.sort_values(["ID", "index"])
    ids = sub["ID"].to_numpy()
    day = to_days(sub["index"].to_numpy(), group="colonised")
    freq = sub["freq"].to_numpy()
    hexes = sub["hex"].to_numpy()
    edges = np.flatnonzero(np.r_[True, ids[1:] != ids[:-1], True])
    return [(day[a:b], freq[a:b], hexes[a])
            for a, b in zip(edges[:-1], edges[1:]) if b - a > 1]


def panels_cd():
    """The caption colours the 1000 most abundant barcodes.

    The stored `hex` assignment is broader (~4000 per mouse, 31842 overall),
    which floods the panel and hides the grey background the published figure
    shows, so barcodes are ranked by peak frequency and the top N kept. Each
    keeps its stored colour, so a lineage looks the same in every mouse.
    """
    from matplotlib.collections import LineCollection

    bc = io.load_barcode_trajectories()
    bc = bc[bc["freq"] > 0]
    peak = (bc[bc["hex"] != "#cccccc"]
            .groupby(["mouse", "ID"])["freq"].max().reset_index())
    top = (peak.sort_values("freq", ascending=False)
           .groupby("mouse").head(N_COLOUR_BARCODES)[["mouse", "ID"]])
    top["coloured"] = True
    bc = bc.merge(top, on=["mouse", "ID"], how="left")
    bc["coloured"] = bc["coloured"].notna()
    print(f"  coloured {int(bc.groupby('mouse')['coloured'].apply(lambda s: s.any()).sum())} "
          f"mice, {N_COLOUR_BARCODES} barcodes each")

    fig, axes = plt.subplots(2, 4, figsize=(15, 6.4), sharey=True)
    for ax, m in zip(axes.flat, COLONISED):
        d = bc[bc["mouse"] == m]
        grey = _polylines(d[~d["coloured"]])
        col = _polylines(d[d["coloured"]])
        ax.add_collection(LineCollection([np.column_stack(t[:2]) for t in grey],
                                         colors="#cccccc", linewidths=0.15,
                                         alpha=0.5, zorder=1,
                                         rasterized=True))
        ax.add_collection(LineCollection([np.column_stack(t[:2]) for t in col],
                                         colors=[t[2] for t in col],
                                         linewidths=0.45, alpha=0.9, zorder=2,
                                         rasterized=True))
        ax.set_yscale("log")
        ax.set_ylim(1e-7, 1)
        ax.set_xlim(0, 16.5 if m in COHORT_1 else 15.5)
        ax.set_title(m, fontweight="bold")
        _time_axis(ax, 16 if m in COHORT_1 else 15)
    for ax in axes[:, 0]:
        ax.set_ylabel("Barcode frequency")
    fig.suptitle("C / D   Barcode composition over time", y=0.98)
    fig.tight_layout()
    style.save(fig, "fig1CD_barcode_dynamics")
    del bc
    return fig


# ── E-F ──────────────────────────────────────────────────────────────────────
def panels_ef():
    fig, axes = plt.subplots(2, 4, figsize=(15, 6.4), sharey=True)
    for ax, m in zip(axes.flat, COLONISED):
        lo = io.load_clone_loess(m)          # applies the m2/m8 C1<->C2 swap
        for clone, grp in lo.groupby("clone"):
            ax.plot(to_days(grp["index"].to_numpy(), group="colonised"),
                    10 ** grp["log10_freq"],
                    color=CLONE_COLORS.get(clone, "#999999"), lw=1.1)
        ax.set_yscale("log")
        ax.set_ylim(1e-7, 1e-1)
        ax.set_title(m, fontweight="bold")
        _time_axis(ax, 16 if m in COHORT_1 else 15)
    for ax in axes[:, 0]:
        ax.set_ylabel("Clone frequency")
    fig.suptitle("E / F   Dominant E. coli clonal clusters", y=0.98)
    fig.tight_layout()
    return style.save(fig, "fig1EF_clonal_clusters")


# ── G-H ──────────────────────────────────────────────────────────────────────
def _composition_panel(ax, mouse, group, upto):
    """One mouse's family composition.

    Two rules, both reconstructed from the published panel:
      * zeros get a pseudocount so a vanishing family drops to the axis floor
        and returns, rather than leaving a gap;
      * a run of zeros becomes NaN so the line breaks instead of being bridged
        by a straight segment;
      * only families whose own mean in this mouse exceeds FAMILY_GATE are drawn.
    """
    fam = io.load_family(mouse)
    means = fam.groupby("Family")["abundance"].mean()
    keep = [f for f in LEGEND_FAMILIES if means.get(f, 0) > FAMILY_GATE]

    for f in keep:
        d = fam[fam["Family"] == f].sort_values("Time")
        a = d["abundance"].to_numpy()
        touching = (a > 0) | (np.r_[0, a[:-1]] > 0) | (np.r_[a[1:], 0] > 0)
        y = np.where(touching, a + PSEUDOCOUNT_16S_PLOT, np.nan)
        ax.plot(to_days(d["Time"].to_numpy(), group=group), y,
                color=FAMILY_COLORS[f], lw=1.2,
                label=f if ax.get_legend_handles_labels()[1].count(f) == 0 else None)
    ax.set_yscale("log")
    ax.set_ylim(1e-6, 1)
    ax.set_title(mouse, fontweight="bold")


def panels_gh():
    """Returns {'G': fig, 'H': fig} — one figure per cohort."""
    figs = {}
    for tag, mice, upto in (("G", COHORT_1, 16), ("H", COHORT_2, 15)):
        fig, axes = plt.subplots(1, 4, figsize=(15, 3.6), sharey=True)
        for ax, m in zip(axes, mice):
            _composition_panel(ax, m, "colonised", upto)
            _time_axis(ax, upto)
        axes[0].set_ylabel("Family composition")
        handles = [plt.Line2D([], [], color=FAMILY_COLORS[f], lw=2, label=f)
                   for f in LEGEND_FAMILIES]
        fig.legend(handles=handles, loc="lower center", ncol=6,
                   bbox_to_anchor=(0.5, -0.22))
        fig.suptitle(f"{tag}   16S community composition", y=1.0)
        fig.tight_layout()
        figs[tag] = style.save(fig, f"fig1{tag}_16S_cohort{1 if tag == 'G' else 2}")
    return figs


# ── I ────────────────────────────────────────────────────────────────────────
def panel_i():
    """UPGMA on the shape-based distances between clone and family trajectories."""
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    for ax, m in zip(axes.flat, COLONISED):
        d = io.load_sbd_distance(m)
        labels = sorted(set(d["series1"]) | set(d["series2"]))
        pos = {l: i for i, l in enumerate(labels)}
        M = np.zeros((len(labels), len(labels)))
        for _, r in d.iterrows():
            i, j = pos[r["series1"]], pos[r["series2"]]
            M[i, j] = M[j, i] = r["dist"]
        Z = linkage(squareform(M, checks=False), method="average")
        dendrogram(Z, labels=labels, orientation="right", ax=ax,
                   color_threshold=0, above_threshold_color="#444444",
                   leaf_font_size=6)
        ax.set_title(m, fontweight="bold")
        ax.set_xticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    fig.suptitle("I   Co-clustering of clonal lineages with bacterial families")
    fig.tight_layout()
    return style.save(fig, "fig1I_coclustering")


if __name__ == "__main__":
    print("Figure 1")
    panel_b()
    panels_ef()
    panels_gh()
    panel_i()
    panels_cd()          # last: loads 10.4 M rows
