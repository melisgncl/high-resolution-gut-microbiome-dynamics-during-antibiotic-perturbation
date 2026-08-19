"""Figure 2 — the ecological shift tracks community diversity.

A  diversity and the Paenibacillaceae transition, per mouse
B  ridgeline distributions of J[i<-j] over time
C  mean negative J against 16S Hill 1D, pooled and per mouse

Run:  python figures/figure2.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from succession import anchors, diversity, io, jacobian, stats, style   # noqa: E402
from succession.config import (BARCODE_3H_EMPTY, COHORT_1, COLONISED,    # noqa: E402
                               MOUSE_COLORS, WINDOW)
from succession.timeaxis import match_nearest, to_days                   # noqa: E402

style.apply()
RA_SCALE = 20.0     # right axis 0-100 % against a 0-5 diversity axis


# ── A ────────────────────────────────────────────────────────────────────────
def panel_a():
    d16 = io.load_diversity_16s()
    dbc = io.load_diversity_barcode()
    anc = anchors.table().set_index("mouse")

    fig, axes = plt.subplots(2, 4, figsize=(15, 6.4), sharex=False, sharey=True)
    for ax, m in zip(axes.flat, COLONISED):
        a = d16[(d16["mouse"] == m) & d16["q1"].notna()]
        ax.plot(to_days(a["index"].to_numpy(), group="colonised"), a["q1"],
                "-o", color="#1F77B4", ms=2.6, lw=1.2, label="Community 16S $^1D$")

        # The 3 h barcode sample is near-empty in m1, m2, m4 and m5 -- a handful
        # of reads rather than low diversity -- so it is excluded here.
        b = dbc[(dbc["mouse"] == m) & dbc["q1"].notna()]
        if m in BARCODE_3H_EMPTY:
            b = b[b["index"] != 1]
        ax.plot(to_days(b["index"].to_numpy(), group="colonised"),
                np.log10(b["q1"]), "--^", color="#7B52AB", ms=2.6, lw=1.1,
                label="Barcode $\\log_{10}(^1D)$")

        fam = io.load_family(m)
        for name, colour, marker in (("Paenibacillaceae", "#2CA02C", "o"),
                                     ("Enterobacteriaceae", "#D62728", "s")):
            f = fam[fam["Family"] == name].sort_values("Time")
            ax.plot(to_days(f["Time"].to_numpy(), group="colonised"),
                    f["abundance"] * 100 / RA_SCALE, marker=marker, ms=2.4,
                    lw=1.2, color=colour, label=name)

        ax.axvline(anc.loc[m, "onset_day"], color="#3B8BD0", ls="--", lw=0.9)
        ax.set_ylim(0, 5)
        ax.set_title(m, fontweight="bold")
        # Figure 2A labels the axis in plain days, unlike Figure 1 and 4A which
        # label 3h first and then whole days. Matching the submitted panel.
        ax.set_xlim(0, 16.5 if m in COHORT_1 else 15.5)
        # A tick at every sample; labels only at the whole days the panel names.
        upto = 16 if m in COHORT_1 else 15
        days = [3 / 24, 6 / 24, 12 / 24] + [float(d) for d in range(1, upto + 1)]
        ax.set_xticks(days)
        ax.set_xticklabels([f"{d:g}" if d in (4.0, 8.0, 12.0, 16.0) else ""
                            for d in days])
        sec = ax.secondary_yaxis("right", functions=(lambda y: y * RA_SCALE,
                                                     lambda y: y / RA_SCALE))
        if ax is not axes.flat[3] and ax is not axes.flat[7]:
            sec.set_yticklabels([])
        else:
            sec.set_ylabel("Family relative abundance (%)")

    for ax in axes[:, 0]:
        ax.set_ylabel("Diversity ($^1D$)")
    axes.flat[0].legend(loc="upper left", fontsize=6.5)
    fig.supxlabel("Time post-colonisation (days)", y=0.02)
    fig.suptitle("A   Diversity and the Paenibacillaceae transition", y=0.99)
    fig.tight_layout()
    return style.save(fig, "fig2A_diversity_transition")


# ── B ────────────────────────────────────────────────────────────────────────
# Ported from the manuscript's 12_ridgeline_jacobian.py, which is the script
# behind the submitted panel. Its normalisation is per MOUSE, not per ridge.
BW = 0.20            #: KDE bandwidth, fixed
RIDGE_SCALE = 2.5    #: a ridge's peak reaches ~2.5 real-day units above its base
COL_PRE = "#E07B6A"  #: pre-onset, collapse
COL_POST = "#5B9BD5" #: post-onset, succession
COL_LINE = "#1A1A1A"
MIN_PAIRS = 4        #: fewer than this and a density is not estimated


XLIM = (-0.42, 0.22)   #: the submitted panel's view, from 12_ridgeline_jacobian.py


def panel_b(warmup_tolerance=None, xlim=XLIM):
    """Ridgelines of the off-diagonal distribution, one ridge per evaluation.

    Three choices here decide how the panel reads, and all three follow
    `12_ridgeline_jacobian.py` rather than the per-ridge normalisation used by
    the overview script in the same folder:

    * **Scaled per mouse, not per ridge.** One scale, `RIDGE_SCALE / median of
      the per-ridge peak densities`, is shared by every ridge in a panel. A
      concentrated distribution therefore draws tall and a diffuse one flat,
      which is the signal the panel exists to show. Normalising each ridge to
      its own maximum would make every row the same height and throw that away.
    * **Positioned at the real day.** A ridge sits at its own time on a
      continuous axis, so the sub-day points crowd near the origin and the
      daily samples spread out. This also aligns the eight panels without any
      padding: a mouse with no evaluation at a time simply has no ridge there
      (m6-m8 start at 12 h, m1-m5 at 6 h; cohort 1 runs to 15.9 d, cohort 2 to
      14.9 d).
    * **Every evaluation, not panel C's subset.** The warm-up rule exists so a
      correlation uses a consistent window width; it is not a claim that the
      early windows are wrong, and dropping them would hide the collapse.
      `warmup_tolerance=WARMUP_TOLERANCE` reproduces the panel-C subset.

    `xlim` defaults to the submitted panel's view, (-0.42, 0.22). The full data
    range runs to -1.68, but 99% of the mass sits inside +/-0.5, so drawing it
    uncut compresses every ridge into a spike. Pass `xlim=None` for the full
    range.
    """
    anc = anchors.table().set_index("mouse")

    dists = {}
    for m in COLONISED:
        state = jacobian.build_state(m)
        if warmup_tolerance is None:
            # Slot-count window: the 5 most recent samples, expanding at the
            # start -- the corrected form of what the submitted panel did.
            raw = jacobian.offdiagonal_by_slot(m, state, n_steps=WINDOW)
        else:
            times = jacobian.evaluation_times(m, state, window=WINDOW,
                                              warmup_tolerance=warmup_tolerance)
            raw = jacobian.offdiagonal(state, times, WINDOW)
        dists[m] = {to_days(t, group="colonised"): v
                    for t, v in raw.items() if v.size >= MIN_PAIRS}

    if xlim is None:
        allv = np.concatenate([v for d in dists.values() for v in d.values()])
        pad = 0.05 * (allv.max() - allv.min())
        xlim = (allv.min() - pad, allv.max() + pad)
    xs = np.linspace(xlim[0], xlim[1], 600)

    fig, axes = plt.subplots(2, 4, figsize=(16, 9), sharex=True)

    for ax, m in zip(axes.flat, COLONISED):
        onset = anc.loc[m, "onset_day"]
        densities = {day: gaussian_kde(v, bw_method=BW)(xs)
                     for day, v in dists[m].items()}
        # One scale per mouse, from the median peak, so ridge heights stay
        # comparable across time within this panel.
        scale = RIDGE_SCALE / np.median([d.max() for d in densities.values()])

        for day in sorted(densities):                 # early first, later on top
            dens = densities[day] * scale
            colour = COL_PRE if day <= onset else COL_POST
            ax.fill_between(xs, day, day + dens, color=colour, lw=0, zorder=day)
            ax.plot(xs, day + dens, color=COL_LINE, lw=0.5, zorder=day + 0.1)

        days = sorted(densities)
        ax.axvline(0, color="#888888", ls="--", lw=0.8, zorder=0)
        ax.set_yticks(days)
        ax.set_yticklabels([f"{d:g}" for d in days], fontsize=6.5)
        ax.set_ylim(min(days) - 0.3, max(days) + RIDGE_SCALE + 0.3)
        ax.set_xlim(*xlim)
        ax.set_box_aspect(1)                          # square plot area
        ax.minorticks_on()
        ax.tick_params(axis="x", which="minor", length=2)
        ax.set_title(m, fontweight="bold")

    for ax in axes[:, 0]:
        ax.set_ylabel("Time post-colonisation (days)")
    for ax in axes[1, :]:
        ax.set_xlabel("Pairwise Jacobian elements $J_{ij}$")
    for ax in axes[:, 0]:
        ax.set_ylabel("Time post-colonisation (days)")
    for ax in axes[1, :]:
        ax.set_xlabel("Pairwise Jacobian elements $J_{ij}$")
    handles = [plt.Rectangle((0, 0), 1, 1, fc="#F1948A"),
               plt.Rectangle((0, 0), 1, 1, fc="#7FB3D5")]
    fig.legend(handles, ["pre-onset (collapse)", "post-onset (succession)"],
               loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.02))
    fig.suptitle("B   Distribution of $J_{ij}$ over time", y=0.99)
    fig.tight_layout()
    return style.save(fig, "fig2B_jacobian_ridgelines")


# ── C ────────────────────────────────────────────────────────────────────────
def figure2c_points():
    """The points behind panel C, one row per surviving evaluation."""
    rows = []
    for m in COLONISED:
        state = jacobian.build_state(m)
        times = jacobian.evaluation_times(m, state, window=WINDOW)
        summ = jacobian.summarise(jacobian.offdiagonal(state, times, WINDOW))
        q1 = diversity.hill_q1_from_taxa(m)
        for _, r in summ.iterrows():
            if not np.isfinite(r["mean_negative"]):
                continue
            q = match_nearest(r["index"], q1["index"].to_numpy(), q1["q1"].to_numpy())
            if q is None:
                continue
            rows.append({"mouse": m, "index": r["index"], "q1": q,
                         "mean_negative": r["mean_negative"]})
    import pandas as pd
    return pd.DataFrame(rows)


def panel_c():
    pts = figure2c_points()
    rho, p, n = stats.spearman(pts["q1"], pts["mean_negative"])
    per = {m: stats.spearman(g["q1"], g["mean_negative"])
           for m, g in pts.groupby("mouse")}
    print(f"  pooled rho = {rho:.3f}  P = {p:.2e}  n = {n}  "
          f"significant in {sum(v[1] < 0.05 for v in per.values())} of 8")

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(11, 4.4),
                                   gridspec_kw={"width_ratios": [1, 1.15]})
    for m in COLONISED:
        d = pts[pts["mouse"] == m]
        axl.scatter(d["q1"], d["mean_negative"], s=22, color=MOUSE_COLORS[m],
                    label=m, lw=0)
    b, a = np.polyfit(pts["q1"], pts["mean_negative"], 1)
    xf = np.linspace(pts["q1"].min(), pts["q1"].max(), 50)
    axl.plot(xf, b * xf + a, color="black", lw=1.3)
    axl.text(0.55, 0.28, f"Spearman $\\rho$ = {rho:.2f}\nP = {p:.1e}\nn = {n}",
             transform=axl.transAxes, fontsize=8.5, va="top")
    axl.set_xlabel("16S Hill $^1D$")
    axl.set_ylabel("Mean negative $J_{ij}$")
    axl.set_title("Pooled, all mice")
    axl.legend(ncol=2, fontsize=6.5, loc="lower right")

    vals = [per[m][0] for m in COLONISED]
    axr.bar(range(8), vals, color="#7F7F7F", edgecolor="black", lw=0.5, width=0.62)
    for i, m in enumerate(COLONISED):
        r_, p_, n_ = per[m]
        axr.text(i, r_ + 0.03, f"P = {p_:.1e}\nn = {n_}", ha="center",
                 fontsize=6.2, va="bottom")
    axr.set_xticks(range(8)); axr.set_xticklabels(COLONISED)
    axr.set_ylim(0, 1.25)
    axr.set_xlabel("Mouse")
    axr.set_ylabel("Spearman $\\rho$")
    n_sig = sum(v[1] < 0.05 for v in per.values())
    axr.set_title(f"Per mouse   ({n_sig}/8 significant at window {WINDOW})")
    fig.suptitle("C   Inhibitory interaction versus community diversity", y=1.0)
    fig.tight_layout()
    return style.save(fig, "fig2C_inhibition_vs_diversity")


if __name__ == "__main__":
    print("Figure 2")
    panel_a()
    panel_c()
    panel_b()
