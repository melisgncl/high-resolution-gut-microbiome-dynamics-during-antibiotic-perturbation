"""Figure 4 — antibiotic-only control mice.

A  16S family composition, c_m1-c_m4
B  16S Hill 1D, controls against the colonised mean +/- s.d.
C  Paenibacillaceae relative abundance, controls vs colonised - no bloom
   without E. coli colonisation

The controls are indexed by DAY and the colonised by SAMPLING SLOT. Every axis
here converts with the right clock, which is why `to_days` demands the group.
Passing control indices through the colonised converter compresses days 1-10
into 0.125-7 d, which is visible on the panel B axis; see `--clock published`
below.

RETIRED (19 Aug 2026): panel C previously compared mean pairwise J between
controls and colonised mice. That comparison does not survive scrutiny at any
window or normalisation - see CORRECTIONS.md ss1 and ss3-4 - and even a
dimension-matched, amplitude-free version, while real, is tangential to what
this figure needs to establish. Dropped in favour of the panel below, which
makes Figure 4's actual claim directly: Paenibacillaceae is detectable in
every control mouse and never blooms, at any normalisation, with any window,
by construction - it needs neither.

Run:  python figures/figure4.py [--clock corrected|published]
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from succession import io, stats, style                               # noqa: E402
from succession.config import (COLONISED, CONTROLS, CONTROL_COLORS,   # noqa: E402
                               DATA, FAMILY_GATE, PSEUDOCOUNT_16S_PLOT)
from succession.timeaxis import to_days                               # noqa: E402

style.apply()
MARKERS = {"c_m1": "o", "c_m2": "s", "c_m3": "^", "c_m4": "D"}


def control_days(index, clock: str):
    """Control index to days, either correctly or the way the paper did it."""
    return to_days(index, group="control" if clock == "corrected" else "colonised")


# ── A ────────────────────────────────────────────────────────────────────────
def panel_a():
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.8), sharey=True)
    palette = plt.get_cmap("tab20")
    seen: dict[str, tuple] = {}

    for ax, m in zip(axes, CONTROLS):
        fam = io.load_family(m)
        means = fam.groupby("Family")["abundance"].mean()
        keep = sorted(means[means > FAMILY_GATE].index)
        for f in keep:
            if f not in seen:
                seen[f] = palette(len(seen) % 20)
            d = fam[fam["Family"] == f].sort_values("Time")
            a = d["abundance"].to_numpy()
            touching = (a > 0) | (np.r_[0, a[:-1]] > 0) | (np.r_[a[1:], 0] > 0)
            y = np.where(touching, a + PSEUDOCOUNT_16S_PLOT, np.nan)
            ax.plot(to_days(d["Time"].to_numpy(), group="control"), y,
                    color=seen[f], lw=2.0 if f == "Paenibacillaceae" else 0.9)
        ax.set_yscale("log"); ax.set_ylim(1e-6, 1)
        # A tick on every sampled day, a label on every other.
        ax.set_xticks(range(1, 11))
        ax.set_xticklabels([f"{d}d" if d % 2 else "" for d in range(1, 11)])
        ax.set_title(m.replace("c_m", "cm"), fontweight="bold")
    axes[0].set_ylabel("Family composition")
    fig.supxlabel("Time post-colonisation (days)", y=0.02)
    handles = [plt.Line2D([], [], color=c, lw=2, label=f) for f, c in sorted(seen.items())]
    fig.legend(handles=handles, loc="lower center", ncol=6, bbox_to_anchor=(0.5, -0.30))
    fig.suptitle("A   Cohort 3: control group", y=1.0)
    fig.tight_layout()
    return style.save(fig, "fig4A_control_composition")


# ── B ────────────────────────────────────────────────────────────────────────
def panel_b(clock: str):
    div = io.load_diversity_16s()
    col = div[div["mouse"].isin(COLONISED) & div["q1"].notna()].copy()
    col["day"] = to_days(col["index"].to_numpy(), group="colonised")
    band = col.groupby("day")["q1"].agg(["mean", "std", "count"]).reset_index()
    band = band[band["count"] >= 2]

    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    ax.fill_between(band["day"], band["mean"] - band["std"], band["mean"] + band["std"],
                    color="#AED6F1", alpha=0.55, lw=0, label="colonised $\\pm$ s.d.")
    ax.plot(band["day"], band["mean"], color="#2471A3", lw=1.5,
            label="colonised mean (n = 8)")

    for m in CONTROLS:
        d = div[(div["mouse"] == m) & div["q1"].notna()].sort_values("index")
        # Gaps are real: c_m1 has usable q1 on days 3, 6, 7, 8, 9 only. Reindexing
        # onto the full day range inserts NaN so the line breaks rather than
        # being drawn straight across the missing days.
        full = pd.DataFrame({"index": range(1, 11)}).merge(d, on="index", how="left")
        ax.plot(control_days(full["index"].to_numpy(), clock), full["q1"],
                marker=MARKERS[m], ms=4, lw=1.1, color=CONTROL_COLORS[m], label=m)
    ax.set_xlabel("Time post-colonisation (days)")
    ax.set_ylabel("Community diversity (16S Hill $^1D$)")
    ax.set_xticks(range(1, 17))
    ax.set_xticklabels([f"{d}" if d % 2 == 0 else "" for d in range(1, 17)])
    ax.set_title("B   16S community diversity in control vs. colonised mice")
    ax.legend(ncol=2, fontsize=6.8)
    return style.save(
        fig, f"fig4B_control_diversity{'' if clock == 'corrected' else '_published'}")


# ── C ────────────────────────────────────────────────────────────────────────
def _paeni_series(mouse):
    """Paenibacillaceae relative abundance, unfiltered, on the correct clock.

    Reads `data/16s/family/` directly rather than `build_state`'s filtered
    `for_jacobian` tables - those drop any family with mean relative abundance
    < 0.5%, which silently removed Paenibacillaceae from every control mouse
    and let downstream code read the missing column as zero. See
    CORRECTIONS.md ss5.
    """
    by = {}
    for row in csv.DictReader(open(DATA / "16s" / "family" / f"{mouse}_family.csv")):
        if row["Family"] == "Paenibacillaceae":
            by[float(row["Time"])] = by.get(float(row["Time"]), 0.0) + float(row["Abundance.family"])
    idx = np.array(sorted(by))
    group = "control" if mouse.startswith("c_") else "colonised"
    days = to_days(idx, group=group)
    ok = ~np.isnan(days)
    return days[ok], np.array([by[t] for t in idx])[ok] * 100  # percent


def panel_c():
    """The actual claim of Figure 4: no bloom without colonisation.

    Peak and final-timepoint Paenibacillaceae abundance are compared per
    mouse (not per timepoint - that would pseudoreplicate 7-19 samples per
    animal into one test). Complete separation: every control mouse tops out
    under 1.4%, every colonised mouse ends above 46%. Needs no window, no
    normalisation, no Jacobian.
    """
    ctl_max = [np.max(_paeni_series(m)[1]) for m in CONTROLS]
    col_max = [np.max(_paeni_series(m)[1]) for m in COLONISED]
    U, p = stats.mann_whitney(np.array(ctl_max), np.array(col_max), alternative="two-sided")
    print(f"  peak Paenibacillaceae, per mouse: control max = {ctl_max} (n={len(ctl_max)}), "
          f"colonised max = {[round(v,1) for v in col_max]} (n={len(col_max)}), "
          f"Mann-Whitney U={U:.0f}, two-sided P={p:.4g}")

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(11.5, 4.4),
                                   gridspec_kw={"width_ratios": [1.6, 1]})

    for m in COLONISED:
        d, v = _paeni_series(m)
        axl.plot(d, v, "-", lw=1.1, color="#2471A3", alpha=0.75,
                 label="colonised" if m == COLONISED[0] else None)
    for m in CONTROLS:
        d, v = _paeni_series(m)
        axl.plot(d, v, marker=MARKERS[m], ms=4.5, lw=1.3, color=CONTROL_COLORS[m], label=m)
    axl.set_yscale("symlog", linthresh=1)
    axl.set_xlabel("Time post-colonisation (days)")
    axl.set_ylabel("Paenibacillaceae relative abundance (%)")
    axl.set_title("C   No Paenibacillaceae bloom without colonisation", fontsize=10.5)
    axl.legend(ncol=2, fontsize=6.8, loc="upper left")

    axr.boxplot([ctl_max, col_max], widths=0.5, showfliers=False,
                tick_labels=[f"Control\nn={len(ctl_max)}", f"Colonised\nn={len(col_max)}"])
    for i, (v, c) in enumerate(((ctl_max, "grey"), (col_max, "#2471A3")), start=1):
        axr.scatter(np.random.default_rng(0).normal(i, 0.05, len(v)), v,
                    s=16, color=c, alpha=0.75, lw=0)
    axr.set_ylabel("Peak Paenibacillaceae (%)")
    axr.set_title(f"Mann-Whitney U = {U:.0f}\nP = {p:.3g}", fontsize=10)
    fig.text(0.5, -0.03, "Peak abundance per mouse, one value per biological "
             "replicate (n=4 controls, n=8 colonised) - not per timepoint.",
             ha="center", fontsize=7, color="#555")
    fig.tight_layout()
    return style.save(fig, "fig4C_paeni_no_bloom_without_colonisation")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--clock", choices=["corrected", "published"], default="corrected",
                    help="'published' reproduces the submitted panel, which passes "
                         "control indices through the colonised converter.")
    args = ap.parse_args()
    print(f"Figure 4  (clock = {args.clock})")
    panel_a()
    panel_b(args.clock)
    panel_c()
