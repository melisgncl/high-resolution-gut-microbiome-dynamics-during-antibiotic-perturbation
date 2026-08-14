"""Figure 4 — antibiotic-only control mice.

A  16S family composition, c_m1-c_m4
B  16S Hill 1D, controls against the colonised mean +/- s.d.
C  mean pairwise J over time, and the control/colonised comparison

The controls are indexed by DAY and the colonised by SAMPLING SLOT. Every axis
here converts with the right clock, which is why `to_days` demands the group.
Passing control indices through the colonised converter compresses days 1-10
into 0.125-7 d, which is visible on the panel B and C axes; see `--clock
published` below.

The panel C *test* is decided by the comparison window, not by the clock. With
`COLONISED_WINDOW` fixed at 1-7 d both clocks give the same one-sided
P = 0.019; widening it to the full series gives P = 0.27 under either clock.
The clock mattered in the submitted figure only because the window was defined
relative to the control span, so a mislabelled axis quietly resized the test.

Controls carry no barcode dimension, so their Jacobian is estimated over 16S
taxa only. The two groups' J values are therefore not built from the same number
of variables, which is stated on the panel.

Run:  python figures/figure4.py [--clock corrected|published]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from succession import io, jacobian, stats, style                    # noqa: E402
from succession.config import (COLONISED, CONTROLS, CONTROL_COLORS,   # noqa: E402
                               FAMILY_GATE, PSEUDOCOUNT_16S_PLOT, WINDOW)
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
def _mean_j(mouse):
    state = jacobian.build_state(mouse)
    times = jacobian.evaluation_times(mouse, state, window=WINDOW)
    summ = jacobian.summarise(jacobian.offdiagonal(state, times, WINDOW))
    return summ[["index", "mean"]]


COLONISED_WINDOW = (1.0, 7.0)
"""Days of the colonised series entering the panel-C test.

Fixed at 1-7 d rather than "every colonised evaluation up to the last control
day". The controls only have evaluations from 4 to 9.9 d, so an open-ended rule
lets the comparison set drift with the clock: 77 colonised points under the
correct clock, 53 under the published one. Naming the window makes the test the
same test in both cases, and confines it to the collapse-and-onset period the
controls actually cover.
"""


def comparison(clock: str, colonised_window=COLONISED_WINDOW) -> dict:
    """The control-vs-colonised test behind panel C, under either clock.

    Separated from the drawing so the notebook can state the numbers, and so
    the published/corrected difference can be shown as a table rather than
    only as two pictures. The clock decides which colonised evaluations are
    in the comparison set: under the colonised converter the controls appear
    to span 0.125-7 d instead of 1-10 d, restricting the comparison to the
    collapse phase where J is most negative.
    """
    col = pd.concat([_mean_j(m).assign(mouse=m) for m in COLONISED])
    col["day"] = to_days(col["index"].to_numpy(), group="colonised")
    ctl = pd.concat([_mean_j(m).assign(mouse=m) for m in CONTROLS])
    ctl["day"] = control_days(ctl["index"].to_numpy(), clock)

    max_ctl = float(ctl["day"].max())
    lo, hi = colonised_window
    a = ctl["mean"].to_numpy()
    b = col.loc[col["day"].between(lo, hi), "mean"].to_numpy()
    U, p = stats.mann_whitney(a, b, alternative="greater")
    return {"clock": clock, "col": col, "ctl": ctl, "a": a, "b": b,
            "U": U, "p": p, "max_ctl": max_ctl,
            "min_ctl": float(ctl["day"].min()), "window": (lo, hi)}


def panel_c(clock: str):
    cmp = comparison(clock)
    col, ctl = cmp["col"], cmp["ctl"]
    a, b, p, max_ctl = cmp["a"], cmp["b"], cmp["p"], cmp["max_ctl"]
    print(f"  clock={clock}: controls span {cmp['min_ctl']:g}-{max_ctl:g} d | "
          f"colonised window {cmp['window'][0]:g}-{cmp['window'][1]:g} d | "
          f"n control = {len(a)}, n colonised = {len(b)} | one-sided P = {p:.4g}")

    band = col.groupby("day")["mean"].agg(["mean", "std", "count"]).reset_index()
    band = band[band["count"] >= 2]

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(11.5, 4.4),
                                   gridspec_kw={"width_ratios": [2, 1]})
    axl.axvspan(ctl["day"].min(), max_ctl, color="grey", alpha=0.12, lw=0)
    axl.axhline(0, color="grey", ls="--", lw=0.5)
    axl.fill_between(band["day"], band["mean"] - band["std"], band["mean"] + band["std"],
                     color="#AED6F1", alpha=0.55, lw=0)
    axl.plot(band["day"], band["mean"], color="#2471A3", lw=1.5, label="colonised")
    for m in CONTROLS:
        d = ctl[ctl["mouse"] == m].sort_values("day")
        axl.plot(d["day"], d["mean"], marker=MARKERS[m], ms=4, lw=1.1,
                 color=CONTROL_COLORS[m], label=m)
    axl.set_xlabel("Time post-colonisation (days)")
    axl.set_ylabel("Mean $J_{ij}$")
    axl.set_title("C   Community interaction in control vs. colonised mice")
    axl.legend(ncol=2, fontsize=6.8)

    axr.boxplot([a, b], widths=0.5, showfliers=False,
                tick_labels=[f"Control\nn={len(a)}",
                             f"Colonised\n(days {cmp['window'][0]:g}-"
                             f"{cmp['window'][1]:g})\nn={len(b)}"])
    for i, (v, c) in enumerate(((a, "grey"), (b, "#2471A3")), start=1):
        axr.scatter(np.random.default_rng(0).normal(i, 0.05, len(v)), v,
                    s=9, color=c, alpha=0.65, lw=0)
    axr.axhline(0, color="grey", ls="--", lw=0.5)
    axr.set_ylabel("Mean $J_{ij}$")
    axr.set_title(f"Mann-Whitney\nP = {p:.3f}")
    fig.text(0.5, -0.04, "Control Jacobians are 16S-only: the two groups are not "
             "built from the same number of variables.", ha="center", fontsize=7,
             color="#555")
    fig.tight_layout()
    return style.save(
        fig, f"fig4C_control_jacobian{'' if clock == 'corrected' else '_published'}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--clock", choices=["corrected", "published"], default="corrected",
                    help="'published' reproduces the submitted panel, which passes "
                         "control indices through the colonised converter.")
    args = ap.parse_args()
    print(f"Figure 4  (clock = {args.clock})")
    panel_a()
    panel_b(args.clock)
    panel_c(args.clock)
