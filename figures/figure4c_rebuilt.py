"""Figure 4C, rebuilt — the colonizer-resident interaction, matched across arms.

Why the published panel could not simply be re-run
--------------------------------------------------
Three faults, only one of which is the clock:

1. CLOCK. Controls are indexed by DAY (1-10), colonised by SAMPLING SLOT
   (1 = 3 h ... 19 = day 16). The submitted comparison window was defined
   relative to the control span, so passing control indices through the
   colonised converter silently narrowed the test to 1-7 d and produced
   P = 0.019. Over the full series it is P ~ 0.27 under either clock.

2. DIMENSION. Control states carry 16S families only (4-7 variables);
   colonised states are dominated by barcode lineages (12 variables). Mean J
   over matrix entries depends on dimension, so the two are not comparable.
   Rebuilding from the unfiltered family tables, the shared subspace across all
   twelve mice is exactly TWO families: Enterobacteriaceae and Paenibacillaceae.
   Controls carry 17-28 families above 0.1%, colonised only 3-6 - the colonizer
   collapses resident diversity, so no wider matched subspace exists.

3. AMPLITUDE. mean |J| says controls are higher (0.060 vs 0.015). Amplitude-free,
   mean |R| says the opposite (0.343 vs 0.522, P = 0.009). The magnitude
   difference is the amplitude term. Only the SIGN survives, so that is what
   this panel reports.

What it shows
-------------
The Enterobacteriaceae-Paenibacillaceae coupling reverses sign between arms, in
both directions, amplitude-free, under the overlap window and the full series.
Under antibiotic alone Enterobacteriaceae suppresses Paenibacillaceae; with the
colonizer present that suppression is gone.

Caveat stated on the figure: "Enterobacteriaceae" is a native population in the
controls and the gavaged K12 in the colonised mice. This compares community
contexts, not the same organism.

Run:  python figures/figure4c_rebuilt.py
"""

from __future__ import annotations

import csv
import collections
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from succession import stats, style                                  # noqa: E402
from succession.config import (COLONISED, CONTROLS, DATA,            # noqa: E402
                               PSEUDOCOUNT)
from succession.timeaxis import to_days                              # noqa: E402

style.apply()

PAIR = ["Enterobacteriaceae", "Paenibacillaceae"]
WIN_DAYS = 5.0          #: identical real duration in both arms — not slots
OVERLAP = (1.0, 10.0)   #: the days both arms actually cover
CTL_C, COL_C = "#B8860B", "#2471A3"


def _family_series(mouse: str):
    by = collections.defaultdict(lambda: collections.defaultdict(float))
    seen = set()
    path = DATA / "16s" / "family" / f"{mouse}_family.csv"
    with open(path) as fh:
        for row in csv.DictReader(fh):
            t = float(row["Time"])
            seen.add(t)
            by[row["Family"]][t] += float(row["Abundance.family"])
    idx = np.array(sorted(seen))
    group = "control" if mouse.startswith("c_") else "colonised"
    days = to_days(idx, group=group)
    ok = ~np.isnan(days)            # m5 carries a pre-gavage index 0
    return idx[ok], days[ok], by


def state(mouse: str, step: float = 0.1):
    idx, days, by = _family_series(mouse)
    grid = np.arange(days[0], days[-1] + step / 2, step)
    cols = []
    for fam in PAIR:
        y = np.log10(np.array([by[fam][t] for t in idx]) + PSEUDOCOUNT)
        v = CubicSpline(days, y, extrapolate=False)(grid)
        bad = np.isnan(v)
        if bad.any():
            good = np.flatnonzero(~bad)
            v[bad] = np.interp(np.flatnonzero(bad), good, v[good])
        cols.append(v)
    return grid, np.column_stack(cols), days, by, idx


def interactions(mouse: str):
    """Directed J and R for the two shared families, window fixed in DAYS."""
    grid, Z, days, _, _ = state(mouse)
    step = grid[1] - grid[0]
    dZ = np.gradient(Z, step, axis=0)
    rows = []
    for te in days:
        if te - WIN_DAYS < grid[0] - 2.0:
            continue
        hi = int(np.searchsorted(grid, te, side="right")) - 1
        lo = max(0, int(np.searchsorted(grid, te - WIN_DAYS, side="left")))
        if hi - lo < 4:
            continue
        zw, dw = Z[lo:hi + 1], dZ[lo:hi + 1]
        J = np.array([[np.cov(dw[:, i], zw[:, j])[0, 1] for j in range(2)]
                      for i in range(2)])
        sz, sd = zw.std(axis=0, ddof=1), dw.std(axis=0, ddof=1)
        R = np.full((2, 2), np.nan)
        for i in range(2):
            for j in range(2):
                if sd[i] > 1e-12 and sz[j] > 1e-12:
                    R[i, j] = J[i, j] / (sd[i] * sz[j])
        if np.isfinite(J).all():
            rows.append((float(te), J[1, 0], J[0, 1], R[1, 0], R[0, 1]))
    return np.array(rows) if rows else np.empty((0, 5))


def _pool(data, mice, col, lo, hi):
    out = []
    for m in mice:
        a = data[m]
        if len(a):
            sel = a[(a[:, 0] >= lo) & (a[:, 0] <= hi), col]
            out += [v for v in sel if np.isfinite(v)]
    return np.array(out)


def main():
    data = {m: interactions(m) for m in COLONISED + CONTROLS}
    abund = {}
    for m in COLONISED + CONTROLS:
        _, days, by = _family_series(m)[0], _family_series(m)[1], _family_series(m)[2]
        abund[m] = (days, np.array([by["Paenibacillaceae"][t]
                                    for t in _family_series(m)[0]]))

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.5))

    # A — Paenibacillaceae abundance, measured rather than assumed zero
    ax = axes[0]
    for m in COLONISED:
        d, v = abund[m]
        ax.plot(d, 100 * v, "-", lw=1.1, color=COL_C, alpha=.75)
    for m in CONTROLS:
        d, v = abund[m]
        ax.plot(d, 100 * v, "-o", ms=3.4, lw=1.3, color=CTL_C)
    ax.set_yscale("symlog", linthresh=1)
    ax.set_xlabel("day"); ax.set_ylabel("Paenibacillaceae (%)")
    ax.set_title("A   Present in every control, never blooming", fontsize=10)
    ax.text(.03, .95, "colonised  53–67%\ncontrols  max 0.36–1.4%,\n4/4 mice",
            transform=ax.transAxes, va="top", fontsize=8.5,
            bbox=dict(fc="white", ec="0.75", alpha=.9, boxstyle="round,pad=0.3"))

    # B — the sign reversal, amplitude-free
    ax = axes[1]
    pairs, labels, ps = [], [], []
    for col, lab in ((3, "R[Paeni<-Entero]"), (4, "R[Entero<-Paeni]")):
        a = _pool(data, CONTROLS, col, 1, 16)
        b = _pool(data, COLONISED, col, 1, 16)
        _, p = stats.mann_whitney(a, b, alternative="two-sided")
        pairs += [a, b]; labels += [f"ctl\nn={len(a)}", f"col\nn={len(b)}"]; ps.append(p)
    bp = ax.boxplot(pairs, widths=.55, showfliers=False, tick_labels=labels,
                    patch_artist=True)
    for i, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(CTL_C if i % 2 == 0 else COL_C); patch.set_alpha(.55)
    for i, v in enumerate(pairs, start=1):
        ax.scatter(np.random.default_rng(0).normal(i, .06, len(v)), v, s=7,
                   color="0.25", alpha=.5, lw=0)
    ax.axhline(0, color="k", lw=.8, ls="--")
    ax.set_ylabel("scale-free interaction  R")
    ax.set_title("B   The coupling reverses sign", fontsize=10)
    for k, p in enumerate(ps):
        ax.text(1.5 + 2 * k, ax.get_ylim()[1] * .92, f"P = {p:.1e}",
                ha="center", fontsize=8.5)
    for k, lab in enumerate(("R[Paeni <- Entero]", "R[Entero <- Paeni]")):
        ax.text(1.5 + 2 * k, -.14, lab, transform=ax.get_xaxis_transform(),
                ha="center", va="top", fontsize=9, fontweight="bold")

    # C — why magnitude is not the story
    ax = axes[2]
    mj_c = np.abs(np.concatenate([_pool(data, CONTROLS, c, 1, 16) for c in (1, 2)]))
    mj_l = np.abs(np.concatenate([_pool(data, COLONISED, c, 1, 16) for c in (1, 2)]))
    mr_c = np.abs(np.concatenate([_pool(data, CONTROLS, c, 1, 16) for c in (3, 4)]))
    mr_l = np.abs(np.concatenate([_pool(data, COLONISED, c, 1, 16) for c in (3, 4)]))
    _, p1 = stats.mann_whitney(mj_c, mj_l, alternative="two-sided")
    _, p2 = stats.mann_whitney(mr_c, mr_l, alternative="two-sided")
    # each arm expressed relative to its own control median, so |J| and |R| are comparable
    ratios = [1.0, np.median(mj_l) / np.median(mj_c), 1.0, np.median(mr_l) / np.median(mr_c)]
    xs = [1, 2, 3.6, 4.6]
    for x, v, c in zip(xs, ratios, [CTL_C, COL_C, CTL_C, COL_C]):
        ax.bar(x, v, width=.82, color=c, alpha=.78)
    ax.axhline(1, color="0.4", lw=.9, ls="--")
    ax.set_xticks([1.5, 4.1])
    ax.set_xticklabels(["$|J|$  raw", "$|R|$  amplitude-free"])
    ax.set_ylabel("median, relative to control")
    ax.set_ylim(0, 1.85)
    ax.set_title("C   Magnitude reverses — so it is amplitude", fontsize=10)
    ax.text(1.5, ratios[1] + 0.12, f"colonised {ratios[1]:.2f}×\nP = {p1:.1e}",
            ha="center", fontsize=8.5)
    ax.text(4.1, ratios[3] + 0.12, f"colonised {ratios[3]:.2f}×\nP = {p2:.1e}",
            ha="center", fontsize=8.5)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=CTL_C, alpha=.78, label="control"),
                       plt.Rectangle((0, 0), 1, 1, fc=COL_C, alpha=.78, label="colonised")],
              fontsize=8, loc="upper left", framealpha=.9)

    fig.text(.5, -.03,
             "Two shared families, window fixed at 5 days in both arms, control indices on the control clock.  "
             "“Enterobacteriaceae” is a native population in controls and the gavaged K12 in colonised mice.",
             ha="center", fontsize=8, color="#555")
    fig.tight_layout()
    out = style.save(fig, "fig4C_rebuilt")
    print("Rebuilt Figure 4C")
    for lab, col in (("R[Paeni<-Entero]", 3), ("R[Entero<-Paeni]", 4)):
        a, b = _pool(data, CONTROLS, col, 1, 16), _pool(data, COLONISED, col, 1, 16)
        _, p = stats.mann_whitney(a, b, alternative="two-sided")
        print(f"  {lab:20} ctl med={np.median(a):+.4f}  col med={np.median(b):+.4f}  P={p:.3g}")
    return out


if __name__ == "__main__":
    main()
