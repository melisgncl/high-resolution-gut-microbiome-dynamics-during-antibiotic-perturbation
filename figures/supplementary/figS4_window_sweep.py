"""Supplementary Figure 4 — the diversity/inhibition correlation across windows.

A  pooled Spearman rho against window width 1-10, with the sample size each
   width leaves. The published caption quotes n = 73-129, which is what the
   evaluation rule produces: 129 at widths 1-3 down to 73 at width 10, losing
   exactly one evaluation per mouse per extra step.
B  per-mouse scatter at the primary width, coloured by day.

Run:  python figures/supplementary/figS4_window_sweep.py
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from _common import (COLONISED, diversity, jacobian, match_nearest, stats,
                     style, to_days, WINDOW, WINDOW_ROBUSTNESS)

style.apply()


def points_at(window: int) -> pd.DataFrame:
    rows = []
    for m in COLONISED:
        state = jacobian.build_state(m)
        times = jacobian.evaluation_times(m, state, window=window)
        summ = jacobian.summarise(jacobian.offdiagonal(state, times, window))
        q1 = diversity.hill_q1_from_taxa(m)
        for _, r in summ.iterrows():
            if not np.isfinite(r["mean_negative"]):
                continue
            q = match_nearest(r["index"], q1["index"].to_numpy(), q1["q1"].to_numpy())
            if q is None:
                continue
            rows.append({"mouse": m, "index": r["index"], "q1": q,
                         "mean_negative": r["mean_negative"],
                         "day": to_days(r["index"], group="colonised")})
    return pd.DataFrame(rows)


def main():
    sweep = []
    for w in range(1, 11):
        p = points_at(w)
        rho, pv, n = stats.spearman(p["q1"], p["mean_negative"])
        nsig = sum(stats.spearman(g["q1"], g["mean_negative"])[1] < 0.05
                   for _, g in p.groupby("mouse"))
        sweep.append({"window": w, "rho": rho, "p": pv, "n": n, "n_sig": nsig})
        print(f"  window {w:2d}: rho = {rho:+.3f}  n = {n:3d}  {nsig}/8 significant")
    sweep = pd.DataFrame(sweep)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 4.4),
                                  gridspec_kw={"width_ratios": [1, 1.25]})
    ax.plot(sweep["window"], sweep["rho"], "o-", color="#2C3E50", lw=1.6)
    for _, r in sweep.iterrows():
        ax.annotate(f"n={int(r['n'])}", (r["window"], r["rho"]),
                    textcoords="offset points", xytext=(0, 8), ha="center", fontsize=6)
    ax.axvline(WINDOW, color="#C0392B", ls="-", lw=1.2, zorder=0,
               label=f"primary window ({WINDOW}) — Fig. 2C, Fig. 3, Fig. 4")
    ax.axvline(WINDOW_ROBUSTNESS, color="#7F7F7F", ls="--", lw=0.9, zorder=0,
               label=f"robustness check ({WINDOW_ROBUSTNESS}) — this figure only")
    for w in sweep["window"]:
        nsig = int(sweep.loc[sweep["window"] == w, "n_sig"].iloc[0])
        ax.annotate(f"{nsig}/8", (w, sweep.loc[sweep["window"] == w, "rho"].iloc[0]),
                    textcoords="offset points", xytext=(0, -16), ha="center",
                    fontsize=6.5, color="#333", zorder=5,
                    bbox=dict(fc="white", ec="none", pad=0.5, alpha=.95))
    ax.set_xticks(range(1, 11))
    ax.set_xlabel("Sliding-window width (sampling slots)")
    ax.set_ylabel(r"Pooled Spearman $\rho$")
    ax.set_title("A   Correlation across window widths, with per-mouse significance")
    ax.legend(fontsize=6.5)

    p = points_at(WINDOW)
    sc = None
    for m in COLONISED:
        d = p[p["mouse"] == m]
        sc = ax2.scatter(d["q1"], d["mean_negative"], c=d["day"], cmap="plasma",
                         s=24, lw=0, vmin=p["day"].min(), vmax=p["day"].max())
    cb = fig.colorbar(sc, ax=ax2, pad=0.015)
    cb.set_label("Day post-colonisation")
    ax2.set_xlabel("16S Hill $^1D$")
    ax2.set_ylabel("Mean negative $J_{ij}$")
    ax2.set_title(f"B   All mice pooled, primary window = {WINDOW}")
    fig.tight_layout()
    return style.save(fig, "figS4_window_sweep")


if __name__ == "__main__":
    print("Supplementary Figure 4")
    main()
