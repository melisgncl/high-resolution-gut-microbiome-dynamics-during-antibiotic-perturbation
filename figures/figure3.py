"""Figure 3 — community stability and directed colonizer/resident interactions.

A  Jacobian eigenvalue spectrum for one mouse, coloured by time
B  Re(lambda) against time, every eigenvalue of every mouse
C  Resident -> Colonizer:  J[clone <- Paenibacillaceae]
D  Colonizer -> Resident:  J[Paenibacillaceae <- clone]

Direction convention: J[i <- j] = cov(dz_i/dt, z_j); i responds, j acts.

Run:  python figures/figure3.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from succession import anchors, io, jacobian, stats, style          # noqa: E402
from succession.config import (COLONISED, MOUSE_COLORS_TAB,         # noqa: E402
                               WINDOW_PRIMARY as WINDOW)
from succession.timeaxis import day_tick_marks, to_days             # noqa: E402

style.apply()
FOCAL = "m1"


def all_eigenvalues() -> pd.DataFrame:
    """Full spectrum, one mouse's-worth at a time — panel A's complex plane.

    Sliding-window (``all_eigenvalues_sliding``), not the stored expanding-window
    matrices: rank no longer tracks time. See CORRECTIONS.md.
    """
    eig = pd.concat([jacobian.all_eigenvalues_sliding(m, WINDOW) for m in COLONISED],
                    ignore_index=True)
    eig["day"] = to_days(eig["index"].to_numpy(), group="colonised")
    return eig


def dominant_eigenvalues() -> pd.DataFrame:
    """One Re(lambda_max) per window per mouse — panel B's trend.

    Not the same table as ``all_eigenvalues``: that one is pseudoreplicated
    (windows x species) and is for the panel-A visual only. This is the
    corrected counterpart of the published Figure 3B statistic. See
    ``jacobian.dominant_eigenvalue`` and ``CORRECTIONS.md``.
    """
    dom = pd.concat([jacobian.dominant_eigenvalue(m, WINDOW) for m in COLONISED],
                    ignore_index=True)
    dom["day"] = to_days(dom["index"].to_numpy(), group="colonised")
    dom = dom.rename(columns={"re_max": "re", "im_max": "im"})
    return dom


# ── A ────────────────────────────────────────────────────────────────────────
def panel_a(eig):
    """Complex plane for one mouse.

    x is Re(lambda) and y is Im(lambda), which is what the published code
    computes; the submitted panel's axis labels are transposed. Stability is
    Re < 0, so the boundary is the vertical line and the shaded half-plane.
    """
    e = eig[eig["mouse"] == FOCAL]
    xr = np.percentile(eig["re"], [1, 99]) * 1.15
    yr = np.abs(np.percentile(eig["im"], [1, 99])).max() * 1.1

    fig, ax = plt.subplots(figsize=(5.4, 4.6))
    ax.axvspan(xr[0], 0, color="#e8f4e8", alpha=0.45, zorder=0)
    ax.axhline(0, color="#aaaaaa", lw=0.5, zorder=1)
    ax.axvline(0, color="#cc2222", lw=0.9, ls="--", zorder=3)
    # RdYlBu_r: early = blue, mid = pale yellow, late = red. This is the
    # colourway of the submitted panel (33b_figure_eigenvalues_zoom.py); the
    # plasma map in 33_figure_eigenvalues.py is a different variant.
    sc = ax.scatter(e["re"], e["im"], c=e["relative_time"], cmap="RdYlBu_r",
                    s=20, alpha=0.8, lw=0, zorder=2)
    cb = fig.colorbar(sc, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label("Time"); cb.set_ticks([0, 0.5, 1])
    cb.set_ticklabels(["Early", "Mid", "Late"])
    ax.text(0.03, 0.96, "Stable\n(Re < 0)", transform=ax.transAxes,
            fontsize=7, color="#336633", va="top")
    ax.text(0.60, 0.96, "Unstable\n(Re > 0)", transform=ax.transAxes,
            fontsize=7, color="#cc2222", va="top")
    ax.set_xlim(*xr); ax.set_ylim(-yr, yr)
    ax.set_xlabel(r"Re($\lambda$)"); ax.set_ylabel(r"Im($\lambda$)")
    ax.set_title(f"A   mouse {FOCAL}")
    return style.save(fig, "fig3A_eigenvalue_spectrum")


# ── B ────────────────────────────────────────────────────────────────────────
def panel_b(eig):
    """``eig`` is ``dominant_eigenvalues()`` — one Re(lambda_max) per window per
    mouse, not the full spectrum. Plotting the full spectrum here would
    pseudoreplicate the trend statistic; see the docstring on that function."""
    rho, p, n = stats.spearman(eig["day"], eig["re"])
    print(f"  Re(lambda_max) vs time [dominant, sliding window {WINDOW}]: "
          f"rho = {rho:.3f}  P = {p:.2e}  n = {n}   "
          f"(published: rho = -0.43, n = 1545, confounded + pseudoreplicated)")
    mean_onset = anchors.table()["onset_day"].mean()

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.axhspan(-10, 0, color="#e8f4e8", alpha=0.45, zorder=0)
    ax.axhline(0, color="#cc2222", ls="--", lw=0.8, zorder=2)
    ax.axvline(mean_onset, color="#2E86C1", ls="-.", lw=1.0, zorder=2,
               label="Paenibacillaceae onset (mean)")
    for m in COLONISED:
        d = eig[eig["mouse"] == m]
        ax.scatter(d["day"], d["re"], s=6, color=MOUSE_COLORS_TAB[m], label=m,
                   alpha=0.85, lw=0, zorder=3)
    # The dominant eigenvalue only, not the full spectrum, so the range is
    # narrower than the old panel and does not need clipping.
    lo, hi = float(eig["re"].min()), float(eig["re"].max())
    pad = 0.06 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)
    ax.text(0.62, 0.93, f"Spearman $\\rho$ = {rho:.2f}\nP = {p:.2e}\nn = {n}",
            transform=ax.transAxes, fontsize=8, va="top")
    frac_pos = float((eig["re"] > 0).mean())
    ax.text(0.62, 0.72, f"Re($\\lambda_{{max}}$) > 0\nat {100*frac_pos:.0f}% of windows",
            transform=ax.transAxes, fontsize=7.5, va="top", color="#cc2222")
    ax.text(0.01, 0.97, "Unstable", transform=ax.transAxes, color="#cc2222",
            fontsize=7.5, va="top")
    ax.text(0.01, 0.04, "Stable", transform=ax.transAxes, color="#336633",
            fontsize=7.5)
    ax.set_xlabel("Time (days post-colonisation)")
    ax.set_ylabel(r"Re($\lambda_{max}$)")
    ax.set_title(r"B   Time vs Re($\lambda_{max}$) — dominant eigenvalue, one per window")
    # A tick at every sampling day, labelled every other one.
    day_tick_marks(ax, eig["day"].to_numpy(), label_every=2)
    ax.legend(ncol=1, fontsize=7, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              borderaxespad=0)
    ax.set_box_aspect(1)              # square plot area, for figure assembly
    return style.save(fig, "fig3B_eigenvalues_vs_time")


# ── C / D ────────────────────────────────────────────────────────────────────
def panels_cd():
    """Directed coefficients between the colonizer's clones and the resident.

    Within a mouse the value is the mean over its clonal clusters; across mice
    it is mean +/- s.d. Enterobacteriaceae is the single 16S taxon standing for
    E. coli in the community model, drawn dashed.
    """
    recs = []
    for m in COLONISED:
        state = jacobian.build_state(m)
        times = jacobian.evaluation_times(m, state, window=WINDOW)
        clones = state.clone_indices()
        paeni = state.index_of("Paenibacillaceae")
        entero = state.index_of("Enterobacteriaceae")
        if paeni is None or not clones:
            continue
        series = {
            ("C", "clones"): jacobian.directed_group(state, times, clones, [paeni], WINDOW),
            ("D", "clones"): jacobian.directed_group(state, times, [paeni], clones, WINDOW),
        }
        if entero is not None:
            series[("C", "entero")] = jacobian.directed(state, times, entero, paeni, WINDOW)
            series[("D", "entero")] = jacobian.directed(state, times, paeni, entero, WINDOW)
        for (panel, kind), d in series.items():
            for t, v in d.items():
                recs.append({"panel": panel, "kind": kind, "mouse": m,
                             "day": to_days(t, group="colonised"), "J": v})
    df = pd.DataFrame(recs)
    band = (df.groupby(["panel", "kind", "day"])["J"]
              .agg(["mean", "std"]).reset_index().fillna({"std": 0.0}))

    # Paenibacillaceae abundance is a cross-mouse mean and carries its own
    # spread, drawn as a band exactly like the clone mean does.
    ra = pd.concat([io.load_family(m).assign(mouse=m) for m in COLONISED])
    ra = ra[ra["Family"] == "Paenibacillaceae"]
    ra["day"] = to_days(ra["Time"].to_numpy(), group="colonised")
    ra = (ra.groupby("day")["abundance"].agg(["mean", "std"])
            .reset_index().fillna({"std": 0.0}))
    ra["pct"] = ra["mean"] * 100
    ra["pct_lo"] = ((ra["mean"] - ra["std"]) * 100).clip(lower=0)
    ra["pct_hi"] = ((ra["mean"] + ra["std"]) * 100).clip(upper=100)

    anc = anchors.table()
    mean_onset, mean_nadir = anc["onset_day"].mean(), anc["nadir_day"].mean()
    titles = {"C": "C   Resident (Paenibacillaceae) $\\rightarrow$ Colonizer (E. coli)",
              "D": "D   Colonizer (E. coli) $\\rightarrow$ Resident (Paenibacillaceae)"}
    ylabs = {"C": "Mean $J_{clone \\leftarrow Paeni}$",
             "D": "Mean $J_{Paeni \\leftarrow clone}$"}

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))
    for ax, panel in zip(axes, ("C", "D")):
        b = band[(band["panel"] == panel) & (band["kind"] == "clones")]
        e = band[(band["panel"] == panel) & (band["kind"] == "entero")]
        lo = float((b["mean"] - b["std"]).min()); hi = float((b["mean"] + b["std"]).max())
        pad = (hi - lo) * 0.06; lo -= pad; hi += pad

        ax.axhline(0, color="grey", ls=":", lw=0.6)
        ax.fill_between(b["day"], b["mean"] - b["std"], b["mean"] + b["std"],
                        color="#E8A33D", alpha=0.3, lw=0)
        ax.plot(b["day"], b["mean"], color="#D2691E", lw=1.6,
                label="E. coli clones (cross-mouse mean)")
        if len(e):
            ax.plot(e["day"], e["mean"], color="#8B1A1A", ls="--", lw=1.2,
                    label="Enterobacteriaceae (16S)")
        # Paenibacillaceae abundance mapped onto the panel height so the right
        # axis reads 0-80 % instead of inheriting the J scale.
        to_panel = lambda pct: lo + pct / 80.0 * (hi - lo)   # noqa: E731
        ax.fill_between(ra["day"], to_panel(ra["pct_lo"]), to_panel(ra["pct_hi"]),
                        color="#AED6F1", alpha=0.45, lw=0,
                        label="Paenibacillaceae ($\\pm$ s.d.)")
        ax.plot(ra["day"], to_panel(ra["pct"]), color="#2E86C1",
                lw=1.6, label="Paenibacillaceae abundance")
        ax.axvline(mean_nadir, color="grey", ls="--", lw=0.8)
        ax.axvline(mean_onset, color="#2E86C1", ls="-.", lw=1.0,
                   label="Paenibacillaceae onset (mean)")
        ax.set_ylim(lo, hi)
        sec = ax.secondary_yaxis(
            "right", functions=(lambda y, lo=lo, hi=hi: (y - lo) / (hi - lo) * 80.0,
                                lambda r, lo=lo, hi=hi: lo + r / 80.0 * (hi - lo)))
        sec.set_ylabel("Paenibacillaceae relative abundance (%)", color="#2E86C1")
        sec.tick_params(colors="#2E86C1")
        ax.set_xlabel("Time (days post-colonisation)")
        ax.set_ylabel(ylabs[panel])
        ax.set_title(titles[panel], fontsize=9.5)
        day_tick_marks(ax, b["day"].to_numpy(), label_every=2)
    axes[0].legend(fontsize=6.2, loc="lower right")
    fig.tight_layout()
    return style.save(fig, "fig3CD_directed_paeni_clone")


if __name__ == "__main__":
    print("Figure 3")
    panel_a(all_eigenvalues())          # full spectrum, panel A's complex plane
    panel_b(dominant_eigenvalues())     # one Re(lambda_max) per window, panel B's trend
    panels_cd()
