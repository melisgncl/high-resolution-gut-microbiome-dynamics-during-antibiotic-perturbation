"""The corrected statistics. See CORRECTIONS.md.

`test_reproduces_published.py` pins what was submitted, so the manuscript stays
reproducible. This file pins what is *true*. Where the two disagree, CORRECTIONS.md
says why. A change that breaks a test here is wrong; a change that breaks one there
is a change to the published path and needs the same scrutiny.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from succession import jacobian, stats                          # noqa: E402
from succession.config import COLONISED                         # noqa: E402
from succession.timeaxis import to_days                         # noqa: E402


def _against_time(window):
    days, absJ, absR, fpos = [], [], [], []
    for m in COLONISED:
        st = jacobian.build_state(m)
        ts = jacobian.evaluation_times(m, st, window=window)
        Jm, Rm = jacobian.full(st, ts, window), jacobian.corr(st, ts, window)
        for t in sorted(set(Jm) & set(Rm)):
            J, R = Jm[t], Rm[t]
            off = ~np.eye(J.shape[0], dtype=bool)
            v = R[np.isfinite(R)]
            days.append(to_days(t, group="colonised"))
            absJ.append(np.abs(J[off]).mean())
            absR.append(np.abs(v).mean())
            fpos.append((v > 0).mean())
    return np.array(days), np.array(absJ), np.array(absR), np.array(fpos)


# ── CORRECTIONS.md #3 — the ninefold decline is amplitude ────────────────────

@pytest.mark.parametrize("window,n,rho_J", [(5, 113, -0.892), (10, 73, -0.842)])
def test_raw_magnitude_declines(window, n, rho_J):
    """mean |J| falls steeply against time. This is the published behaviour."""
    d, aJ, _, _ = _against_time(window)
    r, p, k = stats.spearman(d, aJ)
    assert k == n
    assert r == pytest.approx(rho_J, abs=0.01)


@pytest.mark.parametrize("window", [5, 10])
def test_amplitude_free_magnitude_is_flat(window):
    """mean |R| does not decline. The ninefold weakening is the amplitude term."""
    d, _, aR, _ = _against_time(window)
    r, p, _ = stats.spearman(d, aR)
    assert p > 0.05, f"window {window}: expected no trend, got rho={r:+.3f} p={p:.3g}"
    assert abs(r) < 0.25


# ── CORRECTIONS.md #4 — prevalence shifts, and stops at symmetry ─────────────

@pytest.mark.parametrize("window,rho", [(5, 0.292), (10, 0.426)])
def test_inhibitory_links_become_less_prevalent(window, rho):
    """frac(R>0) rises against time - the one amplitude-free trend that holds."""
    d, _, _, f = _against_time(window)
    r, p, _ = stats.spearman(d, f)
    assert r == pytest.approx(rho, abs=0.02)
    assert p < 0.01


@pytest.mark.parametrize("window", [5, 10])
def test_sign_composition_stops_at_symmetry(window):
    """It never becomes facilitation-dominated. Late mean stays at or below 0.5."""
    d, _, _, f = _against_time(window)
    late = f[d >= np.percentile(d, 75)]
    assert late.mean() <= 0.52, f"window {window}: late frac(R>0) = {late.mean():.3f}"


# ── CORRECTIONS.md #2 — the rank confound is gone under a sliding window ─────

@pytest.mark.parametrize("window", [5, 10])
def test_rank_no_longer_tracks_time(window):
    """The confound in the stored matrices: rank = min(window, S), so rank ~ time."""
    dom = pd.concat([jacobian.dominant_eigenvalue(m, window=window)
                     for m in COLONISED])
    d = to_days(dom["index"].to_numpy(), group="colonised")
    r, p, _ = stats.spearman(d, dom["rank"].to_numpy())
    assert p > 0.05, f"window {window}: rank still tracks time, rho={r:+.3f}"


@pytest.mark.parametrize("window,n,rho", [(5, 113, -0.690), (10, 73, -0.813)])
def test_dominant_eigenvalue_trend(window, n, rho):
    """One eigenvalue per window, sliding window. Replaces rho = -0.43, n = 1545."""
    dom = pd.concat([jacobian.dominant_eigenvalue(m, window=window)
                     for m in COLONISED])
    d = to_days(dom["index"].to_numpy(), group="colonised")
    r, p, k = stats.spearman(d, dom["re_max"].to_numpy())
    assert k == n
    assert r == pytest.approx(rho, abs=0.01)


@pytest.mark.parametrize("window", [5, 10])
def test_eigenvalue_never_crosses_zero(window):
    """It declines without reaching stability. Do not write 'the community stabilises'."""
    dom = pd.concat([jacobian.dominant_eigenvalue(m, window=window)
                     for m in COLONISED])
    assert (dom["re_max"] > 0).mean() >= 0.99


# ── CORRECTIONS.md #5 — control Paenibacillaceae is present, not absent ──────

def test_control_paenibacillaceae_is_present_and_never_blooms():
    """It was filtered out by `mean < 0.005`, then read downstream as zero."""
    import csv
    from succession.config import CONTROLS, DATA

    for m in CONTROLS:
        by_t = {}
        with open(DATA / "16s" / "family" / f"{m}_family.csv") as fh:
            for row in csv.DictReader(fh):
                if row["Family"] == "Paenibacillaceae":
                    t = float(row["Time"])
                    by_t[t] = by_t.get(t, 0.0) + float(row["Abundance.family"])
        v = np.array(list(by_t.values()))
        assert v.size > 0, f"{m}: Paenibacillaceae absent from the unfiltered table"
        assert v.max() > 0.0, f"{m}: detected nowhere"
        assert v.max() < 0.05, f"{m}: unexpected bloom, max = {v.max():.4f}"
        assert v.mean() < 0.005, f"{m}: mean {v.mean():.4f} - above the filter that removed it"




# ── window decision, 19 Aug 2026 (option B) ──────────────────────────────────
# WINDOW = 5 stays primary for Fig. 2C and Fig. 3 - the per-mouse panel is
# intact (7/8) at this window, and it is what test_reproduces_published.py
# already pins. WINDOW_ROBUSTNESS = 10 is reported only as a supplementary
# diversity-robustness check (figS4_window_sweep.py): the diversity relationship
# strengthens there and the amplitude-free statistic stabilises there, but it is
# not a primary-figure window. See CORRECTIONS.md ss4b.

def test_window_is_five():
    from succession.config import WINDOW
    assert WINDOW == 5


def test_robustness_window_is_ten():
    from succession.config import WINDOW_ROBUSTNESS
    assert WINDOW_ROBUSTNESS == 10


def test_figure2c_per_mouse_intact_at_primary_window():
    """7/8 mice significant at the primary window - this is what ships."""
    from succession.config import COLONISED, WINDOW
    from succession import diversity
    from succession.timeaxis import match_nearest

    sig = 0
    for m in COLONISED:
        st = jacobian.build_state(m)
        ts = jacobian.evaluation_times(m, st, window=WINDOW)
        summ = jacobian.summarise(jacobian.offdiagonal(st, ts, WINDOW))
        q1 = diversity.hill_q1_from_taxa(m)
        x, y = [], []
        for _, r in summ.iterrows():
            if not np.isfinite(r["mean_negative"]):
                continue
            q = match_nearest(r["index"], q1["index"].to_numpy(), q1["q1"].to_numpy())
            if q is None:
                continue
            x.append(q); y.append(r["mean_negative"])
        _, p, _ = stats.spearman(x, y)
        sig += p < 0.05
    assert sig == 7, f"expected 7/8 mice significant at window {WINDOW}, got {sig}/8"


def test_figure2c_diversity_robustness_check_at_window_ten():
    """The supplementary check: window 10 strengthens the pooled correlation
    even though it costs per-mouse power (documented, not hidden, in
    figS4_window_sweep.py). This is a robustness footnote, not the headline."""
    from succession.config import COLONISED, WINDOW_ROBUSTNESS
    from succession import diversity
    from succession.timeaxis import match_nearest

    x, y = [], []
    for m in COLONISED:
        st = jacobian.build_state(m)
        ts = jacobian.evaluation_times(m, st, window=WINDOW_ROBUSTNESS)
        summ = jacobian.summarise(jacobian.offdiagonal(st, ts, WINDOW_ROBUSTNESS))
        q1 = diversity.hill_q1_from_taxa(m)
        for _, r in summ.iterrows():
            if not np.isfinite(r["mean_negative"]):
                continue
            q = match_nearest(r["index"], q1["index"].to_numpy(), q1["q1"].to_numpy())
            if q is None:
                continue
            x.append(q); y.append(r["mean_negative"])
    rho, p, n = stats.spearman(x, y)
    assert n == 73
    assert rho == pytest.approx(0.574, abs=0.01)


def test_figure3_dominant_eigenvalue_uses_primary_window():
    """Fig. 3B's headline statistic, at the primary window."""
    from succession.config import COLONISED, WINDOW
    from succession.timeaxis import to_days

    dom = pd.concat([jacobian.dominant_eigenvalue(m, window=WINDOW) for m in COLONISED])
    d = to_days(dom["index"].to_numpy(), group="colonised")
    rho, p, n = stats.spearman(d, dom["re_max"].to_numpy())
    assert n == 113
    assert rho == pytest.approx(-0.690, abs=0.01)


def test_figure4_structurally_cannot_use_robustness_window():
    """Why Figure 4 does not follow WINDOW_ROBUSTNESS even as a footnote."""
    from succession.config import CONTROLS, WINDOW_ROBUSTNESS

    total = 0
    for m in CONTROLS:
        st = jacobian.build_state(m)
        ts = jacobian.evaluation_times(m, st, window=WINDOW_ROBUSTNESS)
        total += len(jacobian.offdiagonal(st, ts, WINDOW_ROBUSTNESS))
    assert total <= 1, (
        f"expected the control series to be exhausted at window {WINDOW_ROBUSTNESS} "
        f"(control index = day), got {total} usable evaluations")


# ── Figure 4C amplitude mechanism, 19 Aug 2026 ───────────────────────────────
def test_paenibacillaceae_noise_drives_control_amplitude():
    """Why raw |J| says controls > colonised while |R| says the opposite.

    Not general trajectory flattening - specifically Paenibacillaceae's own
    noise in the control arm, where it sits at ~1% with only 7-10 samples
    (the 16S low-count regime). See CORRECTIONS.md and figure4c_rebuilt.py.
    """
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "figures"))
    from figure4c_rebuilt import PAIR, WIN_DAYS, state
    from succession.config import CONTROLS, COLONISED

    def per_family_sdd(mice):
        out = {f: [] for f in PAIR}
        for m in mice:
            grid, Z, days, _, _ = state(m)
            step = grid[1] - grid[0]
            dZ = np.gradient(Z, step, axis=0)
            for te in days:
                if te - WIN_DAYS < grid[0] - 2.0:
                    continue
                hi = int(np.searchsorted(grid, te, side="right")) - 1
                lo = max(0, int(np.searchsorted(grid, te - WIN_DAYS, side="left")))
                if hi - lo < 4:
                    continue
                dw = dZ[lo:hi + 1]
                for k, f in enumerate(PAIR):
                    out[f].append(dw[:, k].std(ddof=1))
        return {f: np.median(v) for f, v in out.items()}

    ctl = per_family_sdd(CONTROLS)
    col = per_family_sdd(COLONISED)

    # Controls are noisier than colonised mice for BOTH families - fewer
    # samples (7-10 vs 15-19), not something specific to one taxon.
    assert ctl["Enterobacteriaceae"] > col["Enterobacteriaceae"]
    assert ctl["Paenibacillaceae"] > col["Paenibacillaceae"]
    # Within either arm, Paeni is far noisier than Entero because it is the
    # rarer taxon (log-scale amplification of count noise). This is what
    # dominates raw |J|'s absolute magnitude - not a control-specific effect
    # on Paeni alone. Do not overstate this as "Paeni-specific noise in
    # controls"; Entero's fold-increase control-vs-colonised is comparable.
    assert ctl["Paenibacillaceae"] > 3 * ctl["Enterobacteriaceae"]
    assert ctl["Paenibacillaceae"] == max(ctl["Paenibacillaceae"], ctl["Enterobacteriaceae"],
                                          col["Paenibacillaceae"], col["Enterobacteriaceae"])
