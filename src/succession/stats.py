"""Statistics, with the caveats attached where they matter.

Every timecourse here is autocorrelated, so asymptotic p-values are
anticonservative. `spearman` reports one anyway because the manuscript does, but
`circular_shift_p` gives a null that preserves within-series autocorrelation and
should be preferred when a claim rests on significance.
"""

from __future__ import annotations

import numpy as np
from scipy import stats as _st


def spearman(x, y):
    """Spearman rho with its asymptotic p-value and n."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 4:
        return np.nan, np.nan, int(ok.sum())
    rho, p = _st.spearmanr(x[ok], y[ok])
    return float(rho), float(p), int(ok.sum())


def mann_whitney(a, b, alternative: str = "greater"):
    """Mann-Whitney U. Default tests whether `a` is stochastically greater.

    Figure 4C uses this one-sided form: control interactions greater, i.e. less
    negative, than colonised.
    """
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    U, p = _st.mannwhitneyu(a[np.isfinite(a)], b[np.isfinite(b)],
                            alternative=alternative)
    return float(U), float(p)


def circular_shift_p(x, y, n_perm: int = 1000, seed: int = 0) -> float:
    """Permutation p for a Spearman rho, preserving autocorrelation.

    Shuffling breaks the serial structure and so tests the wrong null for a
    timecourse. Circularly shifting one series keeps it.
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if x.size < 4:
        return np.nan
    obs = abs(_st.spearmanr(x, y)[0])
    rng = np.random.default_rng(seed)
    hits = sum(abs(_st.spearmanr(x, np.roll(y, int(k)))[0]) >= obs
               for k in rng.integers(1, x.size, n_perm))
    return (hits + 1) / (n_perm + 1)
