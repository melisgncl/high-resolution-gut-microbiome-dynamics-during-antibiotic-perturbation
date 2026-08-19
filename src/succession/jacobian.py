"""Time-varying community Jacobian from a sliding window.

The estimator
-------------
State variables ``z`` are log10 abundances on a dense grid (0.1 sampling slots,
ten points per slot). For an evaluation time ``t`` and window width ``w``:

    J[i <- j](t) = cov( dz_i/dt , z_j )   over the half-open window (t - w, t]

``i`` is the target whose rate responds, ``j`` is the driver whose abundance
acts. The window is half-open, so with ``w = 5`` it contains exactly five
sampling slots.

This is a covariance, not a normalised Jacobian: the textbook element would be
``cov(dz_i/dt, z_j) / var(z_j)``. The denominator is dropped here to match the
published analysis. It is not, as is sometimes claimed, close to unity - during
the collapse phase the resident drivers are near-constant, ``var(z_j)`` approaches
zero and the normalised slope diverges. Dropping it is a defensible choice for
stability, but for that reason and not because it makes no difference.

Which evaluation times are used
-------------------------------
Evaluation times come from the stored Jacobian timeseries, then a warm-up rule
discards those whose window runs off the start of the series. See
`config.WARMUP_TOLERANCE`, which documents the effect on the published numbers.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline

from . import io
from .config import MIN_POINTS, PSEUDOCOUNT, WARMUP_TOLERANCE, WINDOW

_STD_FLOOR = 1e-12


@dataclass
class State:
    """The state matrix a Jacobian is estimated from."""

    grid: np.ndarray          #: dense sampling-index grid
    Z: np.ndarray             #: (n_grid, n_species) log10 abundances
    species: list[str]        #: column names, clones first then 16S taxa

    @property
    def start(self) -> float:
        return float(self.grid[0])

    def index_of(self, name: str) -> int | None:
        return self.species.index(name) if name in self.species else None

    def clone_indices(self) -> list[int]:
        return [k for k, s in enumerate(self.species) if s.startswith("C")
                and s[1:].isdigit()]


def build_state(mouse: str, pseudocount: float = PSEUDOCOUNT) -> State:
    """Assemble clonal-cluster and 16S trajectories onto one dense grid.

    Clone trajectories are already LOESS-smoothed log10 frequencies on the dense
    grid. The 16S taxa are measured once per slot, so they are interpolated onto
    the same grid with a cubic spline in log10 space after adding `pseudocount`.
    Only the taxa that appear in the stored Jacobian are used, so the species
    set matches the published model.
    """
    clones = io.load_clone_loess(mouse, relabel=False)
    if clones.empty:                       # controls carry no barcode dimension
        idx, taxa, _ = io.load_taxa_matrix(mouse)
        grid = np.arange(idx[0], idx[-1] + 0.05, 0.1)
        cols, names = [], []
    else:
        wide = clones.pivot(index="index", columns="clone", values="log10_freq")
        wide = wide.sort_index()
        grid = wide.index.to_numpy(float)
        names = [f"C{c}" for c in sorted(int(n[1:]) for n in wide.columns)]
        cols = [wide[n].to_numpy(float) for n in names]

    used = set(io.load_jacobian_timeseries(mouse)["driver"])
    idx, taxa, mat = io.load_taxa_matrix(mouse)
    for k, name in enumerate(taxa):
        if name not in used:
            continue
        spline = CubicSpline(idx, np.log10(mat[:, k] + pseudocount),
                             extrapolate=False)
        vals = spline(grid)
        bad = np.isnan(vals)
        if bad.any():                      # nearest-valid fill at the edges
            good = np.flatnonzero(~bad)
            vals[bad] = np.interp(np.flatnonzero(bad), good, vals[good])
        cols.append(vals)
        names.append(name)

    return State(grid=grid, Z=np.column_stack(cols), species=names)


def evaluation_times(
    mouse: str,
    state: State,
    window: int = WINDOW,
    warmup_tolerance: float | None = WARMUP_TOLERANCE,
) -> np.ndarray:
    """Evaluation times surviving the warm-up rule.

    ``warmup_tolerance=None`` keeps every evaluation time, which is what the
    rule's docstring in `config` describes the consequences of.
    """
    times = io.evaluation_indices(mouse)
    if warmup_tolerance is None:
        return times
    return times[times - window >= state.start - warmup_tolerance]


def _window_slice(grid: np.ndarray, t: float, window: int):
    hi = int(np.searchsorted(grid, t, side="right")) - 1
    lo = max(0, int(np.searchsorted(grid, t - window, side="left")))
    return (lo, hi) if hi - lo >= MIN_POINTS - 1 else None


def offdiagonal(
    state: State,
    times,
    window: int = WINDOW,
) -> dict[float, np.ndarray]:
    """Every off-diagonal J[i <- j] at each evaluation time.

    Returns {time: array of coefficients}. Pairs whose driver abundance or
    target derivative is constant inside the window are dropped, matching the
    published guard of 1e-12 on the standard deviation.
    """
    step = state.grid[1] - state.grid[0]
    dZ = np.gradient(state.Z, step, axis=0)
    n = state.Z.shape[1]

    out: dict[float, np.ndarray] = {}
    for t in np.atleast_1d(times):
        sl = _window_slice(state.grid, float(t), window)
        if sl is None:
            continue
        lo, hi = sl
        zw, dw = state.Z[lo:hi + 1], dZ[lo:hi + 1]
        sd_z, sd_d = zw.std(axis=0), dw.std(axis=0)
        vals = [np.cov(dw[:, i], zw[:, j])[0, 1]
                for i in range(n) for j in range(n)
                if i != j and sd_z[j] >= _STD_FLOOR and sd_d[i] >= _STD_FLOOR]
        if vals:
            out[float(t)] = np.asarray(vals)
    return out


def offdiagonal_by_slot(
    mouse: str,
    state: State,
    n_steps: int = 5,
    min_points: int = 3,
) -> dict[float, np.ndarray]:
    """Off-diagonal J at each sampling slot, windowed by SLOT COUNT.

    This is the estimator behind Figure 2B, ported from the manuscript's
    `12_ridgeline_jacobian.py`, and it is deliberately not the same as
    `offdiagonal`:

    * The window for slot ``k`` spans the ``n_steps`` most recent *samples*,
      ``slots[k - n_steps + 1] .. slots[k]``, rather than a fixed numeric width
      on the index grid. Early on, when samples are hours rather than days
      apart, those are very different spans.
    * It **expands** at the start (``max(0, k - n_steps + 1)``) instead of being
      discarded by a warm-up rule, so the earliest windows are shorter but still
      drawn. Figure 2B describes a distribution; it is not computing the
      correlation that needs a consistent width.

    `offdiagonal` remains the panel-C estimator and reproduces the published
    rho = 0.73 at n = 113. Using this one there would change those numbers.

    Returns {sampling index: array of off-diagonal coefficients}.
    """
    slots, _, _ = io.load_taxa_matrix(mouse)
    step = state.grid[1] - state.grid[0]
    dZ = np.gradient(state.Z, step, axis=0)
    n = state.Z.shape[1]

    out: dict[float, np.ndarray] = {}
    for k in range(1, len(slots)):          # k=0 has no preceding sample
        lo_slot = slots[max(0, k - (n_steps - 1))]
        i_lo = int(np.searchsorted(state.grid, lo_slot, side="left"))
        i_hi = int(np.searchsorted(state.grid, slots[k], side="right"))
        if i_hi - i_lo < min_points:
            continue

        zw, dw = state.Z[i_lo:i_hi], dZ[i_lo:i_hi]
        sd_z, sd_d = zw.std(axis=0), dw.std(axis=0)
        vals = [np.cov(dw[:, i], zw[:, j])[0, 1]
                for i in range(n) for j in range(n)
                if i != j and sd_z[j] >= _STD_FLOOR and sd_d[i] >= _STD_FLOOR]
        if vals:
            out[float(slots[k])] = np.asarray(vals)
    return out


def directed(
    state: State,
    times,
    target: str | int,
    driver: str | int,
    window: int = WINDOW,
) -> dict[float, float]:
    """A single directed coefficient J[target <- driver] over time."""
    i = state.index_of(target) if isinstance(target, str) else target
    j = state.index_of(driver) if isinstance(driver, str) else driver
    if i is None or j is None:
        return {}

    step = state.grid[1] - state.grid[0]
    dZ = np.gradient(state.Z, step, axis=0)
    out = {}
    for t in np.atleast_1d(times):
        sl = _window_slice(state.grid, float(t), window)
        if sl is None:
            continue
        lo, hi = sl
        zw, dw = state.Z[lo:hi + 1], dZ[lo:hi + 1]
        if zw[:, j].std() >= _STD_FLOOR and dw[:, i].std() >= _STD_FLOOR:
            out[float(t)] = float(np.cov(dw[:, i], zw[:, j])[0, 1])
    return out


def directed_group(
    state: State,
    times,
    targets,
    drivers,
    window: int = WINDOW,
) -> dict[float, float]:
    """Mean of J[target <- driver] over every target/driver pair given.

    Figure 3C-D averages over a mouse's clonal clusters before pooling across
    mice, so the per-mouse value is a mean over clones, not over mice.
    """
    ti = [state.index_of(t) if isinstance(t, str) else t for t in targets]
    dj = [state.index_of(d) if isinstance(d, str) else d for d in drivers]
    ti = [k for k in ti if k is not None]
    dj = [k for k in dj if k is not None]
    if not ti or not dj:
        return {}

    step = state.grid[1] - state.grid[0]
    dZ = np.gradient(state.Z, step, axis=0)
    out = {}
    for t in np.atleast_1d(times):
        sl = _window_slice(state.grid, float(t), window)
        if sl is None:
            continue
        lo, hi = sl
        zw, dw = state.Z[lo:hi + 1], dZ[lo:hi + 1]
        vals = [np.cov(dw[:, i], zw[:, j])[0, 1]
                for i in ti for j in dj
                if i != j and zw[:, j].std() >= _STD_FLOOR
                and dw[:, i].std() >= _STD_FLOOR]
        if vals:
            out[float(t)] = float(np.mean(vals))
    return out


def summarise(offdiag: dict[float, np.ndarray]) -> "pd.DataFrame":  # noqa: F821
    """Per-evaluation-time summaries of an off-diagonal distribution.

    ``mean_negative`` is the mean of only the negative coefficients - the
    inhibitory-strength measure Figure 2C correlates against diversity. It is
    not the same as ``mean``, which averages inhibitory and facilitative
    together and lands near zero.
    """
    import pandas as pd

    rows = []
    for t, v in sorted(offdiag.items()):
        neg, pos = v[v < 0], v[v > 0]
        rows.append({
            "index": t,
            "n_pairs": v.size,
            "mean": v.mean(),
            "mean_negative": neg.mean() if neg.size else np.nan,
            "mean_positive": pos.mean() if pos.size else np.nan,
            "frac_negative": float((v < 0).mean()),
        })
    return pd.DataFrame(rows)



def full(
    state: State,
    times,
    window: int = WINDOW,
) -> dict[float, np.ndarray]:
    """The complete n x n Jacobian at each evaluation time, diagonal included.

    ``offdiagonal`` drops the diagonal because Figure 2C is about interactions
    between taxa. Eigenvalues are a different question: J[i <- i] carries the
    self-regulation term that sets where the spectrum sits, and omitting it
    moves every eigenvalue. Use this function, not ``offdiagonal``, for anything
    spectral.

    Pairs failing the ``_STD_FLOOR`` guard are set to 0.0 rather than dropped -
    a matrix cannot have holes. Zero is the honest fill: no detectable response.
    """
    step = state.grid[1] - state.grid[0]
    dZ = np.gradient(state.Z, step, axis=0)
    n = state.Z.shape[1]

    out: dict[float, np.ndarray] = {}
    for t in np.atleast_1d(times):
        sl = _window_slice(state.grid, float(t), window)
        if sl is None:
            continue
        lo, hi = sl
        zw, dw = state.Z[lo:hi + 1], dZ[lo:hi + 1]
        sd_z, sd_d = zw.std(axis=0), dw.std(axis=0)
        J = np.zeros((n, n))
        for i in range(n):
            if sd_d[i] < _STD_FLOOR:
                continue
            for j in range(n):
                if sd_z[j] >= _STD_FLOOR:
                    J[i, j] = np.cov(dw[:, i], zw[:, j])[0, 1]
        if np.isfinite(J).all():
            out[float(t)] = J
    return out


def corr(
    state: State,
    times,
    window: int = WINDOW,
) -> dict[float, np.ndarray]:
    """Scale-free interaction matrix - the amplitude-free counterpart of J.

        R[i <- j] = cov(dz_i/dt, z_j) / ( sd(dz_i/dt) sd(z_j) )

    Why this exists: ``cov(dz_i/dt, z_j)`` estimates A*C, the interaction matrix
    times the state covariance, not A alone. As succession proceeds the
    trajectories flatten, C shrinks, and |J| shrinks with it even if A never
    changes. R divides the amplitude out and is a Pearson correlation, bounded
    in [-1, 1] and invariant to rescaling either series.

    This is not cosmetic. Against TIME, mean |J| gives rho = -0.892 (window 5,
    n = 113, 8/8 mice) while mean |R| gives rho = +0.152, p = 0.11, significant
    in 2/8 mice with opposite signs. Against DIVERSITY the picture differs: the
    prevalence statistic ``frac_positive`` survives at every window tested,
    while mean |R| flips sign between window 5 and window 10. See
    ``CORRECTIONS.md``.

    Diagonal entries are NaN: a self term is not an interaction.
    """
    step = state.grid[1] - state.grid[0]
    dZ = np.gradient(state.Z, step, axis=0)
    n = state.Z.shape[1]

    out: dict[float, np.ndarray] = {}
    for t in np.atleast_1d(times):
        sl = _window_slice(state.grid, float(t), window)
        if sl is None:
            continue
        lo, hi = sl
        zw, dw = state.Z[lo:hi + 1], dZ[lo:hi + 1]
        sd_z, sd_d = zw.std(axis=0, ddof=1), dw.std(axis=0, ddof=1)
        R = np.full((n, n), np.nan)
        for i in range(n):
            if sd_d[i] < _STD_FLOOR:
                continue
            for j in range(n):
                if i != j and sd_z[j] >= _STD_FLOOR:
                    R[i, j] = (np.cov(dw[:, i], zw[:, j])[0, 1]
                               / (sd_d[i] * sd_z[j]))
        out[float(t)] = R
    return out


def summarise_corr(corrmats: dict[float, np.ndarray]) -> "pd.DataFrame":  # noqa: F821
    """Per-time summaries of a scale-free matrix, mirroring ``summarise``.

    ``frac_positive`` is the statistic that survives every robustness check:
    it rises with diversity in all four (cohort, window) configurations tested,
    and it is what licenses the claim that inhibitory links become less
    PREVALENT. ``mean_absolute`` is the magnitude claim's amplitude-free test,
    and it does not survive - see the module docstring of ``corr``.
    """
    import pandas as pd

    rows = []
    for t, R in sorted(corrmats.items()):
        v = R[np.isfinite(R)]
        if v.size == 0:
            continue
        neg = v[v < 0]
        rows.append({
            "index": t,
            "n_pairs": v.size,
            "mean": float(v.mean()),
            "mean_absolute": float(np.abs(v).mean()),
            "mean_negative": float(neg.mean()) if neg.size else np.nan,
            "frac_positive": float((v > 0).mean()),
        })
    return pd.DataFrame(rows)


def effective_rank(J: np.ndarray, tol: float = 1e-10) -> int:
    """Number of singular values above ``tol``.

    Diagnostic for the confound in ``eigenvalues``: a covariance estimator has
    rank <= samples in window, so under an EXPANDING window rank is a
    deterministic function of time and the null-space zeros sit at Re(lambda) = 0,
    above the mostly-negative real spectrum. Under the sliding window used by
    ``full`` the rank is constant and this diagnostic should be flat.
    """
    return int((np.linalg.svd(J, compute_uv=False) > tol).sum())


def dominant_eigenvalue(mouse: str, window: int = WINDOW) -> "pd.DataFrame":  # noqa: F821
    """Re(lambda_max) per time point, from sliding-window matrices.

    The corrected counterpart of ``eigenvalues``. Two differences, each of which
    moves the published number:

    * matrices come from ``full`` (sliding window) rather than from the stored
      expanding-window files, so matrix rank no longer tracks time;
    * one row per time point, not one per eigenvalue, so a mouse contributes
      its number of windows rather than windows x species to any correlation.

    The published Figure 3B trend (rho = -0.43, n = 1545) is both confounded and
    pseudoreplicated. See ``CORRECTIONS.md``.
    """
    import pandas as pd

    state = build_state(mouse)
    times = evaluation_times(mouse, state, window=window)
    rows = []
    for t, J in sorted(full(state, times, window).items()):
        lam = np.linalg.eigvals(J)
        k = int(np.argmax(lam.real))
        rows.append({"mouse": mouse, "index": t,
                     "re_max": float(lam[k].real), "im_max": float(lam[k].imag),
                     "rank": effective_rank(J), "n_species": J.shape[0]})
    return pd.DataFrame(rows)


def eigenvalues(mouse: str) -> "pd.DataFrame":  # noqa: F821
    """Eigenvalues of every stored per-window Jacobian.

    The community is locally stable where every Re(lambda) < 0.

    RETAINED FOR PROVENANCE. This reads the stored matrices in ``data/jacobian``,
    which were estimated with an EXPANDING window, and returns one row per
    eigenvalue. Both properties affect the published Figure 3B trend. Use
    ``dominant_eigenvalue`` for new analysis; see ``CORRECTIONS.md``.
    """
    import pandas as pd

    rows = []
    mats = io.load_jacobian_matrices(mouse)
    idx = np.array([m[0] for m in mats], dtype=float)
    span = idx.max() - idx.min() + 1e-9
    for (t, J, _), norm in zip(mats, (idx - idx.min()) / span):
        for lam in np.linalg.eigvals(J):
            rows.append({"mouse": mouse, "index": t, "relative_time": norm,
                         "re": lam.real, "im": lam.imag})
    return pd.DataFrame(rows)
