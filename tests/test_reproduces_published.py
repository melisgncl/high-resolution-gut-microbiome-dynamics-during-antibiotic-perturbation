"""The estimator must reproduce the published statistics, or nothing else matters.

Numbers below are taken from the submitted manuscript. If a change to
`succession.jacobian` breaks one of these, the change is wrong.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from succession import diversity, io, jacobian, stats          # noqa: E402
from succession.config import COLONISED, WINDOW                # noqa: E402
from succession.timeaxis import match_nearest                  # noqa: E402


def figure2c():
    """Mean negative J against 16S Hill 1D, the published way."""
    xs, ys, per_mouse = [], [], {}
    for m in COLONISED:
        state = jacobian.build_state(m)
        times = jacobian.evaluation_times(m, state, window=WINDOW)
        summ = jacobian.summarise(jacobian.offdiagonal(state, times, WINDOW))
        q1 = diversity.hill_q1_from_taxa(m)

        mx, my = [], []
        for _, row in summ.iterrows():
            if not np.isfinite(row["mean_negative"]):
                continue
            q = match_nearest(row["index"], q1["index"].to_numpy(),
                              q1["q1"].to_numpy())
            if q is None:
                continue
            mx.append(q)
            my.append(row["mean_negative"])
        per_mouse[m] = stats.spearman(mx, my)
        xs += mx
        ys += my
    return stats.spearman(xs, ys), per_mouse


def test_figure2c_pooled():
    (rho, p, n), _ = figure2c()
    assert n == 113, f"published n = 113, got {n}"
    assert rho == pytest.approx(0.73, abs=0.01), f"published rho = 0.73, got {rho:.3f}"
    assert p < 1e-18


def test_figure2c_per_mouse_counts():
    _, per = figure2c()
    assert [per[m][2] for m in COLONISED] == [15, 15, 15, 15, 14, 13, 13, 13]


def test_figure2c_significance_pattern():
    _, per = figure2c()
    sig = [m for m in COLONISED if per[m][1] < 0.05]
    assert len(sig) == 7, f"published: 7 of 8 significant, got {len(sig)}"
    assert "m8" not in sig, "m8 is the published exception"


def test_figure3b_eigenvalue_trend():
    import pandas as pd
    from succession.timeaxis import to_days

    eig = pd.concat([jacobian.eigenvalues(m) for m in COLONISED])
    day = to_days(eig["index"].to_numpy(), group="colonised")
    rho, p, n = stats.spearman(day, eig["re"].to_numpy())
    assert n == 1545
    assert rho == pytest.approx(-0.43, abs=0.01), f"published rho = -0.43, got {rho:.3f}"
    assert p < 1e-60


def test_warmup_rule_sets_sample_size():
    """n must fall by exactly one evaluation per mouse per extra window step."""
    totals = {}
    for w in (1, 4, 5, 10):
        n = 0
        for m in COLONISED:
            state = jacobian.build_state(m)
            n += len(jacobian.evaluation_times(m, state, window=w))
        totals[w] = n
    assert totals == {1: 129, 4: 121, 5: 113, 10: 73}


def test_control_clock_is_not_the_colonised_clock():
    from succession.timeaxis import to_days
    assert to_days(10, group="control") == 10.0
    assert to_days(10, group="colonised") == 7.0
