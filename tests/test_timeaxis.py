"""Unit tests for the two clocks.

`timeaxis` is the module where a mistake has already produced a wrong published
statistic: passing control indices through the colonised converter compresses
days 1-10 into 0.125-7 d and moves the Figure 4C test from P = 0.209 to
P = 0.019. These tests fix the mapping and, more importantly, fix the fact that
`group` cannot be omitted or guessed.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from succession.timeaxis import (day_ticks, match_nearest,      # noqa: E402
                                 sampling_slots, to_days)


# ── the mapping itself ───────────────────────────────────────────────────────

@pytest.mark.parametrize("slot,day", [
    (1, 3 / 24),      # 3 h
    (2, 6 / 24),      # 6 h
    (3, 12 / 24),     # 12 h
    (4, 1.0),         # day 1
    (5, 2.0),
    (18, 15.0),       # last cohort-2 slot
    (19, 16.0),       # last cohort-1 slot
])
def test_colonised_slots(slot, day):
    assert to_days(slot, group="colonised") == pytest.approx(day)


@pytest.mark.parametrize("index", [1, 2, 5, 10])
def test_control_index_is_the_day(index):
    assert to_days(index, group="control") == float(index)


def test_the_two_clocks_disagree_everywhere_they_matter():
    """No control index maps to the same day under both clocks."""
    idx = np.arange(1, 11)
    assert not np.any(to_days(idx, group="control")
                      == to_days(idx, group="colonised"))


def test_published_error_compresses_the_control_span():
    """The exact distortion behind Figure 4C's P = 0.019."""
    idx = np.arange(1, 11)
    right = to_days(idx, group="control")
    wrong = to_days(idx, group="colonised")
    assert (right.min(), right.max()) == (1.0, 10.0)
    assert (wrong.min(), wrong.max()) == (0.125, 7.0)


# ── group is mandatory and never guessed ─────────────────────────────────────

def test_group_has_no_default():
    with pytest.raises(TypeError):
        to_days(4)                      # type: ignore[call-arg]


@pytest.mark.parametrize("bad", ["colonized", "Colonised", "ctrl", "", None])
def test_bad_group_raises_rather_than_falling_back(bad):
    with pytest.raises(ValueError):
        to_days(4, group=bad)           # type: ignore[arg-type]


@pytest.mark.parametrize("call", [
    lambda g: to_days(4, group=g),
    lambda g: day_ticks(g, 16),
    lambda g: sampling_slots([1.0, 2.0], g),
], ids=["to_days", "day_ticks", "sampling_slots"])
def test_every_entry_point_checks_the_group(call):
    """A typo must fail loudly in each function, not just in `to_days`."""
    with pytest.raises(ValueError):
        call("colonized")


# ── shape and dtype behaviour ────────────────────────────────────────────────

def test_scalar_in_scalar_out_array_in_array_out():
    assert isinstance(to_days(4, group="colonised"), float)
    out = to_days([1, 4, 19], group="colonised")
    assert isinstance(out, np.ndarray) and out.shape == (3,)


def test_non_integer_grid_indices_convert():
    """The dense Jacobian grid carries values such as 18.9."""
    assert to_days(18.9, group="colonised") == pytest.approx(15.9)


def test_ambiguous_sub_day_indices_are_nan_not_silently_wrong():
    """A fractional index below slot 4 has no defined meaning."""
    assert np.isnan(to_days(2.5, group="colonised"))


# ── helpers ──────────────────────────────────────────────────────────────────

def test_day_ticks_start_at_3h_for_colonised_and_day_1_for_controls():
    ticks, labels = day_ticks("colonised", 16)
    assert labels[0] == "3h" and ticks[0] == pytest.approx(3 / 24)
    ticks, labels = day_ticks("control", 10)
    assert labels[0] == "1d" and ticks[0] == 1.0


def test_sampling_slots_rounds_the_dense_grid_back():
    assert list(sampling_slots([2.0, 2.4, 18.9], "colonised")) == [2, 2, 19]


def test_match_nearest_pairs_the_endpoint_but_refuses_a_far_one():
    available, values = [17.0, 18.0, 19.0], [170.0, 180.0, 190.0]
    # 18.9 is within the 0.6-slot tolerance of slot 19
    assert match_nearest(18.9, available, values) == 190.0
    # 5.0 is nowhere near any of them
    assert match_nearest(5.0, available, values) is None


def test_match_nearest_on_empty_input():
    assert match_nearest(1.0, [], []) is None
