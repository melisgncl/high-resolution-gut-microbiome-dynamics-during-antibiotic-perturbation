"""Conversion between sampling index and real time.

The two experimental groups are indexed differently, and conflating them is the
single easiest way to produce a wrong statistic from this dataset.

Colonised mice (m1-m8) were sampled three times on the first day and daily
afterwards, so the index is a *sampling slot*:

    slot 1 -> 3 h        slot 2 -> 6 h        slot 3 -> 12 h
    slot 4 -> day 1      slot k -> day k - 3      (k >= 4)

Antibiotic-only control mice (c_m1-c_m4) were sampled once a day from day 1,
with no sub-day timepoints, so the index *is* the day:

    index k -> day k

Applying the colonised mapping to a control index compresses day 1..10 into
0.125..7 d, so every axis and every window boundary built from it is wrong.

The damage is indirect and worth stating precisely. The clock does not change
any measured value; it changes which colonised evaluations fall inside a
comparison window that is defined relative to the control span. In Figure 4C
that silently moved the colonised set between 53 and 113 points, and with it
the one-sided Mann-Whitney P between 0.019 and 0.27. Once the window is named
explicitly (`figure4.COLONISED_WINDOW`), the two clocks give *identical* P at
every window width - which is the real lesson: an open-ended window let a
labelling error masquerade as a result.

`group` is therefore a required argument everywhere in this module. There is no
default, and there is no function that guesses.
"""

from __future__ import annotations

from typing import Iterable, Literal, Sequence

import numpy as np

Group = Literal["colonised", "control"]

#: Sub-day slots for the colonised cohorts, in days.
_SUBDAY = {1: 3 / 24, 2: 6 / 24, 3: 12 / 24}

#: Slot index of day 1 for the colonised cohorts; slot k maps to day k - OFFSET.
_OFFSET = 3


def _check_group(group: str) -> Group:
    if group not in ("colonised", "control"):
        raise ValueError(
            f"group must be 'colonised' or 'control', got {group!r}. "
            "This argument is mandatory: the two groups use different clocks."
        )
    return group  # type: ignore[return-value]


def to_days(index, group: Group):
    """Convert a sampling index to days post-colonisation.

    Parameters
    ----------
    index : scalar or array-like
        Sampling index as stored in the tables. May be non-integer: the dense
        Jacobian grid carries values such as 18.9.
    group : {'colonised', 'control'}
        Required. See the module docstring.

    Examples
    --------
    >>> to_days(1, group="colonised")
    0.125
    >>> to_days(4, group="colonised")
    1.0
    >>> to_days(19, group="colonised")
    16.0
    >>> to_days(1, group="control")
    1.0
    >>> to_days(10, group="control")
    10.0
    """
    _check_group(group)
    scalar = np.isscalar(index)
    idx = np.atleast_1d(np.asarray(index, dtype=float))

    if group == "control":
        out = idx.copy()
    else:
        out = np.where(idx >= _OFFSET + 1, idx - _OFFSET, np.nan)
        for slot, day in _SUBDAY.items():
            out = np.where(np.isclose(idx, slot), day, out)
        # Non-integer indices below slot 4 (none in the shipped data) would be
        # ambiguous rather than silently wrong, so they stay NaN.

    return float(out[0]) if scalar else out


def day_ticks(group: Group, upto: float) -> tuple[list[float], list[str]]:
    """Axis ticks in days, labelled the way the manuscript labels them."""
    _check_group(group)
    if group == "control":
        ticks = [t for t in range(1, int(upto) + 1, 2)]
        return [float(t) for t in ticks], [f"{t}d" for t in ticks]
    # 3 h and 12 h sit 0.375 d apart, so labelling both crowds the axis origin;
    # the published panels label 3 h then move to whole days.
    ticks = [3 / 24] + [float(d) for d in range(2, int(upto) + 1, 2)]
    labels = ["3h"] + [f"{int(d)}d" for d in range(2, int(upto) + 1, 2)]
    return ticks, labels


def day_tick_marks(ax, days, label_every: int = 2, rotation: float = 0,
                   fmt="{:g}") -> None:
    """A tick at every sampled day, a label on only every `label_every`-th.

    Every timepoint gets a mark so a reader can see where samples actually
    fall, without the axis becoming a wall of overlapping text. Ticks are real
    tick marks rather than minor ticks, so they inherit the axis style and
    survive into the EPS.
    """
    d = np.unique(np.asarray(list(days), dtype=float))
    d = d[np.isfinite(d)]
    ax.set_xticks(d)
    ax.set_xticklabels([fmt.format(v) if k % label_every == 0 else ""
                        for k, v in enumerate(d)], rotation=rotation)


def sampling_slots(index: Iterable[float], group: Group) -> np.ndarray:
    """Round a dense-grid index back to the nearest sampling slot."""
    _check_group(group)
    return np.rint(np.asarray(list(index), dtype=float)).astype(int)


def match_nearest(
    target_index: float,
    available: Sequence[float],
    values: Sequence[float],
    tolerance: float = 0.6,
):
    """Value at the sampling index nearest `target_index`, or None.

    The Jacobian is evaluated on a dense grid whose last point is not an integer
    (18.9 for cohort 1, 17.9 for cohort 2), so pairing it with a per-slot
    quantity such as Hill diversity needs a nearest-neighbour match rather than
    an exact lookup. `tolerance` of 0.6 slots is the published choice: it lets
    18.9 pair with slot 19 while refusing anything further than half a slot.
    """
    av = np.asarray(available, dtype=float)
    if av.size == 0:
        return None
    i = int(np.argmin(np.abs(av - target_index)))
    return values[i] if abs(av[i] - target_index) < tolerance else None
