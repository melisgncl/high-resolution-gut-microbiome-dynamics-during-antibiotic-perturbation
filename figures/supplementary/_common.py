"""Shared imports and helpers for the supplementary figures."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from succession import (anchors, diversity, io, jacobian, stats, style)  # noqa: E402,F401
from succession.config import (COHORT_1, COHORT_2, COLONISED, CONTROLS,  # noqa: E402,F401
                               CLONE_COLORS, MOUSE_COLORS,
                               MOUSE_COLORS_TAB, WINDOW, WINDOW_PRIMARY)
from succession.timeaxis import day_ticks, match_nearest, to_days        # noqa: E402,F401


def time_axis(ax, mouse):
    upto = 16 if mouse in COHORT_1 else 15
    ticks, labels = day_ticks("colonised", upto)
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, rotation=45, ha="right")
