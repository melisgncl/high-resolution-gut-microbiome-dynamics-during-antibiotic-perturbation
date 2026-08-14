"""Phase anchors: when does the community stop collapsing and start recovering?

Two definitions appear in this project and they are not interchangeable.

`paeni_onset` is the one the manuscript uses: the first sampling slot at which
Paenibacillaceae reaches 10% of that mouse's own maximum relative abundance. It
marks the point the resident community begins to re-establish.

`diversity_nadir` is argmin of 16S Hill q1. It is reported here for the figures
that draw both, but it is unstable: q1 sits on a plateau near 1.0 for six to
eight consecutive slots during the collapse, so the exact argmin is noise.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import io
from .config import COLONISED, group_of
from .timeaxis import to_days


def paeni_onset(mouse: str, fraction: float = 0.10) -> float:
    fam = io.load_family(mouse)
    p = fam[fam["Family"] == "Paenibacillaceae"].sort_values("Time")
    if p.empty or p["abundance"].max() <= 0:
        return np.nan
    hit = p[p["abundance"] >= fraction * p["abundance"].max()]
    return float(hit["Time"].iloc[0]) if len(hit) else np.nan


def diversity_nadir(mouse: str) -> float:
    div = io.load_diversity_16s()
    d = div[(div["mouse"] == mouse) & div["q1"].notna()]
    return float(d.loc[d["q1"].idxmin(), "index"]) if len(d) else np.nan


def table(mice=None) -> pd.DataFrame:
    """Per-mouse anchors in sampling index and in days."""
    mice = mice or COLONISED
    rows = []
    for m in mice:
        g = group_of(m)
        on, na = paeni_onset(m), diversity_nadir(m)
        rows.append({"mouse": m, "onset_index": on, "nadir_index": na,
                     "onset_day": to_days(on, group=g) if np.isfinite(on) else np.nan,
                     "nadir_day": to_days(na, group=g) if np.isfinite(na) else np.nan})
    return pd.DataFrame(rows)
