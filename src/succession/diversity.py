"""Hill numbers.

Two Hill q1 series exist for the same mice and they are not identical:

* `io.load_diversity_16s` is the pipeline's own table, computed over the full
  family table, and for cohort 1 it has entries at only 15 of the 19 slots.
* `hill_q1_from_taxa` recomputes it from the curated taxon set that feeds the
  Jacobian, which exists at every slot.

Figure 2C uses the second. The difference is small (mean 0.010, max 0.032 for
m1) but present at every timepoint, and it is one reason a reproduction of the
published per-mouse coefficients needs the right source.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import io


def hill_q1(abundances) -> float:
    """exp(Shannon entropy). Zeros are dropped, then values renormalised."""
    v = np.asarray(abundances, dtype=float)
    v = v[v > 0]
    if v.size == 0:
        return np.nan
    p = v / v.sum()
    return float(np.exp(-np.sum(p * np.log(p))))


def hill_q0(abundances) -> float:
    """Richness."""
    v = np.asarray(abundances, dtype=float)
    return float((v > 0).sum())


def hill_q1_from_taxa(mouse: str) -> pd.DataFrame:
    """Hill q1 per sampling slot, from the curated Jacobian taxon set."""
    idx, _, mat = io.load_taxa_matrix(mouse)
    return pd.DataFrame({"index": idx,
                         "q1": [hill_q1(row) for row in mat]})
