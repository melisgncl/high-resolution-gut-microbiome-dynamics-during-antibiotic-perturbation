"""Reproduction code for the gut-colonisation succession manuscript.

Layout
------
`config`    cohorts, palettes, and every constant that changes a published number
`timeaxis`  index -> days, with the group required (the two clocks differ)
`io`        loaders for the shipped tables
`diversity` Hill numbers
`jacobian`  sliding-window interaction estimator and eigenvalues
`stats`     Spearman, Mann-Whitney, autocorrelation-aware permutation
`style`     shared matplotlib theme

Figures live in `figures/` and contain no analysis; this package contains no
plotting beyond `style`.
"""

from . import config, diversity, io, jacobian, stats, style, timeaxis  # noqa: F401

__version__ = "1.0.0"
__all__ = ["config", "diversity", "io", "jacobian", "stats", "style", "timeaxis"]
