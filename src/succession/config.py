"""Cohorts, palettes and the constants that decide published numbers."""

from __future__ import annotations

from pathlib import Path

# ── paths ────────────────────────────────────────────────────────────────────
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
OUT = ROOT / "figures" / "_out"
ASSETS = ROOT / "assets"

# ── cohorts ──────────────────────────────────────────────────────────────────
COHORT_1 = ["m1", "m2", "m3", "m4"]          # slots 1-19, last sample day 16
COHORT_2 = ["m5", "m6", "m7", "m8"]          # slots 1-18 (m5) / 2-18, day 15
COLONISED = COHORT_1 + COHORT_2
CONTROLS = ["c_m1", "c_m2", "c_m3", "c_m4"]  # daily, day 1-10

def group_of(mouse: str) -> str:
    """'colonised' or 'control' — the argument every timeaxis call needs."""
    return "control" if mouse.startswith("c_") else "colonised"


# ── estimator defaults ───────────────────────────────────────────────────────
# Each of these changes a published number. They are arguments everywhere in
# succession.jacobian, not hard-coded, so a reader can see what they do.

WINDOW = 5
"""Sliding-window width in sampling slots. The published main figures use 5."""

PSEUDOCOUNT = 1e-4
"""Added to relative abundance before log10 when building the state matrix.

Not cosmetic: the 16S tables are mostly zeros, so this sets how much weight a
near-absent family carries. Moving it to 1e-6 shifts per-mouse Figure 2C rho
for m3 from 0.74 to 0.78 and for m6 from 0.63 to 0.70.
"""

WARMUP_TOLERANCE = 2.0
"""Slack in the evaluation rule: keep t only if ``t - window >= grid_start - tol``.

An evaluation time owns the window ``(t - window, t]``. Near the start of a
series that window runs off the beginning of the data, making it an expanding
window rather than a sliding one. This rule discards those.

It is what sets the sample size, and the sample size is not incidental:

    window 1-3 -> n = 129      window 5 -> n = 113      window 10 -> n = 73

each extra step of window costing exactly one evaluation per mouse. Turning the
rule off entirely takes pooled Figure 2C rho from 0.73 to 0.60 and the number of
individually significant mice from 7 of 8 to 5 of 8.

Note the tolerance is 2.0, not 0: two of the fifteen windows kept per cohort-1
mouse are still partially clipped. Requiring a genuinely full window instead
gives rho = 0.752 over n = 97, still 7 of 8.
"""

MIN_POINTS = 4
"""Minimum dense-grid points inside a window for it to be evaluated."""

PSEUDOCOUNT_16S_PLOT = 1e-6
"""Added before the log in the Figure 1G-H and 4A composition panels.

Different purpose from PSEUDOCOUNT: here it only places zeros on the axis floor
so a vanishing family draws a spike to 1e-6 rather than leaving a gap. Matches
the lower axis limit.
"""

FAMILY_GATE = 1e-3
"""Per-mouse mean relative abundance required to draw a family in a panel.

Applied per mouse, not per cohort: Akkermansiaceae averages 2e-5 in m1 but
7.9e-3 in m5, so it is floor noise in one cohort and real signal in the other.
"""

# ── palettes ─────────────────────────────────────────────────────────────────
MOUSE_COLORS = {
    "m1": "#C0392B", "m2": "#E67E22", "m3": "#F39C12", "m4": "#D35400",
    "m5": "#8E44AD", "m6": "#2980B9", "m7": "#1ABC9C", "m8": "#2C3E50",
}

MOUSE_COLORS_TAB = {
    "m1": "#1F77B4", "m2": "#17BECF", "m3": "#2CA02C", "m4": "#9467BD",
    "m5": "#D62728", "m6": "#FF7F0E", "m7": "#8C564B", "m8": "#E377C2",
}
"""Used by the eigenvalue scatter, where eight lines share one axis."""

CONTROL_COLORS = {
    "c_m1": "#D62728", "c_m2": "#FF7F0E", "c_m3": "#9467BD", "c_m4": "#2CA02C",
}

LEGEND_FAMILIES = [
    "Enterobacteriaceae", "Paenibacillaceae",
    "Acholeplasmataceae", "Akkermansiaceae", "Bacteroidaceae",
    "Erysipelotrichaceae", "Lactobacillaceae", "Marinifilaceae",
    "Muribaculaceae", "Oscillospiraceae", "Prevotellaceae", "Ruminococcaceae",
]
"""The twelve families named in the published Figure 1 legend.

No abundance threshold reproduces this set: Lachnospiraceae peaks at 22.2% and
is excluded while Marinifilaceae peaks at 2.0% and is included. Applying the
per-mouse FAMILY_GATE to all 49 detected families yields these twelve plus
Lachnospiraceae and Rikenellaceae, so the published set is curated.
"""

FAMILY_COLORS = {
    "Enterobacteriaceae": "#C0392B", "Paenibacillaceae": "#1F5FA8",
    "Acholeplasmataceae": "#E75480", "Akkermansiaceae": "#1ABC9C",
    "Bacteroidaceae": "#8E2D5B", "Erysipelotrichaceae": "#C2189B",
    "Lactobacillaceae": "#E67E22", "Marinifilaceae": "#F5B041",
    "Muribaculaceae": "#7B241C", "Oscillospiraceae": "#E59866",
    "Prevotellaceae": "#D35400", "Ruminococcaceae": "#2471A3",
}

CLONE_COLORS = {
    "C1": "#2CA02C", "C2": "#1F77B4", "C3": "#D62728", "C4": "#FFE0A3",
    "C5": "#9467BD", "C6": "#FF7F0E", "C7": "#17A2A2", "C8": "#7FDBDA",
    "C9": "#E377C2", "C10": "#8C6D31", "C11": "#000000", "C12": "#1F3B73",
    "C13": "#B22222", "C14": "#9E9E9E",
}

# ── data quirks that must be applied at read time ────────────────────────────
CLONE_RELABEL = {"m2": {1: 2, 2: 1}, "m8": {1: 2, 2: 1}}
"""Clusters are numbered by mean frequency at the end of the experiment.

The clustering algorithm's own numbering does not always follow that, and in m2
and m8 C1 and C2 are the wrong way round (m2 end frequencies 7.03e-4 vs 1.22e-3;
m8 2.16e-4 vs 2.32e-3). m4 also ends with C2 above C1 but its labelling is
correct, so the end-frequency test is a flag rather than the rule.
"""

BARCODE_3H_EMPTY = ["m1", "m2", "m4", "m5"]
"""Mice whose 3 h barcode sample is present but essentially empty.

Hill 1D of 5.7, 13.0, 8.0 and 10.5 against 10568, 18984, 19330 and 30774 six
hours later - a handful of reads, not low diversity. m3 is the one cohort-1
mouse with a usable 3 h barcode sample (18172); m6-m8 have none at all.

Excluded from the Figure 2A diversity line, where they would otherwise plot at
0.76-1.11 on a log10 axis whose real signal starts near 4. They are *kept* in
the barcode LOESS that feeds the Jacobian, because removing them would start
cohort 1's dense grid a slot later and change every downstream evaluation count.
"""

N_COLOUR_BARCODES = 1000
"""Barcodes coloured in Figure 1C-D, ranked by peak frequency, per the caption."""
