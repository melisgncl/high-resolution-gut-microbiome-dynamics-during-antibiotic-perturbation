"""Shared matplotlib styling, so every panel matches without repetition."""

from __future__ import annotations

import os

import matplotlib as mpl
import matplotlib.pyplot as plt

from .config import OUT

#: Fixed timestamp for the `%%CreationDate` comment matplotlib writes into every
#: EPS. Without it the line carries the wall clock, so re-running `make figures`
#: rewrites all 24 EPS files with identical content and a new date, and
#: `figures/_out/` — which is committed — comes back dirty on every run for no
#: reason. See https://reproducible-builds.org/specs/source-date-epoch/.
#: Overridable: an environment variable that is already set wins.
SOURCE_DATE_EPOCH = "1735689600"        # 2025-01-01T00:00:00Z


def apply() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 9,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.titlesize": 10,
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.frameon": False,
        "legend.fontsize": 8,
        "figure.dpi": 110,
        "savefig.dpi": 200,
        "savefig.bbox": "tight",
        "savefig.facecolor": "white",
        # vector output that Illustrator can edit
        "ps.fonttype": 42,
        "pdf.fonttype": 42,
    })


def save(fig, name: str, formats=("png", "eps")):
    """Write a panel into figures/_out/ and hand the figure back.

    Both a PNG (for reading) and an EPS (for Illustrator, and for sending to a
    journal) are written by default. Note that EPS has no alpha channel, so
    semi-transparent fills — the s.d. bands, the ridgeline fills — are
    flattened against white rather than blended. The PNG is the faithful one;
    the EPS is the editable one.

    The figure is closed so a full run does not accumulate hundreds of open
    figures, but the object is returned intact. A closed figure still renders
    when it is displayed, which is how `notebooks/reproduce_figures.ipynb`
    shows each panel: it displays the object this function hands back, rather
    than reading the PNG below back off disk.
    """
    os.environ.setdefault("SOURCE_DATE_EPOCH", SOURCE_DATE_EPOCH)
    OUT.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        path = OUT / f"{name}.{fmt}"
        fig.savefig(path, format=fmt)
        print(f"  wrote {path.relative_to(OUT.parents[1])}")
    plt.close(fig)
    return fig


def stability_shading(ax, orientation: str = "vertical") -> None:
    """Green 'stable' region and the Re = 0 boundary used in Figure 3."""
    if orientation == "vertical":
        ax.axvspan(*(ax.get_xlim()[0], 0), color="#e8f4e8", alpha=0.45, zorder=0)
        ax.axvline(0, color="#cc2222", lw=0.9, ls="--", zorder=3)
    else:
        ax.axhspan(ax.get_ylim()[0], 0, color="#e8f4e8", alpha=0.45, zorder=0)
        ax.axhline(0, color="#cc2222", lw=0.9, ls="--", zorder=3)
