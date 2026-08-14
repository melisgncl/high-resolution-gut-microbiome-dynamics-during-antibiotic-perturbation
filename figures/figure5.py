"""Figure 5 — rpsE deletion and the resistance genotype.

A  RpsE protein alignment, reference against the gut isolates
D  CARD/RGI resistance gene survey of the isolate genome

Panels B (the deletion mapped onto the 70S ribosome) and C (colony morphology)
are a structural render and a photograph. Neither has code; both ship as images
under assets/ if available.

Run:  python figures/figure5.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from succession import io, style          # noqa: E402

style.apply()

CLASS_COLORS = {
    "Glycopeptide": "#4878CF", "Lincosamide": "#D65F5F",
    "Fosfomycin": "#6ACC65", "Disinfecting agent": "#B47CC7",
}


# ── A ────────────────────────────────────────────────────────────────────────
def panel_a(wrap: int = 60):
    """Reference vs isolate RpsE, with the deleted tripeptide marked.

    The reference reads ...I-N-R-V-A-K-V-V-K-G-G... and the isolates
    ...I-N-R-V-V-K-G-G..., three residues shorter. Which three are "deleted" is
    ambiguous: removing V-A-K at 20-22 or A-K-V at 21-23 yields the identical
    isolate sequence, because the region is a short repeat. The caption's
    "delta-VAK at 20-22" and the genomics record's "delta(A21-K22-V23)" are the
    same event described from opposite ends of that ambiguity, not a
    contradiction.
    """
    seqs = io.load_rpse_alignment()
    (ref_name, ref), (iso_name, iso) = list(seqs.items())

    # Align by padding the isolate at the first divergence.
    i = next(k for k in range(len(iso)) if ref[k] != iso[k])
    gap = len(ref) - len(iso)
    iso_aln = iso[:i] + "-" * gap + iso[i:]

    n_rows = int(np.ceil(len(ref) / wrap))
    fig, axes = plt.subplots(n_rows, 1, figsize=(11, 1.15 * n_rows))
    axes = np.atleast_1d(axes)

    for r, ax in enumerate(axes):
        s, e = r * wrap, min((r + 1) * wrap, len(ref))
        for k in range(s, e):
            x = k - s
            same = ref[k] == iso_aln[k]
            for y, ch in ((1, ref[k]), (0, iso_aln[k])):
                if not same:
                    ax.add_patch(Rectangle((x - 0.5, y - 0.42), 1, 0.84,
                                           fc="#F5B7B1" if ch != "-" else "#D5D8DC",
                                           ec="none", zorder=1))
                ax.text(x, y, ch, ha="center", va="center", fontsize=6.2,
                        family="monospace", zorder=2,
                        color="#C0392B" if not same else "#222222")
        ax.set_xlim(-1, wrap)
        ax.set_ylim(-0.7, 1.7)
        ax.set_yticks([1, 0])
        ax.set_yticklabels(["reference", "isolates"], fontsize=7)
        ax.set_xticks(range(0, min(wrap, e - s), 10))
        ax.set_xticklabels([str(s + t + 1) for t in range(0, min(wrap, e - s), 10)],
                           fontsize=6.5)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.tick_params(length=0)

    fig.suptitle(f"A   RpsE alignment — {ref_name.split()[0]} ({len(ref)} aa) vs "
                 f"gut isolates ({len(iso)} aa)", y=1.0, fontsize=10)
    fig.tight_layout()
    return style.save(fig, "fig5A_rpsE_alignment")


# ── D ────────────────────────────────────────────────────────────────────────
def panel_d():
    """Antimicrobial resistance gene survey against CARD (RGI, Strict).

    Ported from the manuscript's `32_figure_card_rgi.py`: a count of hits per
    drug class on the left, the hits themselves as a table on the right, and an
    explicit "Spectinomycin - not detected" row.

    That last row is the point of the figure. The detected genes fall in the
    glycopeptide, lincosamide, fosfomycin and disinfecting-agent classes; no
    spectinomycin determinant is present anywhere in the genome, which is what
    rules out an acquired cassette and leaves the chromosomal rpsE deletion as
    the only mechanism.
    """
    hits = io.load_card_hits()
    present = set(hits["drug_class"])
    classes = [c for c in CLASS_COLORS if c in present]        # CLASS_COLORS order
    counts = [int((hits["drug_class"] == c).sum()) for c in classes]

    fig = plt.figure(figsize=(14, 9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.7], wspace=0.08)
    ax, axt = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    axt.axis("off")

    # ── A: hits per drug class, with the negative result as its own row ──
    rows = classes + ["Spectinomycin"]
    y = np.arange(len(rows))[::-1]
    ax.barh(y[:len(classes)], counts, height=0.5,
            color=[CLASS_COLORS.get(c, "#999999") for c in classes])
    for yy, n in zip(y[:len(classes)], counts):
        ax.text(n + 0.15, yy, str(n), va="center", fontsize=11)
    ax.text(0.15, y[-1], "not detected", va="center", fontsize=11,
            style="italic", color="#CC2222")

    ax.set_yticks(y)
    ax.set_yticklabels(rows, fontsize=11)
    ax.get_yticklabels()[-1].set_color("#CC2222")
    ax.get_yticklabels()[-1].set_style("italic")
    ax.set_xlim(0, max(counts) + 2)
    ax.set_ylim(-0.9, len(rows) - 0.3)      # keep the negative row off the axis
    ax.set_xlabel("Number of resistance genes", fontsize=11)
    ax.set_title("A", fontsize=13, fontweight="bold", loc="left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

    # ── B: the hits themselves ──
    tint = {k: f"{c}22" for k, c in CLASS_COLORS.items()}      # 13% alpha hex
    body = [[r["gene"], r["drug_class"], r["mechanism"], f"{r['identity_pct']:.0f}"]
            for _, r in hits.iterrows()]
    colours = [[tint.get(r["drug_class"], "#EEEEEE")] * 4 for _, r in hits.iterrows()]

    tab = axt.table(cellText=body, cellColours=colours,
                    colLabels=["Gene / ARO hit", "Drug class", "Mechanism",
                               "Identity (%)"],
                    colColours=["#333333"] * 4, cellLoc="left", loc="center")
    tab.auto_set_font_size(False)
    tab.set_fontsize(10)
    tab.scale(1, 1.55)
    for (row, _), cell in tab.get_celld().items():
        cell.set_edgecolor("white")
        if row == 0:
            cell.set_text_props(color="white", fontweight="bold")
    axt.set_title("B", fontsize=13, fontweight="bold", loc="left")

    fig.suptitle("Antimicrobial resistance gene survey of "
                 "$\\it{P.\\ macerans}$ I6 (CARD/RGI)",
                 fontsize=14, fontweight="bold", y=0.97)
    fig.text(0.02, 0.03,
             "WGS assembly: P. macerans I6 (HL3J8K); RGI v6 / CARD database; "
             "Strict cut-off. No spectinomycin resistance determinants were "
             "detected.", fontsize=10, color="#555555")
    return style.save(fig, "fig5D_card_rgi")


if __name__ == "__main__":
    print("Figure 5")
    panel_a()
    panel_d()
    print("  B (ribosome render) and C (colony photograph) have no code")
