"""
Figure 5 (revised): rpsE read-level validation vs I6 — all three isolates
Author: Melis Gencel
Date: 2026-06-18
Input:  results/genomics/comparative/{M5,M6,M8}_rpsE_depth.tsv   (real samtools depth vs I6)
        results/genomics/comparative/rpsE_vs_I6_readlevel.log    (read counts)
Output: results/figures/genomics/fig5_rpsE_pileup_vs_I6.png + .pdf

Real numbers (reads mapped to I6 rpsE, contig_1:228319-228816):
  M5: mean 160x, 159/171 spanning reads carry 9-bp deletion (93%)
  M6: mean 133x, 133/139 (96%)
  M8: mean  95x,  95/101 (94%)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os

BASE = r"C:\Users\melis\Desktop\bioinformatics-portfolio\hgt-study"
COMP = os.path.join(BASE, "results", "genomics", "comparative")
OUT_DIR = os.path.join(BASE, "results", "figures", "genomics")

RPSE_START, RPSE_END = 228319, 228816          # rpsE on I6 (contig_1)
# AKV deletion: I6 residues 21-23 -> codons 61-69 bp into the CDS from start.
# rpsE is on the minus strand; the 9-bp deletion maps near the 5' (high-coord) end.
DEL_LEN_BP = 9

ISOLATES = [
    ("M5", "#2166AC", 159, 171, 160),
    ("M6", "#1B7837", 133, 139, 133),
    ("M8", "#B2182B",  95, 101,  95),
]

def load_depth(name):
    pos, dep = [], []
    with open(os.path.join(COMP, f"{name}_rpsE_depth.tsv")) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) == 3:
                pos.append(int(p[1])); dep.append(int(p[2]))
    return np.array(pos), np.array(dep)

fig, axes = plt.subplots(4, 1, figsize=(11, 9),
                         gridspec_kw={'height_ratios': [1, 1, 1, 1.15]})

# ---- Panels A-C: coverage for each isolate ----
for ax, (name, color, ndel, ntot, meancov) in zip(axes[:3], ISOLATES):
    pos, dep = load_depth(name)
    ax.fill_between(pos, dep, color=color, alpha=0.30, zorder=1)
    ax.plot(pos, dep, color=color, lw=0.8, zorder=2)
    # rpsE body
    ax.axvspan(RPSE_START, RPSE_END, color='#FF9933', alpha=0.10, zorder=0)
    ax.axvline(RPSE_START, color='#FF9933', lw=1, ls='--', alpha=0.7)
    ax.axvline(RPSE_END,   color='#FF9933', lw=1, ls='--', alpha=0.7)
    # stats box
    frac = 100.0 * ndel / ntot
    ax.text(0.015, 0.90,
            f"{name}   mean {meancov}×\n{ndel}/{ntot} spanning reads carry the\n9-bp deletion ({frac:.0f}%)",
            transform=ax.transAxes, va='top', ha='left', fontsize=8.5,
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=color, alpha=0.9))
    ax.set_ylabel("Read depth", fontsize=9)
    ax.set_ylim(0, max(dep) * 1.30)
    ax.set_xlim(pos.min(), pos.max())
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(labelbottom=False)

axes[0].set_title("Read coverage across the I6 rpsE locus (reads mapped to I6, contig_1:228,319–228,816)",
                  fontsize=10.5, fontweight='bold')
axes[0].text(RPSE_START + (RPSE_END - RPSE_START) / 2, axes[0].get_ylim()[1] * 0.92,
             "rpsE (498 bp)", ha='center', fontsize=8.5, fontweight='bold', color='#CC6600')

# ---- Panel D: schematic of the deletion + CIGAR ----
ax = axes[3]
ax.set_xlim(0, 10); ax.set_ylim(0, 4.2); ax.axis('off')
ax.set_title("rpsE deletion schematic — isolates vs I6 reference",
             fontsize=10.5, fontweight='bold')

# I6 bar (165 aa)
i6_y = 3.1
ax.add_patch(patches.FancyBboxPatch((0.6, i6_y - 0.18), 8.4, 0.36,
            boxstyle="round,pad=0.04", fc='#A8C8E8', ec='#333'))
ax.text(0.35, i6_y, "I6", fontsize=10, fontweight='bold', ha='right', va='center')
ax.text(4.8, i6_y, "rpsE — 165 aa (498 bp)", ha='center', va='center', fontsize=9)
akv_x = 0.6 + (20 / 165) * 8.4
akv_w = (3 / 165) * 8.4
ax.add_patch(patches.Rectangle((akv_x, i6_y - 0.20), akv_w, 0.40,
            fc='#33AA77', ec='#006633', lw=1.8))
ax.text(akv_x + akv_w / 2, i6_y + 0.42, "A-K-V", ha='center', fontsize=8.5,
        fontweight='bold', color='#006633')

# isolate bar (162 aa)
m_y = 1.7
ax.add_patch(patches.FancyBboxPatch((0.6, m_y - 0.18), 8.25, 0.36,
            boxstyle="round,pad=0.04", fc='#FFB366', ec='#333'))
ax.text(0.35, m_y, "M5\nM6\nM8", fontsize=8.5, fontweight='bold', ha='right', va='center',
        linespacing=1.0, color='#7a3b00')
ax.text(4.75, m_y, "rpsE — 162 aa (489 bp)", ha='center', va='center', fontsize=9)
ax.annotate("", xy=(akv_x, m_y + 0.20), xytext=(akv_x, i6_y - 0.20),
            arrowprops=dict(arrowstyle='<->', color='#CC0000', lw=1.6))
ax.text(akv_x + akv_w / 2 + 0.15, (m_y + i6_y) / 2, "Δ9 bp\n(3 aa)",
        ha='left', va='center', fontsize=8, color='#CC0000', fontweight='bold')

ax.text(5, 0.85, "minimap2 CIGAR signature:  9D261M  (9-bp deletion + 261-bp match)",
        ha='center', fontsize=9.5, fontweight='bold', color='#2166AC')
ax.text(5, 0.40, "Deletion carried by 93–96% of spanning reads in all three isolates",
        ha='center', fontsize=9, color='#555')

plt.tight_layout()
png = os.path.join(OUT_DIR, "fig5_rpsE_pileup_vs_I6.png")
pdf = os.path.join(OUT_DIR, "fig5_rpsE_pileup_vs_I6.pdf")
plt.savefig(png, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(pdf, bbox_inches='tight', facecolor='white')
print("Saved:", png)
print("Saved:", pdf)
