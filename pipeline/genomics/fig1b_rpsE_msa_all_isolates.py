"""
Figure 1b: rpsE MSA — all three isolates (M5, M6, M8) + I6 reference + outgroups
Author: Melis Gencel
Date: 2026-06-18
Input:  results/genomics/comparative/rpsE_msa_output.faa  (M5/I6/Ecoli/Bsubt alignment)
Output: results/figures/genomics/fig1_rpsE_msa_all_isolates.png + .pdf

M6 and M8 share the identical rpsE sequence as M5 (confirmed: Panaroo pan-genome,
direct pairwise comparison, 99.26% FastANI). Their aligned positions are therefore
identical to M5 and are plotted as separate rows without re-running MAFFT.
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from collections import Counter
import os

BASE = r"C:\Users\melis\Desktop\bioinformatics-portfolio\hgt-study"
MSA_FILE = os.path.join(BASE, "results", "genomics", "comparative", "rpsE_msa_output.faa")
OUT_DIR = os.path.join(BASE, "results", "figures", "genomics")

# ---------- parse MAFFT output ----------
sequences = {}
current = None
with open(MSA_FILE) as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            current = line[1:]
            sequences[current] = ''
        elif current:
            sequences[current] += line

m5_aligned = sequences['Paeni_M5_isolate']

# Display order: isolates (M5/M6/M8) → reference (I6) → outgroups (Ecoli, Bsubt)
ROWS = [
    ('M5', 'P. macerans M5 (isolate)', m5_aligned, 'isolate'),
    ('M6', 'P. macerans M6 (isolate)', m5_aligned, 'isolate'),
    ('M8', 'P. macerans M8 (isolate)', m5_aligned, 'isolate'),
    ('I6', 'P. macerans I6 (reference)', sequences['Paeni_I6_reference'], 'reference'),
    ('Ec', 'E. coli K-12',             sequences['Ecoli_K12_NP417762'],   'outgroup'),
    ('Bs', 'B. subtilis 168',          sequences['Bsubtilis_168'],         'outgroup'),
]

labels  = [r[1] for r in ROWS]
seqs    = [r[2] for r in ROWS]
rtypes  = [r[3] for r in ROWS]
n_seqs  = len(ROWS)

SHOW = 45   # columns to display (covers the N-terminal deletion region)
seqs_t = [s[:SHOW] for s in seqs]

# ---------- conservation per column ----------
conservation = []
for ci in range(SHOW):
    col = [s[ci] if ci < len(s) else '-' for s in seqs_t]
    non_gap = [c for c in col if c != '-']
    if not non_gap:
        conservation.append(0.0)
    else:
        most_common = Counter(non_gap).most_common(1)[0][1]
        conservation.append(most_common / n_seqs)

# deletion columns = gap in isolate rows, present in I6
I6_ROW = 3
gap_cols = [
    ci for ci in range(SHOW)
    if (ci < len(seqs_t[0]) and seqs_t[0][ci] == '-') and
       (ci < len(seqs_t[I6_ROW]) and seqs_t[I6_ROW][ci] != '-')
]

# ---------- layout ----------
CELL_W = 1.0
CELL_H = 1.0
LABEL_X = -1.2     # right edge of species labels
BRACKET_X = -10.5  # left edge of grouping label

fig, ax = plt.subplots(figsize=(18, 6.5))

# ---------- cells ----------
for ri, (label, seq, rtype) in enumerate(zip(labels, seqs_t, rtypes)):
    y = (n_seqs - 1 - ri) * CELL_H
    for ci, aa in enumerate(seq):
        x = ci * CELL_W
        cons = conservation[ci]

        if aa == '-':
            bg = '#F4F4F4'
        elif cons >= 1.0:
            bg = '#2B6CB0'
        elif cons >= 5/6:
            bg = '#4A86C8'
        elif cons >= 4/6:
            bg = '#A8C8E8'
        else:
            bg = '#E8E8E8'

        rect = patches.FancyBboxPatch(
            (x + 0.03, y + 0.08), CELL_W * 0.92, CELL_H * 0.82,
            boxstyle="round,pad=0.02",
            facecolor=bg, edgecolor='white', linewidth=0.5)
        ax.add_patch(rect)

        fw     = 'bold'   if cons >= 1.0 else 'normal'
        fc     = 'white'  if cons >= 1.0 else ('#BBBBBB' if aa == '-' else 'black')
        ax.text(x + CELL_W * 0.5, y + CELL_H * 0.5, aa,
                ha='center', va='center', fontsize=8,
                fontweight=fw, color=fc, fontfamily='monospace')

# ---------- deletion highlight ----------
if gap_cols:
    xs = min(gap_cols) * CELL_W - 0.08
    xw = (max(gap_cols) - min(gap_cols) + 1) * CELL_W + 0.16
    yb = -0.12
    yh = n_seqs * CELL_H + 0.12

    # shaded column background
    ax.add_patch(patches.FancyBboxPatch(
        (xs, yb), xw, yh,
        boxstyle="round,pad=0.04",
        facecolor='#FFEAEA', edgecolor='#CC0000',
        linewidth=2.5, zorder=0))

    mid_x = xs + xw / 2
    ax.annotate(
        'Δ(A21-K22-V23)\nConserved Lys absent\nin all 3 isolates',
        xy=(mid_x, yb),
        xytext=(mid_x + 6, -2.5),
        fontsize=9, color='#CC0000', fontweight='bold',
        ha='center', va='top',
        arrowprops=dict(arrowstyle='->', color='#CC0000', lw=1.5))

# ---------- species labels ----------
TYPE_COLOR = {'isolate': '#B22222', 'reference': '#003399', 'outgroup': '#444444'}
for ri, (label, rtype) in enumerate(zip(labels, rtypes)):
    y = (n_seqs - 1 - ri) * CELL_H + CELL_H * 0.5
    ax.text(LABEL_X, y, label, ha='right', va='center', fontsize=9,
            fontstyle='italic', fontweight='bold', color=TYPE_COLOR[rtype])

# ---------- grouping brackets ----------
def draw_bracket(ax, x_tip, y_top, y_bot, color):
    """Vertical bracket: two horizontal stubs + connecting vertical line."""
    stub = 0.6
    ax.plot([x_tip, x_tip + stub], [y_top, y_top], color=color, lw=1.5, solid_capstyle='round')
    ax.plot([x_tip, x_tip + stub], [y_bot, y_bot], color=color, lw=1.5, solid_capstyle='round')
    ax.plot([x_tip, x_tip],        [y_top, y_bot], color=color, lw=1.5, solid_capstyle='round')

# Isolates bracket (rows 0-2)
iso_top = (n_seqs - 1) * CELL_H + CELL_H * 0.9
iso_bot = (n_seqs - 1 - 2) * CELL_H + CELL_H * 0.1
draw_bracket(ax, BRACKET_X + 1.0, iso_top, iso_bot, '#B22222')
ax.text(BRACKET_X + 0.6, (iso_top + iso_bot) / 2, 'Resistant\nisolates',
        ha='right', va='center', fontsize=8, color='#B22222', fontweight='bold',
        rotation=90, multialignment='center')

# Reference bracket (row 3)
ref_y = (n_seqs - 1 - 3) * CELL_H + CELL_H * 0.5
ax.text(BRACKET_X + 0.6, ref_y, 'Reference\n(sensitive?)',
        ha='right', va='center', fontsize=8, color='#003399', fontweight='bold',
        rotation=90, multialignment='center')

# Outgroups bracket (rows 4-5)
out_top = (n_seqs - 1 - 4) * CELL_H + CELL_H * 0.9
out_bot = (n_seqs - 1 - 5) * CELL_H + CELL_H * 0.1
draw_bracket(ax, BRACKET_X + 1.0, out_top, out_bot, '#444444')
ax.text(BRACKET_X + 0.6, (out_top + out_bot) / 2, 'Outgroups\n(sensitive)',
        ha='right', va='center', fontsize=8, color='#444444', fontweight='bold',
        rotation=90, multialignment='center')

# ---------- position numbers (relative to I6 residue count) ----------
i6_pos = 0
for ci in range(SHOW):
    i6_res = seqs_t[I6_ROW][ci] if ci < len(seqs_t[I6_ROW]) else '-'
    if i6_res != '-':
        i6_pos += 1
        if i6_pos % 5 == 0:
            ax.text(ci * CELL_W + CELL_W * 0.5, n_seqs * CELL_H + 0.15,
                    str(i6_pos), ha='center', va='bottom', fontsize=7, color='#666666')

ax.text(LABEL_X, n_seqs * CELL_H + 0.35, 'I6 pos.',
        ha='right', va='bottom', fontsize=7, color='#666666')

# ---------- conservation bar ----------
bar_y = -1.1
for ci in range(SHOW):
    x = ci * CELL_W
    h = conservation[ci] * 0.65
    col = '#2B6CB0' if conservation[ci] >= 1.0 else '#A8C8E8' if conservation[ci] >= 0.75 else '#D0D0D0'
    ax.bar(x + CELL_W * 0.5, h, width=CELL_W * 0.78, bottom=bar_y,
           color=col, edgecolor='none')

ax.text(LABEL_X, bar_y + 0.32, 'Conservation', ha='right', va='center',
        fontsize=7, color='#666666')

# ---------- legend ----------
leg_y = n_seqs * CELL_H + 0.65
leg_items = [
    ('#2B6CB0', 'Identical (all 6)'),
    ('#A8C8E8', 'Conserved (≥4/6)'),
    ('#E8E8E8', 'Variable'),
    ('#F4F4F4', 'Gap / deletion'),
]
for i, (col, txt) in enumerate(leg_items):
    xp = SHOW * CELL_W * 0.40 + i * 8.0
    ax.add_patch(patches.Rectangle((xp, leg_y), 0.85, 0.65,
                                   facecolor=col, edgecolor='#999999', linewidth=0.5))
    ax.text(xp + 1.15, leg_y + 0.32, txt, fontsize=7.5, va='center')

# ---------- axes ----------
ax.set_xlim(BRACKET_X - 0.5, SHOW * CELL_W + 2)
ax.set_ylim(-4.0, n_seqs * CELL_H + 2.2)
ax.set_aspect('equal')
ax.axis('off')

ax.set_title(
    'Ribosomal protein S5 (rpsE) — N-terminal spectinomycin-binding region\n'
    'All three resistant isolates (M5, M6, M8) carry the identical Δ(A21-K22-V23) deletion; '
    'conserved Lys present in I6, E. coli, and B. subtilis',
    fontsize=10.5, fontweight='bold', pad=14)

plt.tight_layout()
out_png = os.path.join(OUT_DIR, 'fig1_rpsE_msa_all_isolates.png')
out_pdf = os.path.join(OUT_DIR, 'fig1_rpsE_msa_all_isolates.pdf')
plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(out_pdf, bbox_inches='tight', facecolor='white')
print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")
