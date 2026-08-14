"""
Figure 1c: rpsE MSA — FULL GENE — all three isolates + reference species
Author: Melis Gencel
Date: 2026-06-18
Input:  results/genomics/comparative/rpsE_msa_output.faa  (169 aligned columns)
Output: results/figures/genomics/fig1_rpsE_msa_full_gene.png + .pdf

Shows the complete rpsE alignment (M5/M6/M8 = 162 aa, I6 = 165 aa, 169 MAFFT columns)
in two horizontal blocks to confirm Δ(A21-K22-V23) is the ONLY difference in the
full gene. Mismatches outside the deletion = 0 (verified in rpsE_alignment_M5_vs_I6.txt).
rpsE is single-copy in all genomes (163x uniform coverage, no tandem duplication).
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
total_cols = len(m5_aligned)   # should be 169

# Row order: isolates → I6 → outgroups
ROWS = [
    ('M5', 'P. macerans M5',  m5_aligned,                          'isolate'),
    ('M6', 'P. macerans M6',  m5_aligned,                          'isolate'),
    ('M8', 'P. macerans M8',  m5_aligned,                          'isolate'),
    ('I6', 'P. macerans I6',  sequences['Paeni_I6_reference'],     'reference'),
    ('Ec', 'E. coli K-12',    sequences['Ecoli_K12_NP417762'],     'outgroup'),
    ('Bs', 'B. subtilis 168', sequences['Bsubtilis_168'],          'outgroup'),
]
labels = [r[1] for r in ROWS]
seqs   = [r[2] for r in ROWS]
rtypes = [r[3] for r in ROWS]
n_seqs = len(ROWS)
I6_ROW = 3

# Conservation per column
conservation = []
for ci in range(total_cols):
    col = [s[ci] if ci < len(s) else '-' for s in seqs]
    non_gap = [c for c in col if c != '-']
    if not non_gap:
        conservation.append(0.0)
    else:
        conservation.append(Counter(non_gap).most_common(1)[0][1] / n_seqs)

# Deletion columns (gap in M5, present in I6)
gap_cols = set(
    ci for ci in range(total_cols)
    if (ci < len(seqs[0]) and seqs[0][ci] == '-') and
       (ci < len(seqs[I6_ROW]) and seqs[I6_ROW][ci] != '-')
)

# ---- layout: two blocks stacked ----
# Block 1: cols 0..84  (85 columns — covers deletion region)
# Block 2: cols 85..168 (84 columns — rest of gene, all conserved)
BLOCK1_END = 85
blocks = [
    (0,          BLOCK1_END,  'Block 1 / 2  (N-terminal, residues 1–82 of I6)'),
    (BLOCK1_END, total_cols,  'Block 2 / 2  (C-terminal, residues 83–165 of I6)'),
]

CELL_W  = 0.72
CELL_H  = 0.85
LABEL_X = -1.0        # right edge of species labels (in cell units)
BLOCK_GAP = 2.5       # vertical gap between blocks (cell units)
TYPE_COLOR = {'isolate': '#B22222', 'reference': '#003399', 'outgroup': '#444444'}

max_block_cols = max(b[1] - b[0] for b in blocks)
fig_w = max_block_cols * CELL_W + 14   # extra for labels
fig_h = (n_seqs * CELL_H + BLOCK_GAP) * len(blocks) + 3

fig, ax = plt.subplots(figsize=(fig_w, fig_h))

def draw_block(ax, block_start, block_end, y_offset, block_label):
    cols = block_end - block_start

    # I6 residue counter up to this block
    i6_before = sum(
        1 for ci in range(block_start)
        if ci < len(seqs[I6_ROW]) and seqs[I6_ROW][ci] != '-'
    )

    # Block label on far left
    ax.text(LABEL_X * CELL_W - 0.3, y_offset + n_seqs * CELL_H / 2,
            block_label, ha='right', va='center', fontsize=7.5,
            color='#777777', rotation=90)

    for ri, (label, seq, rtype) in enumerate(zip(labels, seqs, rtypes)):
        y = y_offset + (n_seqs - 1 - ri) * CELL_H
        for ci_rel in range(cols):
            ci = block_start + ci_rel
            aa = seq[ci] if ci < len(seq) else '-'
            x = ci_rel * CELL_W
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

            # deletion columns: light red shading
            if ci in gap_cols:
                bg_fill = '#FFEAEA' if aa != '-' else '#FFD0D0'
            else:
                bg_fill = None

            rect = patches.FancyBboxPatch(
                (x + 0.03, y + 0.07), CELL_W * 0.92, CELL_H * 0.84,
                boxstyle="round,pad=0.015",
                facecolor=bg_fill if bg_fill and ci in gap_cols else bg,
                edgecolor='white', linewidth=0.4, zorder=2)
            ax.add_patch(rect)

            fw = 'bold'  if cons >= 1.0 else 'normal'
            fc = 'white' if cons >= 1.0 else ('#BBBBBB' if aa == '-' else 'black')
            if ci in gap_cols:
                fc = '#CC0000' if aa != '-' else '#CC000066'
            ax.text(x + CELL_W * 0.5, y + CELL_H * 0.5, aa,
                    ha='center', va='center', fontsize=6.5,
                    fontweight=fw, color=fc, fontfamily='monospace', zorder=3)

        # species label (only first block to avoid repetition)
        label_y = y + CELL_H * 0.5
        ax.text(LABEL_X * CELL_W, label_y, label,
                ha='right', va='center', fontsize=8.5,
                fontstyle='italic', fontweight='bold', color=TYPE_COLOR[rtype])

    # red box around deletion columns in this block
    del_in_block = sorted(ci for ci in gap_cols if block_start <= ci < block_end)
    if del_in_block:
        xs = (del_in_block[0] - block_start) * CELL_W - 0.05
        xw = (del_in_block[-1] - del_in_block[0] + 1) * CELL_W + 0.1
        yb = y_offset - 0.1
        yh = n_seqs * CELL_H + 0.2
        ax.add_patch(patches.FancyBboxPatch(
            (xs, yb), xw, yh,
            boxstyle="round,pad=0.04",
            facecolor='none', edgecolor='#CC0000',
            linewidth=2.2, zorder=4))
        mid_x = xs + xw / 2
        ax.annotate(
            'Δ(A21-K22-V23)\n3-aa deletion\n(only difference\nin full gene)',
            xy=(mid_x, yb),
            xytext=(mid_x + 12, yb - 1.6),
            fontsize=8, color='#CC0000', fontweight='bold',
            ha='center', va='top',
            arrowprops=dict(arrowstyle='->', color='#CC0000', lw=1.5),
            zorder=5)

    # position numbers every 10 I6 residues
    i6_pos = i6_before
    for ci_rel in range(cols):
        ci = block_start + ci_rel
        i6_res = seqs[I6_ROW][ci] if ci < len(seqs[I6_ROW]) else '-'
        if i6_res != '-':
            i6_pos += 1
            if i6_pos % 10 == 0:
                ax.text(ci_rel * CELL_W + CELL_W * 0.5,
                        y_offset + n_seqs * CELL_H + 0.15,
                        str(i6_pos), ha='center', va='bottom',
                        fontsize=6.5, color='#777777')

    # conservation bar
    bar_y = y_offset - 0.85
    for ci_rel in range(cols):
        ci = block_start + ci_rel
        cons = conservation[ci]
        h = cons * 0.55
        col = '#2B6CB0' if cons >= 1.0 else '#A8C8E8' if cons >= 0.75 else '#D0D0D0'
        ax.bar(ci_rel * CELL_W + CELL_W * 0.5, h, width=CELL_W * 0.75,
               bottom=bar_y, color=col, edgecolor='none')
    if block_start == 0:
        ax.text(LABEL_X * CELL_W, bar_y + 0.27, 'Conservation',
                ha='right', va='center', fontsize=7, color='#777777')

# draw both blocks
y_offsets = []
running_y = 0
for b_idx, (bstart, bend, blabel) in enumerate(blocks):
    draw_block(ax, bstart, bend, running_y, blabel)
    y_offsets.append(running_y)
    running_y += n_seqs * CELL_H + BLOCK_GAP

# separator line between blocks
sep_y = y_offsets[1] - BLOCK_GAP * 0.5
ax.axhline(sep_y, color='#CCCCCC', linewidth=0.8, linestyle='--',
           xmin=0.02, xmax=0.98)

# legend
leg_y = running_y - 0.2
leg_items = [
    ('#2B6CB0', 'Identical (all 6)'),
    ('#A8C8E8', 'Conserved (≥4/6)'),
    ('#E8E8E8', 'Variable'),
    ('#F4F4F4', 'Gap / deletion'),
]
for i, (col, txt) in enumerate(leg_items):
    xp = 3 + i * 16.0
    ax.add_patch(patches.Rectangle((xp, leg_y), 1.1, 0.7,
                                   facecolor=col, edgecolor='#999999', linewidth=0.5))
    ax.text(xp + 1.5, leg_y + 0.35, txt, fontsize=8, va='center')

# axes
max_cols = max(b[1] - b[0] for b in blocks)
ax.set_xlim(LABEL_X * CELL_W - 5, max_cols * CELL_W + 1)
ax.set_ylim(-3.5, running_y + 1.2)
ax.invert_yaxis()
ax.axis('off')

ax.set_title(
    'Ribosomal protein S5 (rpsE) — COMPLETE GENE ALIGNMENT (169 columns)\n'
    'Δ(A21-K22-V23) is the sole difference between isolates and I6 across the full 162/165 aa; '
    'rpsE is single-copy in all genomes (163× uniform coverage)',
    fontsize=10.5, fontweight='bold', pad=12)

plt.tight_layout()
out_png = os.path.join(OUT_DIR, 'fig1_rpsE_msa_full_gene.png')
out_pdf = os.path.join(OUT_DIR, 'fig1_rpsE_msa_full_gene.pdf')
plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(out_pdf, bbox_inches='tight', facecolor='white')
print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")
print(f"Total alignment columns: {total_cols}")
print(f"Deletion at columns: {sorted(gap_cols)}")
