"""
Figure 1: rpsE MAFFT MSA with deletion highlighted
Author: Melis Gencel
Date: 2026-03-23
Input: results/genomics/comparative/rpsE_msa_output.faa
Output: results/figures/genomics/fig1_rpsE_msa.png + .pdf
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os

BASE = r"C:\Users\melis\Desktop\bioinformatics-portfolio\hgt-study"
MSA_FILE = os.path.join(BASE, "results", "genomics", "comparative", "rpsE_msa_output.faa")
OUT_DIR = os.path.join(BASE, "results", "figures", "genomics")

# Parse aligned FASTA
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

# Clean names
name_map = {
    'Paeni_M5_isolate': 'P. macerans M5',
    'Paeni_I6_reference': 'P. macerans I6',
    'Ecoli_K12_NP417762': 'E. coli K-12',
    'Bsubtilis_168': 'B. subtilis 168'
}

names = list(sequences.keys())
seqs = [sequences[n] for n in names]
labels = [name_map.get(n, n) for n in names]

# Show first 40 columns (N-terminal region with deletion)
show_cols = 40
seqs_trimmed = [s[:show_cols] for s in seqs]

# Amino acid colour scheme (Taylor)
aa_colors = {
    'A': '#CCFF00', 'R': '#0000FF', 'N': '#CC00FF', 'D': '#FF0000',
    'C': '#FFFF00', 'Q': '#FF00CC', 'E': '#FF0066', 'G': '#FF9900',
    'H': '#0066FF', 'I': '#66FF00', 'L': '#33FF00', 'K': '#6600FF',
    'M': '#00FF00', 'F': '#00FF66', 'P': '#FFCC00', 'S': '#FF3300',
    'T': '#FF6600', 'W': '#00CCFF', 'Y': '#00FFCC', 'V': '#99FF00',
    '-': '#FFFFFF', '*': '#FFFFFF'
}

# Conservation score per column
n_seqs = len(seqs_trimmed)
conservation = []
for col_idx in range(show_cols):
    col = [s[col_idx] if col_idx < len(s) else '-' for s in seqs_trimmed]
    non_gap = [c for c in col if c != '-']
    if len(non_gap) == 0:
        conservation.append(0)
    else:
        from collections import Counter
        counts = Counter(non_gap)
        most_common = counts.most_common(1)[0][1]
        conservation.append(most_common / n_seqs)

# Find the deletion columns (gap in M5 only)
gap_cols = []
for col_idx in range(show_cols):
    m5_res = seqs_trimmed[0][col_idx] if col_idx < len(seqs_trimmed[0]) else ''
    i6_res = seqs_trimmed[1][col_idx] if col_idx < len(seqs_trimmed[1]) else ''
    if m5_res == '-' and i6_res != '-':
        gap_cols.append(col_idx)

# Create figure
fig, ax = plt.subplots(figsize=(14, 3.5))

cell_w = 1.0
cell_h = 1.0

for row_idx, (label, seq) in enumerate(zip(labels, seqs_trimmed)):
    y = (n_seqs - 1 - row_idx) * cell_h
    for col_idx, aa in enumerate(seq):
        x = col_idx * cell_w

        # Background colour
        if aa == '-':
            bg = '#F0F0F0'
        elif conservation[col_idx] == 1.0:
            bg = '#4A86C8'  # fully conserved = blue
        elif conservation[col_idx] >= 0.75:
            bg = '#A8C8E8'  # mostly conserved = light blue
        else:
            bg = '#E8E8E8'  # variable = grey

        rect = patches.FancyBboxPatch((x, y), cell_w * 0.95, cell_h * 0.85,
                                       boxstyle="round,pad=0.02",
                                       facecolor=bg, edgecolor='white', linewidth=0.5)
        ax.add_patch(rect)

        # Residue letter
        fontweight = 'bold' if conservation[col_idx] == 1.0 else 'normal'
        colour = 'white' if conservation[col_idx] == 1.0 else 'black'
        if aa == '-':
            colour = '#999999'
        ax.text(x + cell_w * 0.475, y + cell_h * 0.42, aa,
                ha='center', va='center', fontsize=8, fontweight=fontweight,
                color=colour, fontfamily='monospace')

# Highlight deletion columns with red box
if gap_cols:
    x_start = min(gap_cols) * cell_w - 0.05
    x_width = (max(gap_cols) - min(gap_cols) + 1) * cell_w + 0.1
    y_bottom = -0.15
    y_height = n_seqs * cell_h + 0.1

    rect = patches.FancyBboxPatch((x_start, y_bottom), x_width, y_height,
                                   boxstyle="round,pad=0.05",
                                   facecolor='none', edgecolor='#CC0000',
                                   linewidth=2.5, linestyle='-')
    ax.add_patch(rect)

    # Annotation arrow
    mid_x = x_start + x_width / 2
    ax.annotate('Δ(A21-K22-V23)\nConserved Lys absent\nin isolates',
                xy=(mid_x, y_bottom - 0.05),
                xytext=(mid_x + 5, -1.8),
                fontsize=8, color='#CC0000', fontweight='bold',
                ha='center', va='top',
                arrowprops=dict(arrowstyle='->', color='#CC0000', lw=1.5))

# Species labels
for row_idx, label in enumerate(labels):
    y = (n_seqs - 1 - row_idx) * cell_h + cell_h * 0.42
    style = 'italic'
    ax.text(-0.5, y, label, ha='right', va='center', fontsize=9,
            fontstyle=style, fontweight='bold')

# Position numbers (every 5)
for col_idx in range(0, show_cols, 5):
    ax.text(col_idx * cell_w + cell_w * 0.475, n_seqs * cell_h + 0.15,
            str(col_idx + 1), ha='center', va='bottom', fontsize=7, color='grey')

# Conservation bar at bottom
bar_y = -0.8
for col_idx in range(show_cols):
    x = col_idx * cell_w
    h = conservation[col_idx] * 0.5
    color = '#4A86C8' if conservation[col_idx] == 1.0 else '#A8C8E8' if conservation[col_idx] >= 0.75 else '#D0D0D0'
    ax.bar(x + cell_w * 0.475, h, width=cell_w * 0.8, bottom=bar_y,
           color=color, edgecolor='none')

ax.text(-0.5, bar_y + 0.25, 'Conservation', ha='right', va='center', fontsize=7, color='grey')

# Legend
legend_y = n_seqs * cell_h + 0.7
for i, (color, label_text) in enumerate([('#4A86C8', 'Identical (4/4)'),
                                          ('#A8C8E8', 'Conserved (3/4)'),
                                          ('#E8E8E8', 'Variable'),
                                          ('#F0F0F0', 'Gap')]):
    x_pos = show_cols * cell_w * 0.55 + i * 5
    rect = patches.Rectangle((x_pos, legend_y), 0.8, 0.5, facecolor=color, edgecolor='grey', linewidth=0.5)
    ax.add_patch(rect)
    ax.text(x_pos + 1.2, legend_y + 0.25, label_text, fontsize=7, va='center')

ax.set_xlim(-8, show_cols * cell_w + 1)
ax.set_ylim(-2.5, n_seqs * cell_h + 1.5)
ax.set_aspect('equal')
ax.axis('off')
ax.set_title('Ribosomal protein S5 (rpsE) — N-terminal spectinomycin binding region',
             fontsize=11, fontweight='bold', pad=15)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'fig1_rpsE_msa.png'), dpi=300, bbox_inches='tight',
            facecolor='white')
plt.savefig(os.path.join(OUT_DIR, 'fig1_rpsE_msa.pdf'), bbox_inches='tight',
            facecolor='white')
print("Saved fig1_rpsE_msa.png + .pdf")
