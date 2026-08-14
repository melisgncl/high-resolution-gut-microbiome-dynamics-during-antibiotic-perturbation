"""
Figure 3: Linear genome comparison M5 vs I6
Author: Melis Gencel
Date: 2026-03-23
Input: results/genomics/comparative/M5_novel_regions.bed
       data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.tsv
Output: results/figures/genomics/fig3_genome_comparison.png + .pdf
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import os

BASE = r"C:\Users\melis\Desktop\bioinformatics-portfolio\hgt-study"
BED_FILE = os.path.join(BASE, "results", "genomics", "comparative", "M5_novel_regions.bed")
M5_TSV = os.path.join(BASE, "data", "sequence_data", "TPHP3P_results",
                       "TPHP3P_1_M5__+_", "annotation", "TPHP3P_1_M5__+_.tsv")
OUT_DIR = os.path.join(BASE, "results", "figures", "genomics")

M5_SIZE = 7_468_786
I6_SIZE = 7_113_008

# Parse BED
gap_regions = []
dup_regions = []
with open(BED_FILE) as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split('\t')
        start, stop = int(parts[1]), int(parts[2])
        rtype = parts[3]
        size = int(parts[4])
        entry = {'start': start, 'stop': stop, 'type': rtype, 'size': size}
        if rtype == 'GAP':
            gap_regions.append(entry)
        else:
            dup_regions.append(entry)

# Parse IS element positions from Bakta TSV
is_positions = []
is_keywords = ['transposase', 'IS110', 'IS4', 'IS3', 'IS200', 'IS605']
with open(M5_TSV) as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 8 and parts[1] == 'cds':
            if any(kw.lower() in parts[7].lower() for kw in is_keywords):
                is_positions.append((int(parts[2]) + int(parts[3])) // 2)

# Classify major GAP regions
major_regions = {
    'ICE #1 (108 kb)': (1_285_651, 1_393_969),
    'SP-beta prophage\n(59 kb)': (6_766_334, 6_825_514),
    'Integrative\nelement (48 kb)': (5_322_917, 5_371_421),
    'ICE #2 (25 kb)': (2_911_705, 2_936_697),
    'Low-GC block\n(22 kb)': (5_297_204, 5_322_425),
}

# Colour map for region types
region_colors = {
    'ICE': '#CC3333',
    'prophage': '#7744AA',
    'integrative': '#DD8833',
    'low-GC': '#33AA77',
    'other_gap': '#FFAA66',
    'dup': '#6699CC',
}

# Create figure
fig, ax = plt.subplots(figsize=(16, 5))

# Genome tracks
track_h = 0.6
m5_y = 3.0
i6_y = 1.0
scale = 1e6  # positions in Mb

# Draw M5 genome bar
ax.broken_barh([(0, M5_SIZE / scale)], (m5_y - track_h / 2, track_h),
               facecolor='#E8E8E8', edgecolor='#666666', linewidth=1)
ax.text(-0.15, m5_y, 'M5', fontsize=12, fontweight='bold', ha='right', va='center')
ax.text(-0.15, m5_y - 0.25, f'{M5_SIZE/1e6:.2f} Mb', fontsize=8, ha='right',
        va='center', color='grey')

# Draw I6 genome bar
ax.broken_barh([(0, I6_SIZE / scale)], (i6_y - track_h / 2, track_h),
               facecolor='#E8E8E8', edgecolor='#666666', linewidth=1)
ax.text(-0.15, i6_y, 'I6', fontsize=12, fontweight='bold', ha='right', va='center')
ax.text(-0.15, i6_y - 0.25, f'{I6_SIZE/1e6:.2f} Mb', fontsize=8, ha='right',
        va='center', color='grey')

# Draw GAP regions on M5 track
for reg in gap_regions:
    x = reg['start'] / scale
    w = (reg['stop'] - reg['start']) / scale

    # Determine colour
    is_ice = False
    for name, (ms, me) in major_regions.items():
        if reg['start'] >= ms - 1000 and reg['stop'] <= me + 1000:
            if 'ICE' in name:
                color = region_colors['ICE']
            elif 'prophage' in name:
                color = region_colors['prophage']
            elif 'Integrative' in name:
                color = region_colors['integrative']
            elif 'Low-GC' in name:
                color = region_colors['low-GC']
            is_ice = True
            break
    if not is_ice:
        color = region_colors['other_gap']

    ax.broken_barh([(x, w)], (m5_y - track_h / 2, track_h),
                   facecolor=color, edgecolor='none', alpha=0.85)

# Draw DUP regions as thin ticks
for reg in dup_regions:
    x = (reg['start'] + reg['stop']) / 2 / scale
    ax.plot([x, x], [m5_y - track_h / 2 - 0.1, m5_y + track_h / 2 + 0.1],
            color=region_colors['dup'], linewidth=0.5, alpha=0.6)

# IS element ticks on M5 (above track)
for pos in is_positions:
    x = pos / scale
    ax.plot([x, x], [m5_y + track_h / 2, m5_y + track_h / 2 + 0.15],
            color='#FF6600', linewidth=0.3, alpha=0.5)

# Label major regions
label_offsets = {
    'ICE #1 (108 kb)': (0, 0.8),
    'SP-beta prophage\n(59 kb)': (0, 0.8),
    'Integrative\nelement (48 kb)': (0, 0.8),
    'ICE #2 (25 kb)': (0, 0.8),
    'Low-GC block\n(22 kb)': (-0.2, 0.8),
}

for name, (ms, me) in major_regions.items():
    mid = (ms + me) / 2 / scale
    off_x, off_y = label_offsets[name]
    ax.annotate(name, xy=(mid, m5_y + track_h / 2),
                xytext=(mid + off_x, m5_y + track_h / 2 + off_y),
                fontsize=7, ha='center', va='bottom',
                arrowprops=dict(arrowstyle='->', color='#333333', lw=0.8),
                fontweight='bold', color='#333333')

# Position markers (Mb)
for pos_mb in range(0, 8):
    ax.plot([pos_mb, pos_mb], [i6_y - track_h / 2 - 0.15, i6_y - track_h / 2],
            color='grey', linewidth=0.5)
    ax.text(pos_mb, i6_y - track_h / 2 - 0.25, f'{pos_mb}', fontsize=8,
            ha='center', va='top', color='grey')
ax.text(3.5, i6_y - track_h / 2 - 0.55, 'Position (Mb)', fontsize=9,
        ha='center', va='top', color='grey')

# Connecting lines between M5 and I6 (general synteny)
# Draw a few light grey connecting shapes to show collinearity
for x_mb in np.arange(0, min(M5_SIZE, I6_SIZE) / scale, 0.5):
    ax.plot([x_mb, x_mb], [i6_y + track_h / 2, m5_y - track_h / 2],
            color='#DDDDDD', linewidth=0.3, alpha=0.5)

# rpsE position marker
rpse_pos = 6_425_000  # approximate from Bakta NKFIDM_05734
ax.annotate('rpsE\nΔ(AKV)', xy=(rpse_pos / scale, m5_y - track_h / 2),
            xytext=(rpse_pos / scale, m5_y - track_h / 2 - 0.7),
            fontsize=7, ha='center', va='top', color='#CC0000', fontweight='bold',
            arrowprops=dict(arrowstyle='->', color='#CC0000', lw=1.0))

# Legend
legend_items = [
    (region_colors['ICE'], 'ICE'),
    (region_colors['prophage'], 'Prophage'),
    (region_colors['integrative'], 'Integrative element'),
    (region_colors['low-GC'], 'Low-GC foreign block'),
    (region_colors['other_gap'], 'Other novel region'),
    (region_colors['dup'], 'DUP (IS copies)'),
]
for i, (color, label) in enumerate(legend_items):
    x_leg = 5.5 + (i % 3) * 0.9
    y_leg = 0.0 - (i // 3) * 0.3
    ax.add_patch(patches.Rectangle((x_leg, y_leg), 0.08, 0.2,
                                    facecolor=color, edgecolor='none'))
    ax.text(x_leg + 0.12, y_leg + 0.1, label, fontsize=7, va='center')

ax.set_xlim(-0.6, M5_SIZE / scale + 0.3)
ax.set_ylim(-0.8, 5.5)
ax.axis('off')
ax.set_title('Genome comparison: P. macerans M5 (isolate) vs I6 (reference)',
             fontsize=12, fontweight='bold', pad=10)

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'fig3_genome_comparison.png'), dpi=300,
            bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(OUT_DIR, 'fig3_genome_comparison.pdf'),
            bbox_inches='tight', facecolor='white')
print("Saved fig3_genome_comparison.png + .pdf")
