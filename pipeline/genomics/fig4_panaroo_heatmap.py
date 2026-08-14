"""
Figure 4: Panaroo gene presence/absence heatmap
Author: Melis Gencel
Date: 2026-03-23
Input: results/genomics/comparative/panaroo_output/gene_presence_absence.Rtab
       results/tables/genomics/panaroo_gene_classification.csv
Output: results/figures/genomics/fig4_panaroo_heatmap.png + .pdf
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import csv
import os

BASE = r"C:\Users\melis\Desktop\bioinformatics-portfolio\hgt-study"
RTAB = os.path.join(BASE, "results", "genomics", "comparative",
                     "panaroo_output", "gene_presence_absence.Rtab")
CLASS_CSV = os.path.join(BASE, "results", "tables", "genomics",
                          "panaroo_gene_classification.csv")
OUT_DIR = os.path.join(BASE, "results", "figures", "genomics")

# Read classification to get presence patterns
gene_patterns = {}
gene_funcs = {}
with open(CLASS_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        gene_patterns[row['Gene']] = row['Presence_pattern']
        gene_funcs[row['Gene']] = row['Functional_category']

# Read Rtab
header = None
data = []
gene_names = []
with open(RTAB) as f:
    for i, line in enumerate(f):
        parts = line.strip().split('\t')
        if i == 0:
            header = parts[1:]  # skip 'Gene' column
            continue
        gene_names.append(parts[0])
        data.append([int(x) for x in parts[1:]])

data = np.array(data)

# Sort by pattern: core first, then novel, then lost
order_map = {'core': 0, 'novel': 1, 'lost': 2, 'absent': 3}
sort_keys = [order_map.get(gene_patterns.get(g, 'absent'), 3) for g in gene_names]
sort_idx = np.argsort(sort_keys, kind='stable')
data_sorted = data[sort_idx]
gene_names_sorted = [gene_names[i] for i in sort_idx]

# Find boundaries
patterns_sorted = [gene_patterns.get(gene_names_sorted[i], 'absent') for i in range(len(gene_names_sorted))]
core_end = sum(1 for p in patterns_sorted if p == 'core')
novel_end = core_end + sum(1 for p in patterns_sorted if p == 'novel')
lost_end = novel_end + sum(1 for p in patterns_sorted if p == 'lost')

# Functional category colour strip for novel and lost genes
func_colors = {
    'IS_transposase': '#FF6600',
    'phage': '#7744AA',
    'conjugation_ICE': '#CC3333',
    'integrase_recombinase': '#DD8833',
    'CRISPR': '#33AACC',
    'toxin_antitoxin': '#CC6699',
    'drug_efflux': '#FF3366',
    'transport': '#6699CC',
    'regulatory': '#99CC33',
    'metabolism_enzyme': '#33AA77',
    'hypothetical_DUF': '#BBBBBB',
    'ribosomal_translation': '#9966CC',
    'sporulation': '#CC9933',
    'cell_envelope_motility': '#66CCAA',
    'other': '#DDDDDD',
}

# Create figure with two panels
fig, (ax_heat, ax_func) = plt.subplots(1, 2, figsize=(8, 10),
                                         gridspec_kw={'width_ratios': [4, 0.8],
                                                       'wspace': 0.02})

# Heatmap
cmap = plt.cm.colors.ListedColormap(['#FFFFFF', '#2166AC'])
ax_heat.imshow(data_sorted, aspect='auto', cmap=cmap, interpolation='none')

# Column labels
col_labels = ['M5', 'M6', 'M8', 'I6']
ax_heat.set_xticks(range(4))
ax_heat.set_xticklabels(col_labels, fontsize=10, fontweight='bold')
ax_heat.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)

# Category boundaries
for y_pos, label, count in [(core_end / 2, f'Core\n({core_end})', core_end),
                             ((core_end + novel_end) / 2, f'Novel\n({novel_end - core_end})', novel_end - core_end),
                             ((novel_end + lost_end) / 2, f'Lost\n({lost_end - novel_end})', lost_end - novel_end)]:
    ax_heat.text(-0.8, y_pos, label, fontsize=9, ha='right', va='center', fontweight='bold')

# Draw boundary lines
for y_pos in [core_end, novel_end]:
    ax_heat.axhline(y=y_pos - 0.5, color='#CC0000', linewidth=1.5, linestyle='-')

ax_heat.set_yticks([])
ax_heat.set_ylabel('')

# Functional category strip (only for novel + lost genes)
func_strip = np.zeros((len(gene_names_sorted), 1, 3))
func_strip[:, 0] = [1.0, 1.0, 1.0]  # white default

for i in range(len(gene_names_sorted)):
    gene = gene_names_sorted[i]
    func = gene_funcs.get(gene, 'other')
    color_hex = func_colors.get(func, '#DDDDDD')
    # Convert hex to RGB
    r = int(color_hex[1:3], 16) / 255
    g = int(color_hex[3:5], 16) / 255
    b = int(color_hex[5:7], 16) / 255
    if i >= core_end:  # only colour novel and lost
        func_strip[i, 0] = [r, g, b]

ax_func.imshow(func_strip, aspect='auto', interpolation='none')
ax_func.set_xticks([0])
ax_func.set_xticklabels(['Function'], fontsize=8)
ax_func.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
ax_func.set_yticks([])

# Draw boundary lines on func strip too
for y_pos in [core_end, novel_end]:
    ax_func.axhline(y=y_pos - 0.5, color='#CC0000', linewidth=1.5)

# Legend for functional categories (only show categories present in novel/lost)
novel_lost_funcs = set()
for i in range(core_end, len(gene_names_sorted)):
    gene = gene_names_sorted[i]
    novel_lost_funcs.add(gene_funcs.get(gene, 'other'))

legend_patches = []
for cat in sorted(novel_lost_funcs):
    color = func_colors.get(cat, '#DDDDDD')
    label = cat.replace('_', ' ')
    legend_patches.append(patches.Patch(facecolor=color, edgecolor='grey',
                                         linewidth=0.5, label=label))

ax_func.legend(handles=legend_patches, loc='center left', bbox_to_anchor=(1.5, 0.5),
               fontsize=7, frameon=True, title='Functional\ncategory', title_fontsize=8)

fig.suptitle('Pan-genome gene presence/absence\nP. macerans isolates vs I6 reference',
             fontsize=12, fontweight='bold', y=0.98)

plt.savefig(os.path.join(OUT_DIR, 'fig4_panaroo_heatmap.png'), dpi=300,
            bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(OUT_DIR, 'fig4_panaroo_heatmap.pdf'),
            bbox_inches='tight', facecolor='white')
print("Saved fig4_panaroo_heatmap.png + .pdf")
