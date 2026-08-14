"""
Figure 6: Pan-genome functional breakdown
Author: Melis Gencel
Date: 2026-03-23
Input: results/tables/genomics/panaroo_gene_classification.csv
Output: results/figures/genomics/fig6_pangenome_functional.png + .pdf
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv
from collections import Counter, OrderedDict
import os

BASE = r"C:\Users\melis\Desktop\bioinformatics-portfolio\hgt-study"
CLASS_CSV = os.path.join(BASE, "results", "tables", "genomics",
                          "panaroo_gene_classification.csv")
OUT_DIR = os.path.join(BASE, "results", "figures", "genomics")

# Read classification
rows = []
with open(CLASS_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        rows.append(row)

# Colour map
func_colors = OrderedDict([
    ('metabolism_enzyme', '#33AA77'),
    ('transport', '#6699CC'),
    ('regulatory', '#99CC33'),
    ('hypothetical_DUF', '#BBBBBB'),
    ('IS_transposase', '#FF6600'),
    ('phage', '#7744AA'),
    ('conjugation_ICE', '#CC3333'),
    ('integrase_recombinase', '#DD8833'),
    ('CRISPR', '#33AACC'),
    ('drug_efflux', '#FF3366'),
    ('toxin_antitoxin', '#CC6699'),
    ('sporulation', '#CC9933'),
    ('ribosomal_translation', '#9966CC'),
    ('cell_envelope_motility', '#66CCAA'),
    ('other', '#DDDDDD'),
])

func_labels = {k: k.replace('_', ' ').title() for k in func_colors}
func_labels['IS_transposase'] = 'IS Transposase'
func_labels['CRISPR'] = 'CRISPR'
func_labels['hypothetical_DUF'] = 'Hypothetical / DUF'
func_labels['conjugation_ICE'] = 'Conjugation / ICE'
func_labels['drug_efflux'] = 'Drug Efflux'

# Count by pattern and function
patterns = ['core', 'novel', 'lost']
pattern_labels = ['Core\n(5,825)', 'Novel\n(736)', 'Lost\n(309)']

data = {}
for pattern in patterns:
    subset = [r for r in rows if r['Presence_pattern'] == pattern]
    func_counts = Counter(r['Functional_category'] for r in subset)
    data[pattern] = func_counts

# Create figure — 2 panels
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6),
                                gridspec_kw={'width_ratios': [1, 1.2]})

# Panel A: Stacked bar chart
x_pos = np.arange(len(patterns))
bar_width = 0.6
bottoms = np.zeros(len(patterns))

for cat, color in func_colors.items():
    values = [data[p].get(cat, 0) for p in patterns]
    ax1.bar(x_pos, values, bar_width, bottom=bottoms, color=color,
            edgecolor='white', linewidth=0.3, label=func_labels[cat])
    bottoms += values

ax1.set_xticks(x_pos)
ax1.set_xticklabels(pattern_labels, fontsize=10, fontweight='bold')
ax1.set_ylabel('Number of gene families', fontsize=11)
ax1.set_title('A. Gene families by presence pattern', fontsize=11, fontweight='bold')
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Add count labels on top
for i, pattern in enumerate(patterns):
    total = sum(data[pattern].values())
    ax1.text(i, total + 50, str(total), ha='center', va='bottom',
             fontsize=9, fontweight='bold')

# Panel B: Percentage breakdown for novel genes only
novel_counts = data['novel']
total_novel = sum(novel_counts.values())

# Sort by count descending
sorted_cats = sorted(novel_counts.keys(), key=lambda c: -novel_counts[c])
# Only show categories with >= 2 genes
sorted_cats = [c for c in sorted_cats if novel_counts[c] >= 2]

y_pos = np.arange(len(sorted_cats))
values = [novel_counts[c] for c in sorted_cats]
colors = [func_colors.get(c, '#DDDDDD') for c in sorted_cats]
labels = [func_labels.get(c, c) for c in sorted_cats]

bars = ax2.barh(y_pos, values, color=colors, edgecolor='white', linewidth=0.5)

# Add count and percentage labels
for i, (val, bar) in enumerate(zip(values, bars)):
    pct = val / total_novel * 100
    ax2.text(val + 3, i, f'{val} ({pct:.1f}%)', va='center', fontsize=8)

ax2.set_yticks(y_pos)
ax2.set_yticklabels(labels, fontsize=9)
ax2.invert_yaxis()
ax2.set_xlabel('Number of gene families', fontsize=10)
ax2.set_title('B. Functional categories — novel genes (736)', fontsize=11,
              fontweight='bold')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# Highlight key insight
ax2.text(0.95, 0.95, 'No spectinomycin\nresistance genes\nin novel set',
         transform=ax2.transAxes, fontsize=9, fontweight='bold',
         color='#CC0000', ha='right', va='top',
         bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFEEEE',
                   edgecolor='#CC0000', alpha=0.8))

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'fig6_pangenome_functional.png'), dpi=300,
            bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(OUT_DIR, 'fig6_pangenome_functional.pdf'),
            bbox_inches='tight', facecolor='white')
print("Saved fig6_pangenome_functional.png + .pdf")
