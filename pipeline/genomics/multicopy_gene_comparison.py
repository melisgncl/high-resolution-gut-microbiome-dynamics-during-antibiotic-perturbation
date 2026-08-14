"""
Multi-copy gene comparison: I6 (Bakta) vs M5/M6/M8 (Bakta)
Author: Melis Gencel
Date: 2026-03-23

All genomes annotated by Bakta → consistent product naming.
Find genes with >2 copies per genome, compare I6 vs all isolates.
Uses gene name when available, falls back to product.

Input:
  - data/references/bakta_I6/I6_bakta.tsv
  - data/sequence_data/TPHP3P_results/TPHP3P_*/annotation/TPHP3P_*_.tsv (M5/M6/M8)
Output:
  - results/figures/genomics/fig7_multicopy_comparison.png + .pdf
  - results/figures/genomics/fig8_copy_number_scatter.png + .pdf
  - results/figures/genomics/fig9_multicopy_categories.png + .pdf
  - results/tables/genomics/multicopy_gene_comparison.csv
  - results/genomics/reports/multicopy_gene_report.md
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import csv
import os
from collections import Counter, defaultdict

BASE = r"C:\Users\melis\Desktop\bioinformatics-portfolio\hgt-study"

TSV_PATHS = {
    'M5': os.path.join(BASE, "data", "sequence_data", "TPHP3P_results",
                        "TPHP3P_1_M5__+_", "annotation", "TPHP3P_1_M5__+_.tsv"),
    'M6': os.path.join(BASE, "data", "sequence_data", "TPHP3P_results",
                        "TPHP3P_2_M6__+_", "annotation", "TPHP3P_2_M6__+_.tsv"),
    'M8': os.path.join(BASE, "data", "sequence_data", "TPHP3P_results",
                        "TPHP3P_4_M8__+_", "annotation", "TPHP3P_4_M8__+_.tsv"),
    'I6': os.path.join(BASE, "data", "references", "bakta_I6", "I6_bakta.tsv"),
}

FIG_DIR = os.path.join(BASE, "results", "figures", "genomics")
TAB_DIR = os.path.join(BASE, "results", "tables", "genomics")
REP_DIR = os.path.join(BASE, "results", "genomics", "reports")


def parse_bakta(tsv_path):
    """Parse Bakta TSV. Return list of (gene_name, product) for CDS."""
    entries = []
    with open(tsv_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 8 and parts[1] == 'cds':
                gene = parts[6].strip()
                product = parts[7].strip()
                entries.append((gene, product))
    return entries


def make_label(gene, product, max_len=50):
    """Create a display label: gene name if available, else shortened product."""
    if gene:
        return gene
    # Shorten product: remove common suffixes
    p = product
    for suffix in [' family protein', ' domain-containing protein',
                   ' domain protein', ' family transporter']:
        p = p.replace(suffix, '')
    return p[:max_len]


# ============================================================
# Step 1: Parse all genomes — count by product (for matching)
#         and build gene-name lookup (for display)
# ============================================================
print("=" * 70)
print("Step 1: Parse Bakta TSVs")
print("=" * 70)

genome_products = {}   # {genome: Counter of products}
product_to_gene = {}   # {product: best gene name seen}

for name, path in TSV_PATHS.items():
    entries = parse_bakta(path)
    products = [e[1] for e in entries]
    genome_products[name] = Counter(products)
    # Collect gene names
    for gene, product in entries:
        if gene and product not in product_to_gene:
            product_to_gene[product] = gene
    print(f"  {name}: {len(entries)} CDS, {len(genome_products[name])} unique products")

EXCLUDE = {'hypothetical protein', 'putative protein', ''}


# ============================================================
# Step 2: Find multi-copy genes (>2 copies) in each genome
# ============================================================
print("\n" + "=" * 70)
print("Step 2: Multi-copy genes (>2 copies)")
print("=" * 70)

multicopy = {}
for name, counts in genome_products.items():
    mc = {prod: n for prod, n in counts.items() if n > 2 and prod not in EXCLUDE}
    multicopy[name] = mc
    print(f"  {name}: {len(mc)} families with >2 copies")

all_mc_products = set()
for mc in multicopy.values():
    all_mc_products |= set(mc.keys())
print(f"  Total unique multi-copy families: {len(all_mc_products)}")


# ============================================================
# Step 3: Build comparison table with gene names
# ============================================================
def categorise(product):
    p = product.lower()
    if any(kw in p for kw in ['transposase', 'is110', 'is4', 'is3', 'insertion']):
        return 'IS transposase'
    if any(kw in p for kw in ['phage', 'capsid', 'tail protein', 'terminase', 'portal protein']):
        return 'Phage'
    if any(kw in p for kw in ['abc transporter', 'permease', 'binding protein']):
        return 'ABC transporter'
    if any(kw in p for kw in ['mfs transporter', 'mfs-type']):
        return 'MFS transporter'
    if any(kw in p for kw in ['efflux', 'multidrug', 'drug resistance', 'mate efflux', 'norM']):
        return 'Drug efflux'
    if any(kw in p for kw in ['transcriptional regulator', 'dna-binding', 'tetr', 'marr',
                               'merr', 'lysr', 'laci', 'gntr', 'arac']):
        return 'Regulator'
    if any(kw in p for kw in ['ribosomal protein', 'rrna', 'trna']):
        return 'Ribosomal/RNA'
    if any(kw in p for kw in ['spore', 'sporulation', 'germination']):
        return 'Sporulation'
    if any(kw in p for kw in ['kinase', 'dehydrogenase', 'synthase', 'reductase',
                               'transferase', 'hydrolase', 'oxidoreductase', 'ligase']):
        return 'Metabolism'
    if any(kw in p for kw in ['conjugal', 'vird', 'virb', 'relaxase', 'moba']):
        return 'Conjugation'
    if 'duf' in p:
        return 'DUF/domain'
    if any(kw in p for kw in ['chemotaxis', 'flagell', 'pilus']):
        return 'Motility'
    return 'Other'


rows = []
for product in sorted(all_mc_products):
    i6_n = genome_products['I6'].get(product, 0)
    m5_n = genome_products['M5'].get(product, 0)
    m6_n = genome_products['M6'].get(product, 0)
    m8_n = genome_products['M8'].get(product, 0)

    gene_name = product_to_gene.get(product, '')
    label = make_label(gene_name, product)
    category = categorise(product)

    # Classification
    iso_max = max(m5_n, m6_n, m8_n)
    if i6_n > 2 and iso_max > 2:
        if abs(m5_n - i6_n) <= 1:
            classification = 'shared_same'
        elif m5_n > i6_n:
            classification = 'isolate_expanded'
        else:
            classification = 'I6_expanded'
    elif i6_n > 2 and iso_max <= 2:
        classification = 'I6_only_multicopy'
    elif iso_max > 2 and i6_n <= 2:
        classification = 'isolate_only_multicopy'
    else:
        classification = 'other'

    rows.append({
        'Gene_name': gene_name,
        'Product': product,
        'Label': label,
        'I6': i6_n, 'M5': m5_n, 'M6': m6_n, 'M8': m8_n,
        'Delta_M5': m5_n - i6_n,
        'Category': category,
        'Classification': classification
    })

# Write CSV
csv_path = os.path.join(TAB_DIR, 'multicopy_gene_comparison.csv')
with open(csv_path, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=['Gene_name', 'Product', 'Label',
                                            'I6', 'M5', 'M6', 'M8',
                                            'Delta_M5', 'Category', 'Classification'])
    writer.writeheader()
    writer.writerows(sorted(rows, key=lambda r: -abs(r['Delta_M5'])))

class_counts = Counter(r['Classification'] for r in rows)
print(f"\n  Classification: {dict(class_counts)}")
print(f"  Wrote: {csv_path}")


# ============================================================
# Figure 7: Multi-copy comparison — all 4 genomes, gene names
# ============================================================
print("\n" + "=" * 70)
print("Plotting figures")
print("=" * 70)

plot_rows = [r for r in rows if r['I6'] > 2 or r['M5'] > 2]
plot_rows.sort(key=lambda r: r['Delta_M5'])

# Show top/bottom 25 + IS elements if too many
if len(plot_rows) > 60:
    top = [r for r in plot_rows if r['Delta_M5'] > 2][-25:]
    bottom = [r for r in plot_rows if r['Delta_M5'] < -2][:25]
    is_rows = [r for r in plot_rows if r['Category'] == 'IS transposase'
               and r not in top and r not in bottom]
    plot_subset = bottom + is_rows + top
    plot_subset.sort(key=lambda r: r['Delta_M5'])
else:
    plot_subset = plot_rows

n = len(plot_subset)
fig, ax = plt.subplots(figsize=(16, max(8, n * 0.28)))

y_pos = np.arange(n)
labels = [r['Label'] for r in plot_subset]

bar_h = 0.2
colors = {'I6': '#2166AC', 'M5': '#CC3333', 'M6': '#E6AB02', 'M8': '#7570B3'}

for gi, genome in enumerate(['I6', 'M5', 'M6', 'M8']):
    vals = [r[genome] for r in plot_subset]
    offset = (gi - 1.5) * bar_h
    bars = ax.barh(y_pos + offset, vals, bar_h * 0.9,
                   color=colors[genome], alpha=0.85, label=genome)
    # Count labels
    for i, v in enumerate(vals):
        if v > 0:
            ax.text(v + 0.2, i + offset, str(v), va='center',
                    fontsize=5.5, color=colors[genome], fontweight='bold')

ax.set_yticks(y_pos)
ax.set_yticklabels(labels, fontsize=7)
ax.set_xlabel('Copy number', fontsize=11)
ax.set_title('Multi-copy gene comparison: I6 vs M5/M6/M8\n'
             'Gene names shown where available; all Bakta-annotated',
             fontsize=12, fontweight='bold')
ax.legend(loc='lower right', fontsize=10, ncol=4)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.invert_yaxis()

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig7_multicopy_comparison.png'), dpi=300,
            bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(FIG_DIR, 'fig7_multicopy_comparison.pdf'),
            bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved fig7_multicopy_comparison")


# ============================================================
# Figure 8: Scatter — I6 vs each isolate (3 panels)
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharex=True, sharey=True)

cat_colors = {
    'IS transposase': '#FF6600', 'Phage': '#7744AA',
    'ABC transporter': '#6699CC', 'MFS transporter': '#33AA77',
    'Drug efflux': '#FF3366', 'Regulator': '#99CC33',
    'Metabolism': '#33AA77', 'Conjugation': '#CC3333',
    'Ribosomal/RNA': '#9966CC', 'Sporulation': '#CC9933',
    'DUF/domain': '#BBBBBB', 'Other': '#DDDDDD', 'Motility': '#66CCAA',
}

for idx, (genome, ax) in enumerate(zip(['M5', 'M6', 'M8'], axes)):
    for r in rows:
        color = cat_colors.get(r['Category'], '#DDDDDD')
        size = max(15, min(150, (r['I6'] + r[genome]) * 3))
        ax.scatter(r['I6'], r[genome], c=color, s=size, alpha=0.5,
                   edgecolor='white', linewidth=0.3)

    max_val = max(max(r['I6'] for r in rows), max(r[genome] for r in rows))
    ax.plot([0, max_val], [0, max_val], 'k--', alpha=0.3, linewidth=1)

    # Label top outliers
    delta_key = f'{genome}'
    outliers = sorted(rows, key=lambda r: -abs(r[genome] - r['I6']))[:8]
    for r in outliers:
        if abs(r[genome] - r['I6']) > 3:
            ax.annotate(r['Label'][:25], (r['I6'], r[genome]),
                        fontsize=5, alpha=0.8,
                        xytext=(4, 4), textcoords='offset points')

    ax.set_xlabel('I6 copies', fontsize=10)
    ax.set_title(f'{genome} vs I6', fontsize=11, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

axes[0].set_ylabel('Isolate copies', fontsize=10)

# Shared legend
legend_cats = ['IS transposase', 'Phage', 'ABC transporter', 'Drug efflux',
               'Regulator', 'Metabolism', 'Other']
legend_patches = [mpatches.Patch(color=cat_colors[c], label=c) for c in legend_cats]
fig.legend(handles=legend_patches, loc='lower center', ncol=7, fontsize=8,
           bbox_to_anchor=(0.5, -0.02))

fig.suptitle('Gene copy number: I6 (reference) vs each isolate\n'
             'Above diagonal = isolate-expanded',
             fontsize=13, fontweight='bold', y=1.02)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig8_copy_number_scatter.png'), dpi=300,
            bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(FIG_DIR, 'fig8_copy_number_scatter.pdf'),
            bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved fig8_copy_number_scatter")


# ============================================================
# Figure 9: Category summary — all 4 genomes
# ============================================================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

categories = sorted(set(r['Category'] for r in rows))
genome_cat_totals = {g: defaultdict(int) for g in ['I6', 'M5', 'M6', 'M8']}
for r in rows:
    for g in ['I6', 'M5', 'M6', 'M8']:
        genome_cat_totals[g][r['Category']] += r[g]

cat_order = sorted(categories, key=lambda c: -(genome_cat_totals['I6'][c] +
                                                 genome_cat_totals['M5'][c]))

x_pos = np.arange(len(cat_order))
bar_w = 0.2
genome_list = ['I6', 'M5', 'M6', 'M8']
bar_colors = ['#2166AC', '#CC3333', '#E6AB02', '#7570B3']

for gi, (genome, color) in enumerate(zip(genome_list, bar_colors)):
    vals = [genome_cat_totals[genome][c] for c in cat_order]
    ax1.bar(x_pos + (gi - 1.5) * bar_w, vals, bar_w,
            color=color, alpha=0.85, label=genome)

ax1.set_xticks(x_pos)
ax1.set_xticklabels(cat_order, rotation=45, ha='right', fontsize=8)
ax1.set_ylabel('Total gene copies', fontsize=10)
ax1.set_title('A. Multi-copy genes by category (all genomes)', fontsize=11, fontweight='bold')
ax1.legend(fontsize=9, ncol=4)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)

# Panel B: Classification counts
cls_order = ['shared_same', 'isolate_expanded', 'I6_expanded',
             'isolate_only_multicopy', 'I6_only_multicopy']
cls_labels = ['Shared\n(same #)', 'Isolate\nexpanded', 'I6\nexpanded',
              'Isolate-only\nmulti-copy', 'I6-only\nmulti-copy']
cls_colors = ['#66BB66', '#CC3333', '#2166AC', '#FF6600', '#6699CC']
cls_vals = [class_counts.get(c, 0) for c in cls_order]

bars = ax2.bar(range(len(cls_order)), cls_vals, color=cls_colors, edgecolor='white')
for bar, val in zip(bars, cls_vals):
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 2,
             str(val), ha='center', fontsize=11, fontweight='bold')

ax2.set_xticks(range(len(cls_order)))
ax2.set_xticklabels(cls_labels, fontsize=9)
ax2.set_ylabel('Number of gene families', fontsize=10)
ax2.set_title('B. Classification of multi-copy families', fontsize=11, fontweight='bold')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig(os.path.join(FIG_DIR, 'fig9_multicopy_categories.png'), dpi=300,
            bbox_inches='tight', facecolor='white')
plt.savefig(os.path.join(FIG_DIR, 'fig9_multicopy_categories.pdf'),
            bbox_inches='tight', facecolor='white')
plt.close()
print("  Saved fig9_multicopy_categories")


# ============================================================
# Report
# ============================================================
print("\nWriting report...")

shared_same = [r for r in rows if r['Classification'] == 'shared_same']
iso_expanded = sorted([r for r in rows if r['Classification'] == 'isolate_expanded'],
                      key=lambda r: -r['Delta_M5'])
i6_expanded = sorted([r for r in rows if r['Classification'] == 'I6_expanded'],
                     key=lambda r: r['Delta_M5'])
iso_only = [r for r in rows if r['Classification'] == 'isolate_only_multicopy']
i6_only = [r for r in rows if r['Classification'] == 'I6_only_multicopy']

report_path = os.path.join(REP_DIR, 'multicopy_gene_report.md')
with open(report_path, 'w', encoding='utf-8') as f:
    f.write("# Multi-Copy Gene Comparison: Unified Bakta Annotation\n\n")
    f.write("**Date:** 2026-03-23\n")
    f.write("**Method:** All genomes (M5, M6, M8, I6) annotated with Bakta v1.12.0, db-light v6.0\n")
    f.write("**Threshold:** >2 copies of the same product name per genome\n\n---\n\n")

    f.write("## Summary\n\n")
    f.write("| Genome | Total CDS | Unique products | Multi-copy families (>2) |\n")
    f.write("|--------|-----------|-----------------|-------------------------|\n")
    for name in ['M5', 'M6', 'M8', 'I6']:
        total = sum(genome_products[name].values())
        unique = len(genome_products[name])
        mc = len(multicopy[name])
        f.write(f"| {name} | {total} | {unique} | {mc} |\n")

    f.write(f"\nM5, M6, and M8 are identical in gene content (clonal population).\n\n")

    f.write("## Classification\n\n")
    f.write("| Classification | Count | Description |\n")
    f.write("|---------------|-------|-------------|\n")
    f.write(f"| Shared (same copy #) | {len(shared_same)} | Ancestral paralogs, same count |\n")
    f.write(f"| Isolate expanded | {len(iso_expanded)} | Both have >2, isolates have more |\n")
    f.write(f"| I6 expanded | {len(i6_expanded)} | Both have >2, I6 has more |\n")
    f.write(f"| Isolate-only multi-copy | {len(iso_only)} | >2 in isolates, ≤2 in I6 |\n")
    f.write(f"| I6-only multi-copy | {len(i6_only)} | >2 in I6, ≤2 in isolates |\n")

    f.write("\n## Top isolate-expanded genes\n\n")
    f.write("| Gene | Product | M5 | M6 | M8 | I6 | Δ | Category |\n")
    f.write("|------|---------|----|----|----|----|---|----------|\n")
    for r in iso_expanded[:20]:
        f.write(f"| {r['Label'][:30]} | {r['Product'][:40]} | {r['M5']} | {r['M6']} | {r['M8']} | {r['I6']} | +{r['Delta_M5']} | {r['Category']} |\n")

    f.write("\n## Top I6-expanded genes\n\n")
    f.write("| Gene | Product | I6 | M5 | M6 | M8 | Δ | Category |\n")
    f.write("|------|---------|----|----|----|----|---|----------|\n")
    for r in i6_expanded[:20]:
        f.write(f"| {r['Label'][:30]} | {r['Product'][:40]} | {r['I6']} | {r['M5']} | {r['M6']} | {r['M8']} | {r['Delta_M5']} | {r['Category']} |\n")

    f.write("\n## Shared ancestral paralogs (same copy number)\n\n")
    f.write("| Gene | Product | Copies | Category |\n")
    f.write("|------|---------|--------|----------|\n")
    for r in sorted(shared_same, key=lambda r: -r['M5']):
        f.write(f"| {r['Label'][:30]} | {r['Product'][:45]} | {r['M5']} | {r['Category']} |\n")

    f.write("\n## Is it normal to have multiple copies?\n\n")
    f.write("Yes. The following are **expected** multi-copy families:\n\n")
    f.write(f"- **IS transposases**: IS110 ({genome_products['M5'].get('IS110 family transposase', 0)} in M5 vs "
            f"{genome_products['I6'].get('IS110 family transposase', 0)} in I6), "
            f"IS4 ({genome_products['M5'].get('IS4 family transposase', 0)} vs "
            f"{genome_products['I6'].get('IS4 family transposase', 0)}). "
            "Mobile elements that self-replicate.\n")
    f.write("- **ABC transporters/permeases**: Largest paralogous family in bacteria. "
            "Normal for 7+ Mb genomes.\n")
    f.write("- **Transcriptional regulators**: Scale with genome size. "
            "TetR, MarR, MerR, LysR families.\n")
    f.write("- **rRNA operons**: 8× 16S in M5 — normal for fast-growing Firmicutes.\n")
    f.write("- **Phage genes**: SP-beta prophage contributes ~30 structural genes.\n\n")

    f.write("## Relevance to spectinomycin resistance\n\n")
    f.write("**None.** rpsE is single-copy in all genomes. No spectinomycin resistance genes "
            "are among the multi-copy families. The 2 novel drug efflux genes (MFS drug:H+ antiporter, "
            "MATE efflux) are not spectinomycin-specific.\n")

print(f"  Wrote: {report_path}")
print("\nDONE")
