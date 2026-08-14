"""
Steps 6 + 8: Classify Panaroo gene families + functional descriptions
Author: Melis Gencel
Date: 2026-03-23
Input: results/genomics/comparative/panaroo_output/gene_presence_absence.csv
       results/genomics/comparative/M5_novel_regions.bed
       data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.tsv
Output:
  results/tables/genomics/panaroo_gene_classification.csv
  results/tables/genomics/novel_genes_by_function.csv
  results/tables/genomics/lost_genes_by_function.csv
  results/tables/genomics/shared_gene_functional_summary.csv
"""

import csv
import os
from collections import Counter, defaultdict

BASE = r"C:\Users\melis\Desktop\bioinformatics-portfolio\hgt-study"
PANAROO_CSV = os.path.join(BASE, "results", "genomics", "comparative",
                            "panaroo_output", "gene_presence_absence.csv")
BED_FILE = os.path.join(BASE, "results", "genomics", "comparative", "M5_novel_regions.bed")
M5_TSV = os.path.join(BASE, "data", "sequence_data", "TPHP3P_results",
                       "TPHP3P_1_M5__+_", "annotation", "TPHP3P_1_M5__+_.tsv")
OUT_DIR = os.path.join(BASE, "results", "tables", "genomics")
os.makedirs(OUT_DIR, exist_ok=True)


# ============================================================
# Build M5 locus_tag -> coordinate lookup
# ============================================================
m5_coords = {}
with open(M5_TSV, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 6 and parts[1] == 'cds':
            m5_coords[parts[5]] = (int(parts[2]), int(parts[3]))

# Parse BED for novel GAP regions
gap_regions = []
with open(BED_FILE, 'r') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 4 and parts[3] == 'GAP':
            gap_regions.append((int(parts[1]), int(parts[2])))


def in_novel_region(locus_tag):
    """Check if a locus tag falls within a novel GAP region."""
    if locus_tag not in m5_coords:
        return False
    start, stop = m5_coords[locus_tag]
    for gap_start, gap_stop in gap_regions:
        if start < gap_stop and stop > gap_start:
            return True
    return False


# ============================================================
# Functional category assignment
# ============================================================
def classify_function(annotation):
    """Assign functional category based on annotation keywords."""
    ann_lower = annotation.lower()

    # Order matters — more specific matches first
    if any(kw in ann_lower for kw in ['transposase', 'is110', 'is4', 'is3',
                                       'is200', 'is605', 'is701', 'is1595',
                                       'insertion sequence', 'tn3']):
        return 'IS_transposase'

    if any(kw in ann_lower for kw in ['phage', 'capsid', 'tail fiber',
                                       'tail protein', 'baseplate',
                                       'terminase', 'portal protein',
                                       'head morphogenesis', 'holin']):
        return 'phage'

    if any(kw in ann_lower for kw in ['conjugal', 'conjugation', 'vird',
                                       'virb', 'trbl', 'relaxase',
                                       'type iv secret', 'moba',
                                       'pilus assembly']):
        return 'conjugation_ICE'

    if any(kw in ann_lower for kw in ['integrase', 'recombinase',
                                       'site-specific recombination']):
        return 'integrase_recombinase'

    if any(kw in ann_lower for kw in ['crispr', 'cas1', 'cas2', 'cas3',
                                       'cas9', 'csm', 'cmr']):
        return 'CRISPR'

    if any(kw in ann_lower for kw in ['toxin', 'antitoxin', 'mazf',
                                       'maze', 'vapb', 'vapc',
                                       'parE', 'relE', 'doc']):
        return 'toxin_antitoxin'

    if any(kw in ann_lower for kw in ['efflux', 'multidrug', 'drug resistance',
                                       'acrb', 'mate efflux', 'norM',
                                       'smr family']):
        return 'drug_efflux'

    if any(kw in ann_lower for kw in ['abc transporter', 'mfs transporter',
                                       'permease', 'transporter',
                                       'antiporter', 'symporter']):
        return 'transport'

    if any(kw in ann_lower for kw in ['transcriptional regulator',
                                       'dna-binding', 'tetr', 'laci',
                                       'lysr', 'marr', 'merr', 'gntr',
                                       'arac', 'sigma factor',
                                       'anti-sigma']):
        return 'regulatory'

    if any(kw in ann_lower for kw in ['ribosomal protein', 'rrna',
                                       'trna', 'translation',
                                       'elongation factor']):
        return 'ribosomal_translation'

    if any(kw in ann_lower for kw in ['spore', 'sporulation',
                                       'germination', 'ger']):
        return 'sporulation'

    if any(kw in ann_lower for kw in ['kinase', 'dehydrogenase',
                                       'synthase', 'reductase',
                                       'transferase', 'hydrolase',
                                       'oxidoreductase', 'ligase',
                                       'lyase', 'isomerase', 'mutase',
                                       'epimerase', 'aldolase',
                                       'phosphatase', 'peptidase',
                                       'protease', 'esterase']):
        return 'metabolism_enzyme'

    if any(kw in ann_lower for kw in ['hypothetical', 'uncharacterized',
                                       'putative protein', 'duf']):
        return 'hypothetical_DUF'

    if any(kw in ann_lower for kw in ['membrane protein', 'cell wall',
                                       'peptidoglycan', 'lipopolysaccharide',
                                       'flagell', 'pilus', 'fimbri',
                                       'chemotaxis']):
        return 'cell_envelope_motility'

    return 'other'


# ============================================================
# Parse Panaroo CSV and classify
# ============================================================
print("Parsing Panaroo gene_presence_absence.csv...")

rows = []
with open(PANAROO_CSV, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        m5 = row['TPHP3P_1_M5__+_'].strip()
        m6 = row['TPHP3P_2_M6__+_'].strip()
        m8 = row['TPHP3P_4_M8__+_'].strip()
        i6 = row['I6_for_panaroo'].strip()
        annotation = row['Annotation'].strip()

        in_isolates = bool(m5) or bool(m6) or bool(m8)
        in_i6 = bool(i6)

        if in_isolates and in_i6:
            pattern = 'core'
        elif in_isolates and not in_i6:
            pattern = 'novel'
        elif not in_isolates and in_i6:
            pattern = 'lost'
        else:
            pattern = 'absent'  # shouldn't happen

        func_cat = classify_function(annotation)
        novel_region = in_novel_region(m5) if m5 else False

        rows.append({
            'Gene': row['Gene'],
            'Annotation': annotation,
            'Presence_pattern': pattern,
            'Functional_category': func_cat,
            'M5_locus': m5,
            'M6_locus': m6,
            'M8_locus': m8,
            'I6_locus': i6,
            'In_novel_region': 'yes' if novel_region else 'no'
        })

# Write classification CSV
out_path = os.path.join(OUT_DIR, 'panaroo_gene_classification.csv')
with open(out_path, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['Gene', 'Annotation', 'Presence_pattern',
                                            'Functional_category', 'M5_locus',
                                            'M6_locus', 'M8_locus', 'I6_locus',
                                            'In_novel_region'])
    writer.writeheader()
    writer.writerows(rows)

print(f"Wrote {len(rows)} gene families to {out_path}")

# ============================================================
# Summary statistics
# ============================================================
pattern_counts = Counter(r['Presence_pattern'] for r in rows)
print(f"\nPresence patterns: {dict(pattern_counts)}")

# Functional breakdown by pattern
for pattern in ['core', 'novel', 'lost']:
    subset = [r for r in rows if r['Presence_pattern'] == pattern]
    func_counts = Counter(r['Functional_category'] for r in subset)
    print(f"\n{pattern.upper()} ({len(subset)} genes):")
    for cat, cnt in func_counts.most_common():
        pct = cnt / len(subset) * 100
        print(f"  {cat}: {cnt} ({pct:.1f}%)")


# ============================================================
# Step 8: Detailed functional tables
# ============================================================

# Novel genes by function
novel = [r for r in rows if r['Presence_pattern'] == 'novel']
novel_out = os.path.join(OUT_DIR, 'novel_genes_by_function.csv')
with open(novel_out, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Functional_category', 'Count', 'In_novel_region',
                      'Representative_annotations'])
    func_groups = defaultdict(list)
    for r in novel:
        func_groups[r['Functional_category']].append(r)

    for cat in sorted(func_groups, key=lambda c: -len(func_groups[c])):
        genes = func_groups[cat]
        n_in_novel = sum(1 for g in genes if g['In_novel_region'] == 'yes')
        # Get unique annotations (up to 5)
        unique_anns = list(set(g['Annotation'] for g in genes))[:5]
        writer.writerow([cat, len(genes), n_in_novel,
                          '; '.join(unique_anns)])

print(f"\nWrote: {novel_out}")

# Lost genes by function
lost = [r for r in rows if r['Presence_pattern'] == 'lost']
lost_out = os.path.join(OUT_DIR, 'lost_genes_by_function.csv')
with open(lost_out, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Functional_category', 'Count', 'Representative_annotations'])
    func_groups = defaultdict(list)
    for r in lost:
        func_groups[r['Functional_category']].append(r)

    for cat in sorted(func_groups, key=lambda c: -len(func_groups[c])):
        genes = func_groups[cat]
        unique_anns = list(set(g['Annotation'] for g in genes))[:5]
        writer.writerow([cat, len(genes), '; '.join(unique_anns)])

print(f"Wrote: {lost_out}")

# Shared gene functional summary (comparing core vs novel vs lost)
summary_out = os.path.join(OUT_DIR, 'shared_gene_functional_summary.csv')
all_cats = sorted(set(r['Functional_category'] for r in rows))
with open(summary_out, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Functional_category', 'Core', 'Novel', 'Lost',
                      'Pct_of_novel', 'Enriched_in'])
    for cat in all_cats:
        core_n = sum(1 for r in rows if r['Presence_pattern'] == 'core' and r['Functional_category'] == cat)
        novel_n = sum(1 for r in rows if r['Presence_pattern'] == 'novel' and r['Functional_category'] == cat)
        lost_n = sum(1 for r in rows if r['Presence_pattern'] == 'lost' and r['Functional_category'] == cat)

        total_core = pattern_counts['core']
        total_novel = pattern_counts['novel']

        pct_core = core_n / total_core * 100 if total_core else 0
        pct_novel = novel_n / total_novel * 100 if total_novel else 0

        if pct_novel > pct_core * 2:
            enriched = 'novel'
        elif pct_core > pct_novel * 2:
            enriched = 'core'
        else:
            enriched = 'neither'

        writer.writerow([cat, core_n, novel_n, lost_n,
                          f"{pct_novel:.1f}%", enriched])

print(f"Wrote: {summary_out}")
print("\nDone — gene classification complete.")
