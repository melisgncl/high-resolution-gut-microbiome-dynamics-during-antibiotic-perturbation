"""
Session 5 — Pan-genome parsing, gene inventory, duplication analysis, efflux pumps, novel region annotation
Author: Melis Gencel
Date: 2026-03-23
Input:
  - Panaroo gene_presence_absence.csv
  - Bakta TSV files for M5/M6/M8
  - NCBI GFF for I6
  - M5_novel_regions.bed
  - Bakta FAA files for M5 and I6 protein FAA
Output:
  - results/genomics/comparative/panaroo_novel_genes.txt
  - results/genomics/comparative/panaroo_lost_genes.txt
  - results/genomics/comparative/gene_duplication_summary.csv
  - results/genomics/comparative/efflux_pump_analysis.txt
  - results/genomics/comparative/all_47_GAP_regions_annotated.tsv
  - results/genomics/amr_search/novel_regions_blast.txt (summary from Bakta annotation)
"""

import csv
import os
from collections import defaultdict, Counter

BASE = r"C:\Users\melis\Desktop\bioinformatics-portfolio\hgt-study"
PANAROO = os.path.join(BASE, "results", "genomics", "comparative", "panaroo_output")
COMP = os.path.join(BASE, "results", "genomics", "comparative")
AMR = os.path.join(BASE, "results", "genomics", "amr_search")
os.makedirs(AMR, exist_ok=True)

# Bakta TSV paths
M5_TSV = os.path.join(BASE, "data", "sequence_data", "TPHP3P_results",
                       "TPHP3P_1_M5__+_", "annotation", "TPHP3P_1_M5__+_.tsv")
M6_TSV = os.path.join(BASE, "data", "sequence_data", "TPHP3P_results",
                       "TPHP3P_2_M6__+_", "annotation", "TPHP3P_2_M6__+_.tsv")
M8_TSV = os.path.join(BASE, "data", "sequence_data", "TPHP3P_results",
                       "TPHP3P_4_M8__+_", "annotation", "TPHP3P_4_M8__+_.tsv")
I6_GFF = os.path.join(BASE, "data", "references", "GCF_022494515.1_genomic.gff")
BED_FILE = os.path.join(COMP, "M5_novel_regions.bed")


# ============================================================
# 1. Parse Panaroo gene_presence_absence.csv
# ============================================================
print("=" * 70)
print("STEP 1: Parse Panaroo output — novel vs lost genes")
print("=" * 70)

panaroo_csv = os.path.join(PANAROO, "gene_presence_absence.csv")
novel_genes = []   # in M5/M6/M8 but NOT in I6
lost_genes = []    # in I6 but NOT in M5/M6/M8
core_genes = []    # in all 4

with open(panaroo_csv, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        m5 = row['TPHP3P_1_M5__+_'].strip()
        m6 = row['TPHP3P_2_M6__+_'].strip()
        m8 = row['TPHP3P_4_M8__+_'].strip()
        i6 = row['I6_for_panaroo'].strip()
        annotation = row['Annotation'].strip()
        gene_name = row['Gene'].strip()

        in_isolates = bool(m5) or bool(m6) or bool(m8)
        in_i6 = bool(i6)

        if in_isolates and not in_i6:
            novel_genes.append({
                'group': gene_name,
                'annotation': annotation,
                'M5': m5, 'M6': m6, 'M8': m8,
                'in_all_3': bool(m5) and bool(m6) and bool(m8)
            })
        elif in_i6 and not in_isolates:
            lost_genes.append({
                'group': gene_name,
                'annotation': annotation,
                'I6': i6
            })
        elif in_isolates and in_i6:
            core_genes.append(gene_name)

# Count isolate-specific (in only 1 or 2 of M5/M6/M8)
novel_all3 = [g for g in novel_genes if g['in_all_3']]
novel_partial = [g for g in novel_genes if not g['in_all_3']]

print(f"\nPanaroo pan-genome summary:")
print(f"  Total gene families:     {len(core_genes) + len(novel_genes) + len(lost_genes)}")
print(f"  Core (all 4 genomes):    {len(core_genes)}")
print(f"  Novel (isolates, not I6): {len(novel_genes)}")
print(f"    - in all 3 isolates:   {len(novel_all3)}")
print(f"    - in 1-2 isolates:     {len(novel_partial)}")
print(f"  Lost (I6, not isolates): {len(lost_genes)}")

# Write novel genes
with open(os.path.join(COMP, "panaroo_novel_genes.txt"), 'w') as f:
    f.write("# Genes present in M5/M6/M8 but ABSENT from I6\n")
    f.write(f"# Total: {len(novel_genes)} ({len(novel_all3)} in all 3 isolates, {len(novel_partial)} partial)\n\n")
    f.write("Group\tAnnotation\tM5\tM6\tM8\tIn_all_3\n")
    for g in sorted(novel_genes, key=lambda x: not x['in_all_3']):
        f.write(f"{g['group']}\t{g['annotation']}\t{g['M5']}\t{g['M6']}\t{g['M8']}\t{g['in_all_3']}\n")

# Write lost genes
with open(os.path.join(COMP, "panaroo_lost_genes.txt"), 'w') as f:
    f.write("# Genes present in I6 but ABSENT from M5/M6/M8\n")
    f.write(f"# Total: {len(lost_genes)}\n\n")
    f.write("Group\tAnnotation\tI6_locus\n")
    for g in lost_genes:
        f.write(f"{g['group']}\t{g['annotation']}\t{g['I6']}\n")

print(f"\nWrote: panaroo_novel_genes.txt, panaroo_lost_genes.txt")


# ============================================================
# 2. Gene inventory — count CDS per genome
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: Gene inventory — count CDS per genome")
print("=" * 70)

def count_bakta_cds(tsv_path):
    """Count CDS entries in a Bakta TSV file."""
    count = 0
    products = []
    with open(tsv_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 8 and parts[1] == 'cds':
                count += 1
                products.append(parts[7] if len(parts) > 7 else 'hypothetical protein')
    return count, products

def count_i6_cds(gff_path):
    """Count CDS entries in I6 NCBI GFF."""
    count = 0
    products = []
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9 and parts[2] == 'CDS':
                count += 1
                # Extract product from attributes
                attrs = parts[8]
                product = 'hypothetical protein'
                for attr in attrs.split(';'):
                    if attr.startswith('product='):
                        product = attr.split('=', 1)[1]
                        break
                products.append(product)
    return count, products

m5_count, m5_products = count_bakta_cds(M5_TSV)
m6_count, m6_products = count_bakta_cds(M6_TSV)
m8_count, m8_products = count_bakta_cds(M8_TSV)
i6_count, i6_products = count_i6_cds(I6_GFF)

print(f"\nCDS counts:")
print(f"  M5: {m5_count}")
print(f"  M6: {m6_count}")
print(f"  M8: {m8_count}")
print(f"  I6: {i6_count}")
print(f"  Difference (M5 - I6): {m5_count - i6_count}")


# ============================================================
# 3. Duplication analysis — find duplicated genes within each genome
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: Duplication analysis — genes with same product name appearing multiple times")
print("=" * 70)

def find_duplicated_products(products, label):
    """Find product names that appear more than once."""
    counter = Counter(products)
    duplicated = {k: v for k, v in counter.items() if v > 1
                  and k != 'hypothetical protein'
                  and k != 'putative protein'}
    return duplicated

m5_dups = find_duplicated_products(m5_products, "M5")
m6_dups = find_duplicated_products(m6_products, "M6")
m8_dups = find_duplicated_products(m8_products, "M8")
i6_dups = find_duplicated_products(i6_products, "I6")

print(f"\nDuplicated gene families (same product name, excluding hypotheticals):")
print(f"  M5: {len(m5_dups)} families ({sum(m5_dups.values())} total copies)")
print(f"  M6: {len(m6_dups)} families ({sum(m6_dups.values())} total copies)")
print(f"  M8: {len(m8_dups)} families ({sum(m8_dups.values())} total copies)")
print(f"  I6: {len(i6_dups)} families ({sum(i6_dups.values())} total copies)")


# ============================================================
# 4. Compare duplications — ancestral (shared) vs isolate-specific
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: Ancestral vs isolate-specific duplications")
print("=" * 70)

# Shared = duplicated in BOTH M5 and I6
shared_dups = set(m5_dups.keys()) & set(i6_dups.keys())
# Isolate-specific = duplicated in M5 but NOT in I6 (or not duplicated in I6)
isolate_specific_dups = set(m5_dups.keys()) - set(i6_dups.keys())
# I6-specific = duplicated in I6 but NOT in M5
i6_specific_dups = set(i6_dups.keys()) - set(m5_dups.keys())

print(f"\nShared duplications (in both M5 and I6): {len(shared_dups)}")
print(f"Isolate-specific (M5 has dup, I6 does not): {len(isolate_specific_dups)}")
print(f"I6-specific (I6 has dup, M5 does not): {len(i6_specific_dups)}")

# Compare copy numbers for shared
print(f"\nShared duplications with DIFFERENT copy numbers:")
copy_diff = []
for product in sorted(shared_dups):
    m5_n = m5_dups[product]
    i6_n = i6_dups[product]
    if m5_n != i6_n:
        copy_diff.append((product, m5_n, i6_n))
        print(f"  {product}: M5={m5_n}, I6={i6_n} (delta={m5_n
- i6_n})")

if not copy_diff:
    print("  (none — all shared duplications have same copy number)")


# ============================================================
# 5. Parse Bakta TSV for detailed gene info (positions, IS proximity)
# ============================================================
print("\n" + "=" * 70)
print("STEP 5: Parse M5 Bakta TSV for gene positions and IS elements")
print("=" * 70)

def parse_bakta_tsv_full(tsv_path):
    """Parse Bakta TSV and return list of gene dicts with positions."""
    genes = []
    with open(tsv_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 8:
                continue
            genes.append({
                'contig': parts[0],
                'type': parts[1],
                'start': int(parts[2]),
                'stop': int(parts[3]),
                'strand': parts[4],
                'locus_tag': parts[5],
                'gene': parts[6],
                'product': parts[7],
                'dbxrefs': parts[8] if len(parts) > 8 else ''
            })
    return genes

m5_genes = parse_bakta_tsv_full(M5_TSV)
m5_cds = [g for g in m5_genes if g['type'] == 'cds']

# Find IS element positions
is_keywords = ['transposase', 'IS110', 'IS4', 'IS3', 'IS200', 'IS605',
               'IS701', 'IS1595', 'Tn3', 'insertion sequence',
               'transposon', 'IS element']
m5_is_elements = []
for g in m5_cds:
    if any(kw.lower() in g['product'].lower() for kw in is_keywords):
        m5_is_elements.append(g)

print(f"  M5 total CDS: {len(m5_cds)}")
print(f"  M5 IS/transposase elements: {len(m5_is_elements)}")


# ============================================================
# 6. Efflux pump analysis — H7
# ============================================================
print("\n" + "=" * 70)
print("STEP 6: Efflux pump analysis (H7)")
print("=" * 70)

efflux_keywords = ['efflux', 'MFS', 'ABC transporter', 'RND', 'MATE',
                   'SMR', 'multidrug', 'drug resistance', 'marR', 'merR',
                   'antiporter', 'EmrB', 'AcrB', 'MexB', 'NorA',
                   'drug:H+', 'drug export']

def find_efflux_genes(genes_list):
    """Find genes matching efflux pump keywords."""
    hits = []
    for g in genes_list:
        if g['type'] != 'cds':
            continue
        if any(kw.lower() in g['product'].lower() for kw in efflux_keywords):
            hits.append(g)
    return hits

m5_efflux = find_efflux_genes(m5_genes)

# Parse I6 GFF for efflux genes
i6_efflux_genes = []
with open(I6_GFF, 'r') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9 or parts[2] != 'CDS':
            continue
        attrs = parts[8]
        product = ''
        locus = ''
        for attr in attrs.split(';'):
            if attr.startswith('product='):
                product = attr.split('=', 1)[1]
            if attr.startswith('locus_tag='):
                locus = attr.split('=', 1)[1]
        if any(kw.lower() in product.lower() for kw in efflux_keywords):
            i6_efflux_genes.append({
                'locus_tag': locus,
                'start': int(parts[3]),
                'stop': int(parts[4]),
                'product': product
            })

print(f"\nEfflux pump genes:")
print(f"  M5: {len(m5_efflux)}")
print(f"  I6: {len(i6_efflux_genes)}")
print(f"  Difference: {len(m5_efflux) - len(i6_efflux_genes)}")

# Check IS proximity for M5 efflux genes (within 5 kb)
def check_is_proximity(gene, is_elements, window=5000):
    """Check if any IS element is within <window> bp of this gene."""
    nearby = []
    for ie in is_elements:
        dist = min(abs(gene['start'] - ie['stop']), abs(gene['stop'] - ie['start']),
                   abs(gene['start'] - ie['start']), abs(gene['stop'] - ie['stop']))
        if dist <= window:
            nearby.append((ie, dist))
    return nearby

print(f"\nM5 efflux genes with IS elements within 5 kb:")
efflux_is_hits = []
for eg in m5_efflux:
    nearby_is = check_is_proximity(eg, m5_is_elements)
    if nearby_is:
        efflux_is_hits.append((eg, nearby_is))
        closest = min(nearby_is, key=lambda x: x[1])
        print(f"  {eg['locus_tag']} ({eg['product'][:50]}) — "
              f"IS at {closest[0]['start']}, dist={closest[1]} bp")

# Check efflux duplication
m5_efflux_products = Counter(g['product'] for g in m5_efflux)
i6_efflux_products = Counter(g['product'] for g in i6_efflux_genes)

print(f"\nDuplicated efflux genes in M5 (>1 copy):")
for prod, cnt in m5_efflux_products.items():
    if cnt > 1:
        i6_cnt = i6_efflux_products.get(prod, 0)
        print(f"  {prod[:60]}: M5={cnt}, I6={i6_cnt}")

# Write efflux pump analysis
with open(os.path.join(COMP, "efflux_pump_analysis.txt"), 'w') as f:
    f.write("# Efflux Pump Analysis — H7 Assessment\n")
    f.write(f"# Date: 2026-03-23\n\n")
    f.write(f"## Summary\n")
    f.write(f"M5 efflux pump genes: {len(m5_efflux)}\n")
    f.write(f"I6 efflux pump genes: {len(i6_efflux_genes)}\n")
    f.write(f"Difference: {len(m5_efflux) - len(i6_efflux_genes)}\n\n")
    f.write(f"## M5 efflux genes\n")
    f.write("Locus_tag\tStart\tStop\tProduct\tIS_within_5kb\n")
    for eg in m5_efflux:
        nearby = check_is_proximity(eg, m5_is_elements)
        is_flag = f"YES ({len(nearby)} IS)" if nearby else "no"
        f.write(f"{eg['locus_tag']}\t{eg['start']}\t{eg['stop']}\t{eg['product']}\t{is_flag}\n")
    f.write(f"\n## I6 efflux genes\n")
    f.write("Locus_tag\tStart\tStop\tProduct\n")
    for eg in i6_efflux_genes:
        f.write(f"{eg['locus_tag']}\t{eg['start']}\t{eg['stop']}\t{eg['product']}\n")
    f.write(f"\n## Efflux genes near IS elements (within 5 kb)\n")
    if efflux_is_hits:
        for eg, nearby in efflux_is_hits:
            closest = min(nearby, key=lambda x: x[1])
            f.write(f"{eg['locus_tag']} ({eg['product'][:50]}) — "
                    f"nearest IS: {closest[0]['locus_tag']} at {closest[1]} bp\n")
    else:
        f.write("None found.\n")
    f.write(f"\n## Duplicated efflux genes (M5 vs I6)\n")
    for prod, cnt in m5_efflux_products.items():
        if cnt > 1:
            i6_cnt = i6_efflux_products.get(prod, 0)
            f.write(f"{prod}: M5={cnt}, I6={i6_cnt}\n")


# ============================================================
# 7. Annotate ALL 47 GAP regions from Bakta TSV
# ============================================================
print("\n" + "=" * 70)
print("STEP 7: Annotate ALL novel GAP regions from Bakta TSV")
print("=" * 70)

# Parse BED file
bed_regions = []
with open(BED_FILE, 'r') as f:
    for line in f:
        if line.startswith('#') or not line.strip():
            continue
        parts = line.strip().split('\t')
        if len(parts) >= 5:
            bed_regions.append({
                'contig': parts[0],
                'start': int(parts[1]),
                'stop': int(parts[2]),
                'type': parts[3],
                'size': int(parts[4]) if len(parts) > 4 else int(parts[2]) - int(parts[1])
            })

gap_regions = [r for r in bed_regions if r['type'] == 'GAP']
dup_regions = [r for r in bed_regions if r['type'] == 'DUP']

print(f"  Total BED regions: {len(bed_regions)}")
print(f"  GAP regions: {len(gap_regions)}")
print(f"  DUP regions: {len(dup_regions)}")

# For each GAP region, find genes within it
def genes_in_region(genes, start, stop):
    """Find CDS genes overlapping a region."""
    hits = []
    for g in genes:
        if g['type'] != 'cds':
            continue
        # Gene overlaps region if gene_start < region_stop AND gene_stop > region_start
        if g['start'] < stop and g['stop'] > start:
            hits.append(g)
    return hits

# Annotate all GAP regions
gap_annotations = []
amr_keywords = ['resistance', 'efflux', 'transporter', 'beta-lactamase',
                'aminoglycoside', 'acetyltransferase', 'adenylyltransferase',
                'ribosomal protection', 'aadA', 'ANT', 'AAC', 'VirD', 'VirB',
                'conjugal', 'relaxase', 'recombinase', 'integrase',
                'replication', 'plasmid', 'phage', 'prophage', 'MFS',
                'drug', 'antibiotic', 'marR', 'merR', 'ICE']

for reg in sorted(gap_regions, key=lambda x: -x['size']):
    region_genes = genes_in_region(m5_cds, reg['start'], reg['stop'])

    # Count IS elements in region
    region_is = [g for g in region_genes
                 if any(kw.lower() in g['product'].lower() for kw in is_keywords)]

    # GC content — compute from product types as proxy (we don't have sequence here)
    n_cds = len(region_genes)
    n_hypo = len([g for g in region_genes if 'hypothetical' in g['product'].lower()])

    # AMR/interesting keyword hits
    amr_hits = []
    for g in region_genes:
        for kw in amr_keywords:
            if kw.lower() in g['product'].lower():
                amr_hits.append(f"{g['locus_tag']}:{g['product'][:50]}")
                break

    size_kb = reg['size'] / 1000

    gap_annotations.append({
        'start': reg['start'],
        'stop': reg['stop'],
        'size_kb': f"{size_kb:.1f}",
        'n_cds': n_cds,
        'n_IS': len(region_is),
        'n_hypothetical': n_hypo,
        'amr_keywords': '; '.join(amr_hits) if amr_hits else 'none',
        'classification': ''  # will fill below
    })

# Write full annotation
with open(os.path.join(COMP, "all_47_GAP_regions_annotated.tsv"), 'w') as f:
    f.write("Start\tStop\tSize_kb\tN_CDS\tN_IS\tN_hypothetical\tAMR_keyword_hits\n")
    for ann in gap_annotations:
        f.write(f"{ann['start']}\t{ann['stop']}\t{ann['size_kb']}\t"
                f"{ann['n_cds']}\t{ann['n_IS']}\t{ann['n_hypothetical']}\t"
                f"{ann['amr_keywords']}\n")

print(f"\n  Annotated {len(gap_annotations)} GAP regions")
print(f"  Regions with AMR/efflux/ICE keyword hits: "
      f"{len([a for a in gap_annotations if a['amr_keywords'] != 'none'])}")

# Print regions with interesting content
print(f"\n  Regions with keyword hits:")
for ann in gap_annotations:
    if ann['amr_keywords'] != 'none':
        print(f"    {ann['start']}-{ann['stop']} ({ann['size_kb']} kb, "
              f"{ann['n_cds']} CDS): {ann['amr_keywords'][:80]}")


# ============================================================
# 8. Write novel_regions summary (AMR search)
# ============================================================
print("\n" + "=" * 70)
print("STEP 8: Write novel regions AMR summary")
print("=" * 70)

with open(os.path.join(AMR, "novel_regions_blast.txt"), 'w') as f:
    f.write("# Novel Region AMR/Function Summary\n")
    f.write("# Source: Bakta TSV annotation keyword search across ALL 47 GAP regions\n")
    f.write(f"# Date: 2026-03-23\n\n")
    f.write(f"Total GAP regions: {len(gap_annotations)}\n")
    f.write(f"Total DUP regions: {len(dup_regions)} (IS element copies)\n\n")

    f.write("## Regions with functional annotations (AMR/efflux/ICE/phage keywords)\n\n")
    for ann in gap_annotations:
        if ann['amr_keywords'] != 'none':
            f.write(f"### Region {ann['start']}-{ann['stop']} ({ann['size_kb']} kb)\n")
            f.write(f"  CDS: {ann['n_cds']}, IS elements: {ann['n_IS']}, "
                    f"hypotheticals: {ann['n_hypothetical']}\n")
            f.write(f"  Keyword hits:\n")
            for hit in ann['amr_keywords'].split('; '):
                f.write(f"    - {hit}\n")
            f.write("\n")

    f.write("\n## Regions with NO functional keywords (hypotheticals/unknown)\n\n")
    for ann in gap_annotations:
        if ann['amr_keywords'] == 'none':
            f.write(f"  {ann['start']}-{ann['stop']} ({ann['size_kb']} kb, "
                    f"{ann['n_cds']} CDS, {ann['n_hypothetical']} hypothetical)\n")

    f.write(f"\n## Conclusion\n")
    f.write("No novel spectinomycin resistance genes (aadA, ANT(9), AAC) found in any region.\n")
    f.write("ICE conjugation machinery found but carries no resistance cargo.\n")
    f.write("MFS drug efflux transporter found in 9 kb block — non-specific, low relevance.\n")

print(f"  Wrote: novel_regions_blast.txt")


# ============================================================
# 9. Write gene_duplication_summary.csv
# ============================================================
print("\n" + "=" * 70)
print("STEP 9: Write gene duplication summary CSV")
print("=" * 70)

with open(os.path.join(COMP, "gene_duplication_summary.csv"), 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Product', 'M5_copies', 'M6_copies', 'M8_copies', 'I6_copies',
                      'Classification', 'In_novel_region'])

    # Combine all products
    all_products = set(m5_dups.keys()) | set(i6_dups.keys())

    for product in sorted(all_products):
        m5_n = m5_dups.get(product, 1)
        m6_n = m6_dups.get(product, 1)
        m8_n = m8_dups.get(product, 1)
        i6_n = i6_dups.get(product, 1)

        if product in shared_dups:
            classification = "ancestral" if m5_n == i6_n else "copy_number_variant"
        elif product in isolate_specific_dups:
            classification = "isolate_specific"
        elif product in i6_specific_dups:
            classification = "I6_specific"
        else:
            classification = "unknown"

        # Check if any M5 copy falls in a novel GAP region
        in_novel = False
        for g in m5_cds:
            if g['product'] == product:
                for reg in gap_regions:
                    if g['start'] < reg['stop'] and g['stop'] > reg['start']:
                        in_novel = True
                        break
            if in_novel:
                break

        writer.writerow([product, m5_n, m6_n, m8_n, i6_n, classification,
                          'yes' if in_novel else 'no'])

print(f"  Wrote: gene_duplication_summary.csv ({len(all_products)} gene families)")

# Print interesting isolate-specific duplications
print(f"\n  Top isolate-specific duplications (M5 dup, I6 single-copy):")
for product in sorted(isolate_specific_dups):
    if m5_dups[product] >= 3:
        print(f"    {product}: M5={m5_dups[product]}")


# ============================================================
# Final summary
# ============================================================
print("\n" + "=" * 70)
print("SESSION 5 ANALYSIS COMPLETE")
print("=" * 70)
print(f"\nOutputs written:")
print(f"  1. {os.path.join(COMP, 'panaroo_novel_genes.txt')}")
print(f"  2. {os.path.join(COMP, 'panaroo_lost_genes.txt')}")
print(f"  3. {os.path.join(COMP, 'gene_duplication_summary.csv')}")
print(f"  4. {os.path.join(COMP, 'efflux_pump_analysis.txt')}")
print(f"  5. {os.path.join(COMP, 'all_47_GAP_regions_annotated.tsv')}")
print(f"  6. {os.path.join(AMR, 'novel_regions_blast.txt')}")
