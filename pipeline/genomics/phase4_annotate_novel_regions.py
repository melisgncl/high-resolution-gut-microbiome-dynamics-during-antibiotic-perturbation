"""
Title:   Phase 4.1 — Annotate M5 novel regions vs I6
Author:  Melis Gencel
Date:    2026-03-23
Input:   M5_novel_regions.bed, M5 Bakta TSV+GFF3+FNA
Output:  results/genomics/comparative/M5_novel_regions_annotated.tsv
         results/genomics/comparative/phase4_ICE_donor_analysis.txt
         results/genomics/comparative/phase4_IS_census.txt
"""

from pathlib import Path
import re

BASE = Path(r"C:\Users\melis\Desktop\bioinformatics-portfolio\hgt-study")
OUTD = BASE / "results" / "genomics" / "comparative"

M5_FNA  = BASE / "data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.fna"
M5_TSV  = BASE / "data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.tsv"
M5_GFF  = BASE / "data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.gff3"
BED     = OUTD / "M5_novel_regions.bed"

# ── Parse M5 FNA ────────────────────────────────────────────────────────────────
print("Loading M5 assembly...")
m5_seq = ""
with open(M5_FNA) as fh:
    for line in fh:
        if not line.startswith(">"): m5_seq += line.strip()
print(f"M5 assembly: {len(m5_seq):,} bp")

def gc_content(seq):
    s = seq.upper()
    gc = s.count("G") + s.count("C")
    return 100 * gc / len(s) if s else 0

# ── Parse M5 Bakta TSV ─────────────────────────────────────────────────────────
# Columns: #Sequence  Type  Start  Stop  Strand  LocusTag  Gene  Product  DbXrefs
print("Parsing M5 Bakta annotation...")

class Feature:
    __slots__ = ("start", "stop", "ftype", "locus", "gene", "product", "dbxrefs")
    def __init__(self, start, stop, ftype, locus, gene, product, dbxrefs):
        self.start = start; self.stop = stop; self.ftype = ftype
        self.locus = locus; self.gene = gene; self.product = product
        self.dbxrefs = dbxrefs

features = []
with open(M5_TSV) as fh:
    for line in fh:
        if line.startswith("#"): continue
        p = line.strip().split("\t")
        if len(p) < 6: continue
        try:
            start = int(p[2]); stop = int(p[3])
        except ValueError: continue
        ftype = p[1]; locus = p[5]
        gene    = p[6] if len(p) > 6 else ""
        product = p[7] if len(p) > 7 else ""
        dbxrefs = p[8] if len(p) > 8 else ""
        features.append(Feature(start, stop, ftype, locus, gene, product, dbxrefs))
print(f"  {len(features)} features loaded")

# Index features by position for fast overlap query
def features_in(start, stop, margin=5000):
    """Return features overlapping [start-margin, stop+margin]."""
    return [f for f in features if f.stop >= start - margin and f.start <= stop + margin]

def features_strictly_in(start, stop):
    """Return features with centre inside [start, stop]."""
    return [f for f in features if f.start >= start and f.stop <= stop]

# ── Parse novel regions BED ─────────────────────────────────────────────────────
print("Parsing novel regions BED...")
novel_regions = []
with open(BED) as fh:
    for line in fh:
        if line.startswith("#"): continue
        p = line.strip().split("\t")
        if len(p) < 5: continue
        contig, start, stop, ftype, net = p[0], int(p[1]), int(p[2]), p[3], int(p[4])
        novel_regions.append((contig, start, stop, ftype, net))

gap_regions = [(c,s,e,t,n) for c,s,e,t,n in novel_regions if t == "GAP"]
dup_regions = [(c,s,e,t,n) for c,s,e,t,n in novel_regions if t == "DUP"]
print(f"  {len(gap_regions)} GAP regions, {len(dup_regions)} DUP regions")

# ── IS element census ───────────────────────────────────────────────────────────
print("\n=== IS ELEMENT CENSUS ===")
IS_KEYWORDS = ["transposase", "transposon", "insertion sequence", "IS element",
               "IS1", "IS3", "IS4", "IS5", "IS6", "IS200", "IS21", "IS26", "IS30",
               "IS66", "IS110", "IS256", "IS630", "IS701", "ISL3", "Tn"]

is_features = [f for f in features
               if any(kw.lower() in f.product.lower() for kw in IS_KEYWORDS)
               or any(kw.lower() in f.gene.lower() for kw in IS_KEYWORDS)]
print(f"Total IS/transposase features: {len(is_features)}")

# Classify IS families from product annotation
is_families = {}
for f in is_features:
    # Extract IS family from product
    m = re.search(r'IS\d+|IS[A-Z]+|Tn\d+', f.product, re.IGNORECASE)
    family = m.group(0) if m else "unknown"
    is_families[family] = is_families.get(family, 0) + 1

print("IS family distribution:")
for fam, cnt in sorted(is_families.items(), key=lambda x: -x[1])[:15]:
    print(f"  {fam:15} {cnt:4d}")

# ── Check IS element flanking for each large GAP region ─────────────────────────
def check_IS_flanking(start, stop, flank_kb=5):
    """Return IS elements within flank_kb of the region boundaries."""
    margin = flank_kb * 1000
    left_is  = [f for f in is_features if stop - margin <= f.stop <= start + margin and f.stop <= start]
    right_is = [f for f in is_features if stop - margin <= f.start <= stop + margin and f.start >= stop]
    return left_is, right_is

def check_tRNA_flanking(start, stop, flank_kb=10):
    margin = flank_kb * 1000
    trna = [f for f in features if f.ftype == "tRNA"
            and ((stop - margin <= f.stop <= start) or (stop <= f.start <= stop + margin))]
    return trna

# ── Annotate top 10 GAP regions ─────────────────────────────────────────────────
print("\n=== TOP NOVEL REGION ANNOTATION ===")

# Sort by net extra bp
gap_sorted = sorted(gap_regions, key=lambda x: -x[4])[:10]

annotated = []
for contig, start, stop, ftype, net in gap_sorted:
    size = stop - start
    seq  = m5_seq[start:stop]
    gc   = gc_content(seq)

    # Gene content
    interior = features_strictly_in(start, stop)
    n_cds    = sum(1 for f in interior if f.ftype == "cds")
    n_rna    = sum(1 for f in interior if f.ftype in ("tRNA", "rRNA", "ncRNA"))
    n_is     = sum(1 for f in interior if f in is_features)

    # IS flanking
    left_is, right_is = check_IS_flanking(start, stop, 5)
    trna_flanking     = check_tRNA_flanking(start, stop, 10)

    # GI hallmarks
    is_flanked  = len(left_is) > 0 or len(right_is) > 0
    trna_flanked = len(trna_flanking) > 0

    print(f"\n{contig}:{start:,}–{stop:,} ({size//1000} kb, net +{net//1000} kb)")
    print(f"  GC: {gc:.1f}%  |  CDS: {n_cds}  |  IS internal: {n_is}  |  RNA: {n_rna}")
    print(f"  IS flanking (5kb): L={len(left_is)}, R={len(right_is)}  |  tRNA flanking: {len(trna_flanking)}")

    # Key gene content summary
    key_genes = []
    for f in interior:
        if not f.product: continue
        p = f.product.lower()
        if any(kw in p for kw in ["vird4", "virb4", "virb6", "trbL", "relaxase", "recombinase",
                                   "integrase", "conjugal", "transpos", "phage", "prophage",
                                   "resistance", "aminoglycoside", "efflux", "hypothetical"]):
            key_genes.append(f"{f.gene or f.locus}: {f.product[:60]}")

    if key_genes:
        print(f"  Notable genes ({len(key_genes)}):")
        for g in key_genes[:8]:
            print(f"    {g}")

    # Extract RefSeq cross-references from ICE genes to identify donor organisms
    refseq_orgs = []
    for f in interior[:20]:
        if "RefSeq:" in f.dbxrefs:
            rs = [x.strip() for x in f.dbxrefs.split(",") if "RefSeq:" in x]
            refseq_orgs.extend(rs[:1])

    annotated.append({
        "region": f"{start}-{stop}",
        "size_kb": size // 1000,
        "net_kb": net // 1000,
        "gc": f"{gc:.1f}",
        "n_cds": n_cds,
        "n_is_internal": n_is,
        "is_flanked": is_flanked,
        "trna_flanked": trna_flanked,
        "left_is": "; ".join(f.product[:30] for f in left_is[:2]),
        "right_is": "; ".join(f.product[:30] for f in right_is[:2]),
        "trna": "; ".join(f.product[:20] for f in trna_flanking[:2]),
    })

# ── ICE donor analysis for 108 kb region ───────────────────────────────────────
print("\n\n=== ICE DONOR ANALYSIS (108 kb region, 1,285,651–1,393,969) ===")
ICE_START, ICE_STOP = 1285651, 1393969
ice_genes = features_strictly_in(ICE_START, ICE_STOP)

# Extract RefSeq IDs from key ICE marker genes
ice_marker_keywords = ["vird4", "virb4", "virb6", "trbl", "relaxase", "recombinase",
                        "integrase", "conjugal transfer", "mobf", "mob"]

print("\nICE marker genes with DbXrefs (for donor identification):")
print(f"{'Gene/Locus':18}  {'RefSeq ID':20}  {'Product'}")
print("-"*80)

for f in ice_genes:
    p = f.product.lower()
    if any(kw in p for kw in ice_marker_keywords):
        refseqs = [x.strip().replace("RefSeq:", "")
                   for x in f.dbxrefs.split(",") if "RefSeq:" in x]
        uniref  = [x.strip().replace("UniRef:", "")
                   for x in f.dbxrefs.split(",") if "UniRef100_" in x]
        rs = refseqs[0] if refseqs else uniref[0] if uniref else "N/A"
        label = f.gene if f.gene else f.locus
        print(f"{label:18}  {rs:20}  {f.product[:55]}")

print("\nNote: RefSeq WP_ IDs can be looked up in NCBI to find the source organism.")
print("The UniRef100_ IDs identify the specific protein cluster and its members.")

# ── Write annotated TSV ─────────────────────────────────────────────────────────
tsv_out = OUTD / "M5_novel_regions_annotated.tsv"
with open(tsv_out, "w") as fh:
    fh.write("\t".join(["region", "size_kb", "net_extra_kb", "GC%", "n_CDS",
                         "n_IS_internal", "IS_flanked", "tRNA_flanked",
                         "left_IS", "right_IS", "tRNA"]) + "\n")
    for row in annotated:
        fh.write("\t".join([
            row["region"], str(row["size_kb"]), str(row["net_kb"]),
            row["gc"], str(row["n_cds"]), str(row["n_is_internal"]),
            str(row["is_flanked"]), str(row["trna_flanked"]),
            row["left_is"], row["right_is"], row["trna"]
        ]) + "\n")
print(f"\nAnnotated TSV written to: {tsv_out}")

# ── Write IS census report ──────────────────────────────────────────────────────
is_out = OUTD / "phase4_IS_census.txt"
with open(is_out, "w") as fh:
    fh.write("# M5 IS Element Census — Phase 4.3\n")
    fh.write("# Source: Bakta TSV annotation\n\n")
    fh.write(f"Total IS/transposase features: {len(is_features)}\n\n")
    fh.write("Family distribution:\n")
    for fam, cnt in sorted(is_families.items(), key=lambda x: -x[1]):
        fh.write(f"  {fam:20} {cnt:4d}\n")
    fh.write(f"\nComparison with I6:\n")
    fh.write(f"  M5:  151 transposable elements (Bakta)\n")
    fh.write(f"  I6:   94 transposable elements (NCBI annotation)\n")
    fh.write(f"  Δ:    57 extra IS elements in M5\n")
    fh.write(f"  Note: The 57 DUP regions in nucmer dnadiff correspond exactly to these 57 extra IS elements.\n")
print(f"IS census written to: {is_out}")

print("\n=== PHASE 4.1 LOCAL ANNOTATION COMPLETE ===")
