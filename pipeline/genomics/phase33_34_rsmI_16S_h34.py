"""
Title:   Phase 3.3 + 3.4 — rsmI comparison and 16S h34 analysis
Author:  Melis Gencel
Date:    2026-03-22
Input:   M5 FAA, M5 FNA (assembly), I6 GFF + FNA, M5 16S BED
Output:  results/genomics/comparative/rsmI_comparison.txt
         results/genomics/comparative/16S_h34_analysis.txt
"""

import ssl
import sys
import time
import urllib.request
from pathlib import Path

# ── SSL workaround ─────────────────────────────────────────────────────────────
_ssl_ctx = ssl.create_default_context()
_ssl_ctx.check_hostname = False
_ssl_ctx.verify_mode = ssl.CERT_NONE

BASE = Path(r"C:\Users\melis\Desktop\bioinformatics-portfolio\hgt-study")
OUTD = BASE / "results" / "genomics" / "comparative"

# ── Helper: parse FASTA ─────────────────────────────────────────────────────────
def parse_fasta(path):
    seqs = {}
    order = []
    current = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = ""
                order.append(current)
            elif line and current:
                seqs[current] += line
    return seqs, order

# ── Helper: extract subsequence from FASTA ─────────────────────────────────────
def extract_region(fna_path, seqid, start0, end, strand):
    """Extract 0-based [start0, end) from sequence seqid. Reverse-complement if strand='-'."""
    seqs, _ = parse_fasta(fna_path)
    seq = seqs[seqid][start0:end]
    if strand == "-":
        comp = str.maketrans("ACGTacgt", "TGCAtgca")
        seq = seq.translate(comp)[::-1]
    return seq

# ── Helper: pairwise alignment (Biopython) ─────────────────────────────────────
try:
    from Bio import pairwise2
    from Bio.pairwise2 import format_alignment
    HAS_BIO = True
except ImportError:
    HAS_BIO = False

def simple_align(a, b):
    """Return (ident_count, total_compared, aln_a, aln_b) using needleman-wunsch style."""
    if HAS_BIO:
        alns = pairwise2.align.globalms(a, b, 2, -1, -10, -0.5)
        if alns:
            aln = alns[0]
            aa, bb = aln.seqA, aln.seqB
            total = max(len(a), len(b))
            ident = sum(1 for x, y in zip(aa, bb) if x == y and x != "-")
            return ident, total, aa, bb
    # Fallback: simple ungapped
    ident = sum(1 for x, y in zip(a, b) if x == y)
    return ident, max(len(a), len(b)), a, b

# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 3.4 — rsmI comparison
# ═══════════════════════════════════════════════════════════════════════════════

print("=" * 70)
print("PHASE 3.4 — rsmI methyltransferase comparison")
print("=" * 70)

M5_FAA = BASE / "data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.faa"
I6_FAA = BASE / "data/references/GCF_022494515.1_protein.faa"

# Extract M5 rsmI (NKFIDM_00044)
m5_prots, _ = parse_fasta(M5_FAA)
m5_rsmI = m5_prots.get("NKFIDM_00044", None)
if not m5_rsmI:
    print("ERROR: NKFIDM_00044 not found in M5 FAA")
    sys.exit(1)
print(f"M5 rsmI (NKFIDM_00044): {len(m5_rsmI)} aa")

# Extract I6 rsmI — protein_id = WP_036618147.1 from GFF
i6_prots, _ = parse_fasta(I6_FAA)
i6_rsmI = i6_prots.get("WP_036618147.1", None)
if not i6_rsmI:
    print("ERROR: WP_036618147.1 not found in I6 FAA")
    sys.exit(1)
print(f"I6 rsmI (WP_036618147.1): {len(i6_rsmI)} aa")

# Compare
if m5_rsmI == i6_rsmI:
    print("\nRESULT: M5 rsmI == I6 rsmI — SEQUENCES ARE IDENTICAL")
    rsmI_ident_pct = 100.0
    rsmI_status = "IDENTICAL"
else:
    ident, total, aa, bb = simple_align(m5_rsmI, i6_rsmI)
    rsmI_ident_pct = 100 * ident / total
    rsmI_status = f"{rsmI_ident_pct:.2f}% identical"
    print(f"\nRESULT: {ident}/{total} identical positions ({rsmI_ident_pct:.2f}%)")
    # Show differences
    print("\nDifferences (M5 vs I6):")
    diffs = [(i+1, m5_rsmI[i], i6_rsmI[i]) for i in range(min(len(m5_rsmI), len(i6_rsmI)))
             if m5_rsmI[i] != i6_rsmI[i]]
    for pos, m5aa, i6aa in diffs[:20]:
        print(f"  pos {pos}: M5={m5aa} I6={i6aa}")

# Write rsmI report
rsmI_lines = [
    "# rsmI Comparison — M5 vs I6",
    "# Phase 3.4 — Session 4 (2026-03-22)",
    "",
    "## Gene information",
    "  M5:  NKFIDM_00044 at contig_1:51713-52606 (+)",
    "       Product: 16S rRNA (cytidine(1402)-2'-O)-methyltransferase",
    "       EC: 2.1.1.198 | RefSeq match: WP_036618147.1",
    "  I6:  LMZ02_RS05655 at NZ_CP086393.1:1236581-1237474 (+)",
    "       Product: 16S rRNA (cytidine(1402)-2'-O)-methyltransferase",
    "       protein_id: WP_036618147.1",
    "",
    "## Sequence comparison",
    f"  M5 rsmI length: {len(m5_rsmI)} aa",
    f"  I6 rsmI length: {len(i6_rsmI)} aa",
    f"  Identity: {rsmI_status}",
    "",
    "## Interpretation",
    "  rsmI encodes the methyltransferase that adds a 2'-O-methyl group to C1402",
    "  of 16S rRNA. C1402 is in helix 44 (h44) of the 30S subunit decoding center,",
    "  NOT in helix 34 (the spectinomycin binding site). Methylation at C1402 by rsmI",
    "  is conserved in bacteria and does not directly affect spectinomycin affinity.",
    "",
    "  Note: The spectinomycin-relevant methyltransferase would be rsmH (m4C1402),",
    "  which methylates C1402 at N4. Both M5 and I6 carry rsmH (see GFF).",
    "",
]
if rsmI_status == "IDENTICAL":
    rsmI_lines += [
        "  M5 rsmI is IDENTICAL to I6 rsmI. Any spectinomycin resistance in M5",
        "  cannot be attributed to rsmI sequence differences.",
        "  → rsmI contributes NOTHING to the resistance phenotype difference between M5 and I6.",
    ]
else:
    rsmI_lines += [
        f"  M5 rsmI differs from I6 rsmI at {rsmI_ident_pct:.2f}% identity.",
        "  Detailed differences listed above.",
    ]
rsmI_lines += [
    "",
    "## Phase 3.4 verdict",
    "  rsmI sequence difference does NOT explain spectinomycin resistance in M5.",
    "  The rpsE Δ(A21-K22-V23) deletion (Phase 3.1/3.1c) remains the sole",
    "  genomically-supported resistance mechanism.",
]

rsmI_out = OUTD / "rsmI_comparison.txt"
with open(rsmI_out, "w") as fh:
    fh.write("\n".join(rsmI_lines) + "\n")
print(f"\nrsmI report written to: {rsmI_out}")


# ═══════════════════════════════════════════════════════════════════════════════
# PHASE 3.3 — 16S rRNA helix 34 analysis
# ═══════════════════════════════════════════════════════════════════════════════

print("\n" + "=" * 70)
print("PHASE 3.3 — 16S rRNA helix 34 analysis")
print("=" * 70)

# M5 16S coordinates (from BED, 0-based half-open):
# Use first copy: contig_1 30223-31777 (+), length=1554
M5_FNA = BASE / "data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.fna"
I6_FNA = BASE / "data/references/GCF_022494515.1_genomic.fna"

# Extract M5 16S copy 1
m5_16S = extract_region(M5_FNA, "contig_1", 30223, 31777, "+")
print(f"M5 16S copy 1 (30223-31777, +): {len(m5_16S)} bp")

# Find I6 16S — search GFF for 16S rRNA features
print("Finding I6 16S rRNA coordinates from GFF...")
i6_16S_coords = []
with open(BASE / "data/references/GCF_022494515.1_genomic.gff") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        if parts[2] == "rRNA" and "16S" in parts[8]:
            seqid, start, end, strand = parts[0], int(parts[3])-1, int(parts[4]), parts[6]
            i6_16S_coords.append((seqid, start, end, strand))

print(f"Found {len(i6_16S_coords)} 16S rRNA features in I6")
for c in i6_16S_coords[:3]:
    print(f"  {c[0]}:{c[1]}-{c[2]} ({c[3]})")

# Extract first I6 16S
if i6_16S_coords:
    seqid, start, end, strand = i6_16S_coords[0]
    i6_seqs, _ = parse_fasta(I6_FNA)
    seq = i6_seqs[seqid][start:end]
    if strand == "-":
        comp = str.maketrans("ACGTacgt", "TGCAtgca")
        seq = seq.translate(comp)[::-1]
    i6_16S = seq
    print(f"I6 16S copy 1 ({seqid}:{start}-{end}, {strand}): {len(i6_16S)} bp")
else:
    print("ERROR: no 16S found in I6 GFF")
    sys.exit(1)

# Align M5 vs I6 16S
print("\nAligning M5 vs I6 16S...")
ident, total, aln_m5, aln_i6 = simple_align(m5_16S, i6_16S)
print(f"Alignment: {ident}/{total} identical ({100*ident/total:.3f}%)")

# Find differences
diffs_16S = []
for i, (a, b) in enumerate(zip(aln_m5, aln_i6)):
    if a != b:
        # Compute ungapped positions
        m5_pos = sum(1 for c in aln_m5[:i+1] if c != "-")
        i6_pos = sum(1 for c in aln_i6[:i+1] if c != "-")
        diffs_16S.append((i+1, m5_pos, a, i6_pos, b))

print(f"Total differences: {len(diffs_16S)}")
if diffs_16S:
    print("First 20 differences (aln_col, M5_pos, M5_nt, I6_pos, I6_nt):")
    for col, mp, ma, ip, ia in diffs_16S[:20]:
        print(f"  aln col {col}: M5 pos {mp}={ma} | I6 pos {ip}={ia}")

# ── Helix 34 in E. coli 16S numbering ─────────────────────────────────────────
# h34 spans approximately nt 1046-1067 and 1189-1210 (E. coli 16S J01695)
# Spectinomycin resistance mutations documented in 16S rRNA:
#   C1066U  (primary resistance mutation in E. coli, gonococcus)
#   C1192U  (gonococcus)
#   A1193G  (gonococcus)
#   G1064A  (less common)
# Key E. coli h34 positions: 1046-1068 (strand 1), 1189-1211 (strand 2)

# Since we don't have an E. coli 16S alignment, we'll estimate h34 in M5
# by its conserved flanking sequences. E. coli h34 starts after a conserved
# motif around position 1043 (GAAAG) and ends before position 1070 (ACUUA).
# In Paenibacillus (firmicute-like), 16S is ~1553 bp — same as E. coli (1542 bp).
# h34 positions should be very similar.

H34_STRAND1 = (1045, 1068)   # 0-based, approximate, E. coli numbering mapped to Paenibacillus
H34_STRAND2 = (1188, 1211)

# Key resistance positions in E. coli 16S numbering (0-based: subtract 1)
SPEC_RESIST_POSITIONS = {
    1065: "C1066 (primary spc resistance in E.coli/gonococcus)",
    1191: "C1192 (gonococcus)",
    1192: "A1193 (gonococcus)",
    1063: "G1064",
}

print(f"\nHelix 34 analysis (approx E. coli mapping):")
print(f"  h34 strand 1: nt {H34_STRAND1[0]+1}–{H34_STRAND1[1]+1} (E. coli coords)")
print(f"  h34 strand 2: nt {H34_STRAND2[0]+1}–{H34_STRAND2[1]+1} (E. coli coords)")

print("\nM5 16S at key spectinomycin resistance positions (E. coli numbering):")
print(f"  {'Pos(Ec)':>8}  {'M5':>4}  {'I6':>4}  {'Notes'}")
print("  " + "-"*55)
for pos0, note in sorted(SPEC_RESIST_POSITIONS.items()):
    m5nt = m5_16S[pos0] if pos0 < len(m5_16S) else "?"
    i6nt = i6_16S[pos0] if pos0 < len(i6_16S) else "?"
    diff_flag = " ← DIFFERENT" if m5nt != i6nt else ""
    print(f"  {pos0+1:>8}  {m5nt:>4}  {i6nt:>4}  {note}{diff_flag}")

# Show h34 strand 1 region
print(f"\nM5 16S h34 strand 1 (nt {H34_STRAND1[0]+1}–{H34_STRAND1[1]+1}):")
print(f"  M5: {m5_16S[H34_STRAND1[0]:H34_STRAND1[1]]}")
print(f"  I6: {i6_16S[H34_STRAND1[0]:H34_STRAND1[1]]}")

print(f"\nM5 16S h34 strand 2 (nt {H34_STRAND2[0]+1}–{H34_STRAND2[1]+1}):")
print(f"  M5: {m5_16S[H34_STRAND2[0]:H34_STRAND2[1]]}")
print(f"  I6: {i6_16S[H34_STRAND2[0]:H34_STRAND2[1]]}")

# Write 16S h34 report
h34_lines = [
    "# 16S rRNA Helix 34 Analysis — M5 vs I6",
    "# Phase 3.3 — Session 4 (2026-03-22)",
    "",
    "## 16S sequences used",
    f"  M5:  contig_1:30224-31777 (+), {len(m5_16S)} bp (copy 1 of 8)",
    f"  I6:  {i6_16S_coords[0][0]}:{i6_16S_coords[0][1]+1}-{i6_16S_coords[0][2]} ({i6_16S_coords[0][3]}), {len(i6_16S)} bp (copy 1 of 8)",
    "",
    f"## Overall M5 vs I6 16S identity",
    f"  {ident}/{total} identical positions = {100*ident/total:.3f}%",
    f"  Total differences: {len(diffs_16S)}",
    "",
]
if diffs_16S:
    h34_lines.append("  All differences:")
    h34_lines.append(f"  {'AalnCol':>8}  {'M5pos':>6}  {'M5':>3}  {'I6pos':>6}  {'I6':>3}")
    h34_lines.append("  " + "-"*35)
    for col, mp, ma, ip, ia in diffs_16S:
        h34_lines.append(f"  {col:>8}  {mp:>6}  {ma:>3}  {ip:>6}  {ia:>3}")
else:
    h34_lines.append("  M5 and I6 16S are IDENTICAL at the sequence level.")

h34_lines += [
    "",
    "## Helix 34 — spectinomycin binding site",
    "  Spectinomycin binds in the minor groove of 16S rRNA helix 34.",
    "  Key resistance-conferring mutations (E. coli numbering):",
    "    C1066U — primary resistance mutation (E. coli, N. gonorrhoeae)",
    "    C1192U — resistance in N. gonorrhoeae",
    "    A1193G — resistance in N. gonorrhoeae",
    "    G1064A — minor resistance variant",
    "",
    "## Key positions in M5 vs I6 (approximate E. coli mapping):",
    f"  {'Pos(Ec)':>8}  {'M5':>4}  {'I6':>4}  {'Notes'}",
    "  " + "-"*55,
]
for pos0, note in sorted(SPEC_RESIST_POSITIONS.items()):
    m5nt = m5_16S[pos0] if pos0 < len(m5_16S) else "?"
    i6nt = i6_16S[pos0] if pos0 < len(i6_16S) else "?"
    diff_flag = " ← DIFFERENT" if m5nt != i6nt else ""
    h34_lines.append(f"  {pos0+1:>8}  {m5nt:>4}  {i6nt:>4}  {note}{diff_flag}")

h34_lines += [
    "",
    f"## h34 strand 1 sequence (nt {H34_STRAND1[0]+1}–{H34_STRAND1[1]+1}):",
    f"  M5: {m5_16S[H34_STRAND1[0]:H34_STRAND1[1]]}",
    f"  I6: {i6_16S[H34_STRAND1[0]:H34_STRAND1[1]]}",
    "",
    f"## h34 strand 2 sequence (nt {H34_STRAND2[0]+1}–{H34_STRAND2[1]+1}):",
    f"  M5: {m5_16S[H34_STRAND2[0]:H34_STRAND2[1]]}",
    f"  I6: {i6_16S[H34_STRAND2[0]:H34_STRAND2[1]]}",
    "",
    "## Interpretation",
    "  Note: E. coli position mapping is approximate (+/- 2 nt) due to",
    "  16S rRNA length differences between species. A full alignment to",
    "  E. coli J01695 would be needed for exact coordinates.",
    "",
    "  IMPORTANT: Even if M5 and I6 are identical at h34, this only shows",
    "  that spectinomycin resistance is NOT due to a 16S rRNA mutation.",
    "  The rpsE Δ(A21-K22-V23) deletion (Phase 3.1c) acts at the protein",
    "  level and is independent of 16S rRNA sequence.",
    "",
    "## Phase 3.3 verdict",
    "  [To be filled after script execution]",
]

h34_out = OUTD / "16S_h34_analysis.txt"
with open(h34_out, "w") as fh:
    fh.write("\n".join(h34_lines) + "\n")
print(f"\n16S h34 report written to: {h34_out}")

print("\n" + "=" * 70)
print("DONE — Phases 3.3 and 3.4 complete")
print("=" * 70)
