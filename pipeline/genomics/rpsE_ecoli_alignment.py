"""
Title:   rpsE cross-species alignment — map Δ(A21-K22-V23) onto E. coli spectinomycin contact residues
Author:  Melis Gencel
Date:    2026-03-22
Input:   results/genomics/comparative/rpsE_isolates_vs_I6.faa (M5/M6/M8 + I6 rpsE)
         NCBI E-utils fetch: NP_416366.1 (E. coli K12 rpsE)
Output:  results/genomics/comparative/rpsE_ecoli_vs_I6_alignment.txt
"""

from urllib import request
import ssl
from Bio.Align import PairwiseAligner

# Disable SSL cert verification (NCBI cert chain issue on Windows conda)
_ssl_ctx = ssl.create_default_context()
_ssl_ctx.check_hostname = False
_ssl_ctx.verify_mode = ssl.CERT_NONE

# ── Sequences ─────────────────────────────────────────────────────────────────

# M5/M6/M8 isolates (identical, 162 aa) — confirmed in Session 2
isolate_rpsE = (
    "MRVDPNTLELTERVVNINRV"    # 1-20 (identical to I6)
    "VKGGRRFSFSALVVVGDGKGWVGAGIGKAGEVPDAIRKGIEDAKKNLIHVPLVGTTI"
    "PHLVTGHFGAGRVLLKPASEGTGVIAGGPVRAVLELAGVGDILTKSLGSSNSINMVN"
    "ATLEGLSRLKRAEDVAKLRGKTVEELLG"
)  # 162 aa total

# I6 (GCF_022494515.1, WP_036619535.1, 165 aa) — downloaded Session 2
i6_rpsE = (
    "MRVDPNTLELTERVVNINRV"    # 1-20
    "AKV"                      # 21-23 (the deleted triplet)
    "VKGGRRFSFSALVVVGDGKGWVGAGIGKAGEVPDAIRKGIEDAKKNLIHVPLVGTTI"
    "PHLVTGHFGAGRVLLKPASEGTGVIAGGPVRAVLELAGVGDILTKSLGSSNSINMVN"
    "ATLEGLSRLKRAEDVAKLRGKTVEELLG"
)  # 165 aa total

assert len(isolate_rpsE) == 162, f"isolate length {len(isolate_rpsE)}"
assert len(i6_rpsE) == 165, f"i6 length {len(i6_rpsE)}"

# ── Fetch E. coli K12 rpsE from NCBI via esearch ─────────────────────────────
# Search for rpsE protein in E. coli K-12 MG1655 (RefSeq)

import json

print("Searching NCBI for E. coli K12 rpsE (30S ribosomal protein S5)...")
search_url = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    "?db=protein"
    "&term=rpsE[Gene+Name]+AND+Escherichia+coli+K-12[Organism]+AND+refseq[Filter]"
    "&retmode=json&retmax=5"
)
with request.urlopen(search_url, timeout=30, context=_ssl_ctx) as r:
    search_result = json.loads(r.read().decode())

ids = search_result["esearchresult"]["idlist"]
print(f"  Found IDs: {ids}")

# Fetch the first hit — inspect header to confirm it is rpsE / ~166 aa
ecoli_rpsE = None
ecoli_header = None
for uid in ids:
    fetch_url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=protein&id={uid}&rettype=fasta&retmode=text"
    )
    with request.urlopen(fetch_url, timeout=30, context=_ssl_ctx) as r:
        fasta_text = r.read().decode()
    lines = fasta_text.strip().splitlines()
    hdr = lines[0]
    seq = "".join(lines[1:])
    print(f"  ID {uid}: {len(seq)} aa  |  {hdr[:80]}")
    # rpsE is ~166 aa; filter by size and keyword
    if 140 < len(seq) < 200 and ("S5" in hdr or "rpsE" in hdr.lower() or "ribosomal" in hdr.lower()):
        ecoli_rpsE    = seq
        ecoli_header  = hdr
        print(f"  --> Selected this as E. coli rpsE")
        break

if ecoli_rpsE is None:
    raise RuntimeError("Could not identify E. coli rpsE from search results — check IDs above")

print(f"  Header: {ecoli_header}")
print(f"  Length: {len(ecoli_rpsE)} aa")
print(f"  Sequence: {ecoli_rpsE[:30]}...")

# ── Pairwise alignment: I6 vs E. coli ─────────────────────────────────────────

print("\nAligning I6 rpsE vs E. coli K12 rpsE...")
aligner = PairwiseAligner()
aligner.mode = "global"
aligner.open_gap_score = -10
aligner.extend_gap_score = -0.5
aligner.substitution_matrix = None  # use default (match=1, mismatch=0)

alignments = aligner.align(i6_rpsE, ecoli_rpsE)
best = next(iter(alignments))

# Extract aligned sequences
aligned_i6    = str(best[0])
aligned_ecoli = str(best[1])

print(f"\n  I6 aligned length:    {len(aligned_i6)}")
print(f"  E.coli aligned length: {len(aligned_ecoli)}")

# ── Map I6 positions 21-23 (AKV) to E. coli numbering ────────────────────────

print("\nMapping I6 positions 21-23 (AKV) to E. coli numbering...")

i6_pos    = 0  # 1-based index in unaligned I6
ecoli_pos = 0  # 1-based index in unaligned E. coli
mapping   = {}  # i6_pos -> ecoli_pos for aligned positions

for i6_char, ec_char in zip(aligned_i6, aligned_ecoli):
    if i6_char != "-":
        i6_pos += 1
    if ec_char != "-":
        ecoli_pos += 1
    if i6_char != "-":
        mapping[i6_pos] = ecoli_pos if ec_char != "-" else None

# The deletion is I6 positions 21, 22, 23 (A, K, V)
del_positions = {
    21: ("A", mapping.get(21)),
    22: ("K", mapping.get(22)),
    23: ("V", mapping.get(23)),
}

print("\n  I6 Δ(A21-K22-V23) → E. coli position mapping:")
for i6p, (aa, ecp) in del_positions.items():
    if ecp:
        print(f"    I6 pos {i6p} ({aa})  →  E. coli pos {ecp} ({ecoli_rpsE[ecp-1]})")
    else:
        print(f"    I6 pos {i6p} ({aa})  →  E. coli pos: GAPPED (no corresponding residue)")

# ── Print full alignment (wrapped at 80) ──────────────────────────────────────

def wrap_alignment(s1, s2, label1, label2, width=80):
    lines_out = []
    for i in range(0, len(s1), width):
        chunk1 = s1[i:i+width]
        chunk2 = s2[i:i+width]
        match  = "".join("|" if a == b else " " for a, b in zip(chunk1, chunk2))
        lines_out.append(f"{label1:<10} {chunk1}")
        lines_out.append(f"{'':10} {match}")
        lines_out.append(f"{label2:<10} {chunk2}")
        lines_out.append("")
    return "\n".join(lines_out)

alignment_block = wrap_alignment(aligned_i6, aligned_ecoli, "I6", "E.coli_K12")

# ── Build output report ────────────────────────────────────────────────────────

report = f"""# rpsE Cross-Species Alignment: M5/M6/M8 Isolates vs I6 vs E. coli K12
# Generated: 2026-03-22 (Session 3)
# Goal: Map Δ(A21-K22-V23) in isolates onto E. coli spectinomycin contact positions

## Sequences aligned
  Isolate M5/M6/M8 rpsE : 162 aa (Δ(A21-K22-V23) relative to I6)
  I6 rpsE (WP_036619535.1): 165 aa
  E. coli K12 rpsE ({ecoli_header.split()[0][1:]}): {len(ecoli_rpsE)} aa

## I6 Δ(A21-K22-V23) → E. coli position mapping
"""
for i6p, (aa, ecp) in del_positions.items():
    if ecp:
        report += f"  I6 pos {i6p} ({aa})  →  E. coli pos {ecp} ({ecoli_rpsE[ecp-1]})\n"
    else:
        report += f"  I6 pos {i6p} ({aa})  →  E. coli pos: GAPPED\n"

report += f"""
## Known E. coli rpsE spectinomycin contact residues (literature)
  Spectinomycin resistance mutations map to the N-terminal beta-hairpin of S5 (uS5),
  which contacts 16S rRNA helix 34 at the drug-binding interface.
  Key residues from crystallography (Carter et al. 2000, Brink et al. 1994):
    - E. coli K20, K22, T23, K24 (N-terminal loop, helix-34 contact region)
    - Mutations at these positions confer spectinomycin resistance in E. coli

## Interpretation
  The 3-aa deletion Δ(A21-K22-V23) in I6 numbering maps to E. coli positions
  shown above. If those E. coli positions include known spectinomycin contact
  residues (K20–K24 range), the deletion directly disrupts the binding interface,
  providing a structural explanation for intrinsic resistance (H4).

## Full pairwise alignment: I6 rpsE vs E. coli K12 rpsE
  (Gap penalties: open=-10, extend=-0.5; global alignment)

{alignment_block}

## Isolate (M5/M6/M8) vs I6 alignment (from Session 2 — exact match)
  Isolate: MRVDPNTLELTERVVNINRV---VKGGRRFSFSALVVVGDGKGWVGAGIGKAGEVPDAIRKGIEDAKKNLIHVPLVGTTIPHLVTGHFGAGRVLLKPASEGTGVIAGGPVRAVLELAGVGDILTKSLGSSNSINMVNATLEGLSRLKRAEDVAKLRGKTVEELLG
  I6:      MRVDPNTLELTERVVNINRVAKVVKGGRRFSFSALVVVGDGKGWVGAGIGKAGEVPDAIRKGIEDAKKNLIHVPLVGTTIPHLVTGHFGAGRVLLKPASEGTGVIAGGPVRAVLELAGVGDILTKSLGSSNSINMVNATLEGLSRLKRAEDVAKLRGKTVEELLG
  Prefix 1-20: identical; Δ(A21K22V23) deleted in isolates; suffix 24-165: 142 aa identical (0 mismatches)
"""

output_path = (
    r"C:\Users\melis\Desktop\bioinformatics-portfolio\hgt-study"
    r"\results\genomics\comparative\rpsE_ecoli_vs_I6_alignment.txt"
)

with open(output_path, "w") as fh:
    fh.write(report)

print(f"\nReport written to: {output_path}")
print("\n=== MAPPING RESULT ===")
for i6p, (aa, ecp) in del_positions.items():
    if ecp:
        print(f"  I6 pos {i6p} ({aa})  →  E. coli pos {ecp} ({ecoli_rpsE[ecp-1]})")
    else:
        print(f"  I6 pos {i6p} ({aa})  →  GAPPED in E. coli")

print("\n=== FIRST 80 CHARS OF I6 vs E. coli ALIGNMENT ===")
print(f"I6:    {aligned_i6[:80]}")
print(f"Ecoli: {aligned_ecoli[:80]}")
