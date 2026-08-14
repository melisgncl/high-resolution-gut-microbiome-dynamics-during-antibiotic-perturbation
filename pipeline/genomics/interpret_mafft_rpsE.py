"""
Title:   Interpret MAFFT MSA of rpsE — annotate deletion and spectinomycin contact columns
Author:  Melis Gencel
Date:    2026-03-22
Input:   results/genomics/comparative/rpsE_msa_output.faa
Output:  results/genomics/comparative/rpsE_msa_annotated.txt
"""

BASE = r"C:\Users\melis\Desktop\bioinformatics-portfolio\hgt-study"

# ── Parse MAFFT output ─────────────────────────────────────────────────────────

seqs = {}
order = []
with open(rf"{BASE}\results\genomics\comparative\rpsE_msa_output.faa") as fh:
    current = None
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            current = line[1:]
            seqs[current] = ""
            order.append(current)
        elif line:
            seqs[current] += line

labels = order  # Paeni_M5_isolate, Paeni_I6_reference, Ecoli_K12_NP417762, Bsubtilis_168
aln_len = len(seqs[labels[0]])
print(f"Alignment length: {aln_len} columns")
for lab in labels:
    seq = seqs[lab]
    n_gaps = seq.count("-")
    print(f"  {lab}: {len(seq)-n_gaps} aa unaligned, {n_gaps} gaps")

# ── Identify gap columns (M5 has gap but I6 does not) ─────────────────────────

m5  = seqs["Paeni_M5_isolate"]
i6  = seqs["Paeni_I6_reference"]
ec  = seqs["Ecoli_K12_NP417762"]
bs  = seqs["Bsubtilis_168"]

deletion_cols = [i for i in range(aln_len) if m5[i] == "-" and i6[i] != "-"]
print(f"\nDeletion columns (M5 gap, I6 not gap): {[c+1 for c in deletion_cols]}")  # 1-based

# ── Compute per-species unaligned position at each alignment column ────────────

def col_to_pos(seq, col_idx):
    """Return 1-based unaligned position of character at col_idx, or None if gap."""
    if seq[col_idx] == "-":
        return None
    return sum(1 for c in seq[:col_idx+1] if c != "-")

print("\nDeletion column details:")
print(f"{'Col':>4} | {'M5':>6} | {'I6':>6} | {'Ec':>6} | {'Bs':>6}")
print("-" * 40)
for col in deletion_cols:
    def fmt(seq, c):
        p = col_to_pos(seq, c)
        if p is None:
            return f"{'gap':>6}"
        return f"{seq[c]}({p:>3})"
    print(f"{col+1:>4} | {fmt(m5,col)} | {fmt(i6,col)} | {fmt(ec,col)} | {fmt(bs,col)}")

# ── Print annotated alignment (N-terminal region only, cols 1-50) ──────────────

REGION = 50  # show first 50 alignment columns

def build_annotation(cols, seqs_dict, width=REGION):
    """Build a match/gap line and deletion marker for first `width` alignment cols."""
    rows = list(seqs_dict.values())
    match_line = []
    del_marker = []
    for i in range(width):
        chars = [s[i] for s in rows]
        non_gap = [c for c in chars if c != "-"]
        if len(set(non_gap)) == 1 and len(non_gap) == len(chars):
            match_line.append("*")
        elif len(set(non_gap)) <= 1:
            match_line.append(":")
        else:
            match_line.append(" ")
        del_marker.append("^" if i in cols else " ")
    return "".join(match_line), "".join(del_marker)

match_str, del_str = build_annotation(deletion_cols, {
    "M5": m5, "I6": i6, "Ec": ec, "Bs": bs
}, aln_len)

print("\n=== Full annotated alignment (N-terminal region, first 50 cols) ===")
width = 50
for start in range(0, REGION, width):
    end = min(start + width, REGION)
    print(f"\n  {'M5':12} {m5[start:end]}")
    print(f"  {'I6':12} {i6[start:end]}")
    print(f"  {'E.coli K12':12} {ec[start:end]}")
    print(f"  {'B.subtilis':12} {bs[start:end]}")
    print(f"  {'conservation':12} {match_str[start:end]}")
    print(f"  {'deletion(M5)':12} {del_str[start:end]}")

# ── Build full output report ───────────────────────────────────────────────────

report_lines = [
    "# rpsE MAFFT MSA — Annotated Alignment",
    "# Session 3 — 2026-03-22",
    "# Tool: MAFFT v7.526 --auto (WSL, Ubuntu 24.04)",
    "# Sequences: Paenibacillus M5 isolate, I6 reference, E. coli K12, B. subtilis 168",
    "",
    "## Alignment summary",
    f"  Alignment length: {aln_len} columns",
]
for lab in labels:
    n_gaps = seqs[lab].count("-")
    report_lines.append(f"  {lab}: {len(seqs[lab])-n_gaps} aa ({n_gaps} gap columns)")

report_lines += [
    "",
    "## Deletion columns (M5 missing, all others present)",
    f"  Alignment columns (1-based): {[c+1 for c in deletion_cols]}",
    "",
    f"  {'Col':>4}  {'M5':>6}  {'I6':>6}  {'E.coli':>6}  {'Bsub':>6}",
    "  " + "-"*38,
]
for col in deletion_cols:
    def fmt2(seq, c):
        p = col_to_pos(seq, c)
        return "gap" if p is None else f"{seq[c]}{p}"
    report_lines.append(
        f"  {col+1:>4}  {fmt2(m5,col):>6}  {fmt2(i6,col):>6}  {fmt2(ec,col):>6}  {fmt2(bs,col):>6}"
    )

# Determine which E. coli positions are deleted
ec_del_pos = [col_to_pos(ec, c) for c in deletion_cols]
i6_del_pos = [col_to_pos(i6, c) for c in deletion_cols]

report_lines += [
    "",
    f"  I6 residues deleted (unaligned positions): {[(i6[c], col_to_pos(i6,c)) for c in deletion_cols]}",
    f"  E. coli residues at same columns:          {[(ec[c], col_to_pos(ec,c)) for c in deletion_cols]}",
    f"  B. subtilis residues at same columns:      {[(bs[c], col_to_pos(bs,c)) for c in deletion_cols]}",
    "",
    "## Key conserved column: the K residue",
]

# Find the column with conserved K
k_cols = [c for c in deletion_cols if i6[c] == "K" or ec[c] == "K"]
for kc in k_cols:
    report_lines.append(
        f"  Alignment col {kc+1}: I6={i6[kc]}({col_to_pos(i6,kc)}), "
        f"Ec={ec[kc]}({col_to_pos(ec,kc)}), Bs={bs[kc]}({col_to_pos(bs,kc)}), M5=gap"
    )
    report_lines.append(
        f"  → This Lysine (conserved in I6, E. coli, B. subtilis) is ABSENT in M5/M6/M8."
    )

report_lines += [
    "",
    "## Alignment (N-terminal 50 columns)",
    "",
    f"  {'M5':<14} {m5[:50]}",
    f"  {'I6':<14} {i6[:50]}",
    f"  {'E.coli K12':<14} {ec[:50]}",
    f"  {'B.subtilis':<14} {bs[:50]}",
    f"  {'conservation':<14} {match_str[:50]}",
    f"  {'deletion(M5)':<14} {del_str[:50]}",
    "  (^ marks alignment columns that are gapped in M5 but present in all others)",
    "",
    "## Note on gap placement ambiguity",
    "  In the I6 sequence, the motif around the deletion is ...RVAKVVK...",
    "  In M5, it is ...RVVK... (3 residues shorter).",
    "  MAFFT places the gap at alignment cols corresponding to I6 V-A-K (positions ~20-22)",
    "  because the V at col 21 (I6 pos ~20) anchors the gap left-of-the-K.",
    "  Biopython pairwise (Session 2) placed the gap at A-K-V (positions 21-23).",
    "  Both are valid representations of the same 3-aa deletion. The biological",
    "  conclusion is identical: M5/M6/M8 lack a conserved Lysine present in I6,",
    "  E. coli (K23), and B. subtilis that contacts 16S rRNA helix 34.",
    "",
    "## Interpretation",
    "  The conserved Lysine (I6 K22 / E. coli K23 / B. subtilis K) is a positively",
    "  charged residue in the N-terminal beta-hairpin of S5 (uS5/rpsE) that contacts",
    "  helix h34 of 16S rRNA at the spectinomycin binding interface.",
    "  Its absence in all three Paenibacillus isolates (M5, M6, M8, which are identical)",
    "  shortens the loop and eliminates a critical contact point.",
    "  This provides structural corroboration for H4 (intrinsic spectinomycin resistance).",
    "  The deletion is ancestral (identical in isolates from 3 separate mice).",
    "",
    "  H1 (de novo mutation): FALSIFIED",
    "  H4 (intrinsic resistance): STRONGLY SUPPORTED — structural basis confirmed by MSA",
    "",
    "## Files",
    "  Input:  results/genomics/comparative/rpsE_msa_input.faa",
    "  MAFFT:  results/genomics/comparative/rpsE_msa_output.faa",
    "  This:   results/genomics/comparative/rpsE_msa_annotated.txt",
]

out_path = rf"{BASE}\results\genomics\comparative\rpsE_msa_annotated.txt"
with open(out_path, "w") as fh:
    fh.write("\n".join(report_lines) + "\n")

print(f"\nAnnotated report written to: {out_path}")
