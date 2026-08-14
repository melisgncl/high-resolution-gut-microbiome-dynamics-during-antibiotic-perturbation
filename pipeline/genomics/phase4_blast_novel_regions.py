"""
Title:   Phase 4 — BLAST representative windows from unexplored novel regions
Author:  Melis Gencel
Date:    2026-03-23
Purpose: NCBI BLAST top unexplored novel regions to identify their origin.
         The 108 kb region is already characterized (ICE). Here we BLAST the
         59 kb, 48 kb, 25 kb, and 22 kb blocks.
Input:   M5 FNA assembly
Output:  results/genomics/comparative/phase4_blast_results.txt
"""

import ssl
import time
from pathlib import Path
from io import StringIO

BASE = Path(r"C:\Users\melis\Desktop\bioinformatics-portfolio\hgt-study")
OUTD = BASE / "results" / "genomics" / "comparative"

_ssl_ctx = ssl.create_default_context()
_ssl_ctx.check_hostname = False
_ssl_ctx.verify_mode = ssl.CERT_NONE

# Load M5 assembly
m5_seq = ""
FNA = BASE / "data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.fna"
with open(FNA) as fh:
    for line in fh:
        if not line.startswith(">"):
            m5_seq += line.strip()
print(f"Assembly loaded: {len(m5_seq):,} bp")

# Regions to BLAST (top unexplored — skip 108 kb ICE which is already characterized)
# Format: (label, start, stop, note)
REGIONS = [
    ("59kb_6766kb",  6766334, 6825514, "59 kb block at 6.77 Mb"),
    ("48kb_5322kb",  5322917, 5371421, "48 kb block at 5.32 Mb"),
    ("25kb_2911kb",  2911705, 2936697, "25 kb block at 2.91 Mb"),
    ("22kb_5297kb",  5297204, 5322425, "22 kb block at 5.30 Mb"),
    ("18kb_2547kb",  2546997, 2565301, "18 kb block at 2.55 Mb"),
    # Also include 108 kb ICE for donor identification
    ("108kb_ICE",    1285651, 1393969, "108 kb ICE — identify donor organism"),
]

def extract_window(seq, start, stop, window=4000):
    """Extract a window centred in the region, max `window` bp."""
    mid = (start + stop) // 2
    w2  = window // 2
    s   = max(0, mid - w2)
    e   = min(len(seq), s + window)
    return seq[s:e], s, e

# Run BLAST via Biopython NCBI qblast
try:
    from Bio.Blast import NCBIWWW, NCBIXML
    HAS_BIO = True
except ImportError:
    HAS_BIO = False
    print("ERROR: Biopython not available")

report_lines = [
    "# Phase 4 — BLAST Results for M5 Novel Regions",
    "# Date: 2026-03-23",
    "# Tool: NCBI qblast (blastn vs nt)",
    "# Window: 4 kb centred in each novel region",
    "",
]

if HAS_BIO:
    for label, start, stop, note in REGIONS:
        window_seq, ws, we = extract_window(m5_seq, start, stop)
        print(f"\nBLASTing {label} ({note})")
        print(f"  Window: {ws:,}–{we:,} ({len(window_seq)} bp)")
        fasta_query = f">M5_{label}\n{window_seq}\n"

        try:
            print("  Submitting to NCBI qblast...")
            result_handle = NCBIWWW.qblast(
                "blastn", "nt", fasta_query,
                hitlist_size=5,
                expect=1e-10,
                word_size=28,
                perc_ident=70,
            )
            blast_records = list(NCBIXML.parse(result_handle))
            record = blast_records[0] if blast_records else None

            report_lines += [
                f"## {label} — {note}",
                f"   Window: contig_1:{ws+1}-{we} ({len(window_seq)} bp centred in {start}-{stop})",
                "",
            ]

            if record and record.alignments:
                print(f"  Hits: {len(record.alignments)}")
                report_lines.append(f"   Top hits (blastn vs nt, E<1e-10):")
                report_lines.append(f"   {'Rank':4}  {'%Id':5}  {'Cov':5}  {'E-val':10}  Organism/Description")
                report_lines.append("   " + "-"*75)
                for i, aln in enumerate(record.alignments[:5]):
                    hsp = aln.hsps[0]
                    pct_id  = 100 * hsp.identities / hsp.align_length
                    pct_cov = 100 * hsp.align_length / len(window_seq)
                    desc = aln.title[:80].replace("\n", " ")
                    eval_str = f"{hsp.expect:.1e}"
                    report_lines.append(f"   {i+1:4}  {pct_id:5.1f}  {pct_cov:5.1f}  {eval_str:10}  {desc}")
                    print(f"  Hit {i+1}: {pct_id:.1f}% id, {eval_str} — {desc[:60]}")
            else:
                report_lines.append("   No significant hits (blastn vs nt, E<1e-10)")
                print("  No hits")

            report_lines.append("")

        except Exception as e:
            msg = f"  BLAST error for {label}: {e}"
            print(msg)
            report_lines += [
                f"## {label} — {note}",
                f"   BLAST FAILED: {e}",
                "",
            ]

        # Rate limiting: NCBI allows 3 requests/second with API key, 1/sec without
        print("  Waiting 5 s before next query...")
        time.sleep(5)

out = OUTD / "phase4_blast_results.txt"
with open(out, "w") as fh:
    fh.write("\n".join(report_lines) + "\n")
print(f"\nBLAST report written to: {out}")
