"""
Title:   Build multi-FASTA input for rpsE MAFFT MSA
Author:  Melis Gencel
Date:    2026-03-22
Input:   NCBI E-utils (E. coli NP_417762.1, B. subtilis rpsE)
         results/genomics/comparative/rpsE_isolates_vs_I6.faa (M5 + I6)
Output:  results/genomics/comparative/rpsE_msa_input.faa
         (5 sequences: M5, I6, E. coli K12, B. subtilis, H. influenzae)
"""

import ssl
import json
import time
from urllib import request

_ssl_ctx = ssl.create_default_context()
_ssl_ctx.check_hostname = False
_ssl_ctx.verify_mode = ssl.CERT_NONE

BASE = r"C:\Users\melis\Desktop\bioinformatics-portfolio\hgt-study"


def ncbi_search_rpsE(organism, retmax=5):
    """Return list of protein UIDs for rpsE in given organism (RefSeq)."""
    q = f"rpsE[Gene+Name]+AND+{organism.replace(' ', '+')}[Organism]+AND+refseq[Filter]"
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=protein&term={q}&retmode=json&retmax={retmax}"
    )
    with request.urlopen(url, timeout=30, context=_ssl_ctx) as r:
        return json.loads(r.read().decode())["esearchresult"]["idlist"]


def ncbi_fetch_fasta(uid):
    """Fetch FASTA for a protein UID. Respects NCBI rate limit (3/sec)."""
    time.sleep(0.4)
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=protein&id={uid}&rettype=fasta&retmode=text"
    )
    with request.urlopen(url, timeout=30, context=_ssl_ctx) as r:
        return r.read().decode().strip()


def ncbi_search_rpsE_with_delay(organism, retmax=5):
    """Search + delay."""
    time.sleep(0.4)
    return ncbi_search_rpsE(organism, retmax)


def pick_rpsE(ids, size_range=(140, 200), keywords=("S5", "rpsE", "ribosomal")):
    """Fetch and return the first hit matching expected rpsE size and keywords."""
    for uid in ids:
        fasta = ncbi_fetch_fasta(uid)
        lines = fasta.splitlines()
        hdr = lines[0]
        seq = "".join(lines[1:])
        print(f"  ID {uid}: {len(seq)} aa | {hdr[:80]}")
        if size_range[0] < len(seq) < size_range[1] and any(
            kw.lower() in hdr.lower() for kw in keywords
        ):
            print(f"  --> Selected")
            return hdr, seq
    raise RuntimeError(f"No suitable rpsE found among IDs: {ids}")


# ── Load Paenibacillus sequences from existing file ────────────────────────────
print("Loading M5 and I6 from rpsE_isolates_vs_I6.faa ...")
m5_hdr = m5_seq = i6_hdr = i6_seq = None
with open(rf"{BASE}\results\genomics\comparative\rpsE_isolates_vs_I6.faa") as fh:
    current_hdr = None
    current_seq = []
    for line in fh:
        line = line.rstrip()
        if line.startswith(">"):
            if current_hdr and current_seq:
                seq = "".join(current_seq)
                if "M5" in current_hdr:
                    m5_hdr, m5_seq = current_hdr, seq
                elif "I6" in current_hdr:
                    i6_hdr, i6_seq = current_hdr, seq
            current_hdr = line
            current_seq = []
        else:
            current_seq.append(line)
    # last record
    if current_hdr and current_seq:
        seq = "".join(current_seq)
        if "M5" in current_hdr:
            m5_hdr, m5_seq = current_hdr, seq
        elif "I6" in current_hdr:
            i6_hdr, i6_seq = current_hdr, seq

print(f"  M5: {len(m5_seq)} aa")
print(f"  I6: {len(i6_seq)} aa")

# ── Fetch E. coli K12 rpsE ────────────────────────────────────────────────────
print("\nFetching E. coli K12 rpsE ...")
ec_ids = ncbi_search_rpsE_with_delay("Escherichia coli K-12")
ec_hdr, ec_seq = pick_rpsE(ec_ids)
print(f"  E. coli: {len(ec_seq)} aa")

# ── Fetch B. subtilis rpsE ────────────────────────────────────────────────────
print("\nFetching B. subtilis rpsE ...")
bs_ids = ncbi_search_rpsE_with_delay("Bacillus subtilis 168")
bs_hdr, bs_seq = pick_rpsE(bs_ids)
print(f"  B. subtilis: {len(bs_seq)} aa")

# ── Write multi-FASTA ─────────────────────────────────────────────────────────
out_path = rf"{BASE}\results\genomics\comparative\rpsE_msa_input.faa"

def fmt_fasta(header, seq, width=60):
    lines = [header]
    for i in range(0, len(seq), width):
        lines.append(seq[i:i+width])
    return "\n".join(lines)

# Clean headers for MAFFT (short, no special chars)
records = [
    (">Paeni_M5_isolate",       m5_seq),
    (">Paeni_I6_reference",     i6_seq),
    (">Ecoli_K12_NP417762",     ec_seq),
    (">Bsubtilis_168",          bs_seq),
]

with open(out_path, "w") as fh:
    for hdr, seq in records:
        fh.write(fmt_fasta(hdr, seq) + "\n\n")

print(f"\nMulti-FASTA written: {out_path}")
print("Sequences:")
for hdr, seq in records:
    print(f"  {hdr}: {len(seq)} aa")
