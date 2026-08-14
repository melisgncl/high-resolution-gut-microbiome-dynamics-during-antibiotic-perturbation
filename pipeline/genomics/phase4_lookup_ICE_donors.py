"""
Title:   Phase 4 — Look up ICE marker gene WP_ accessions to identify donor organisms
Author:  Melis Gencel
Date:    2026-03-23
Method:  NCBI Entrez efetch on WP_ protein IDs → get source organism
Output:  results/genomics/comparative/phase4_ICE_donor_lookup.txt
"""

import ssl, time, urllib.request, json
from pathlib import Path

BASE = Path(r"C:\Users\melis\Desktop\bioinformatics-portfolio\hgt-study")
OUTD = BASE / "results" / "genomics" / "comparative"

_ssl_ctx = ssl.create_default_context()
_ssl_ctx.check_hostname = False
_ssl_ctx.verify_mode = ssl.CERT_NONE

def entrez_fetch_organism(accession):
    """Fetch protein summary from NCBI eutils to get source organism."""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    # Use esummary for protein accession
    url = f"{base_url}esummary.fcgi?db=protein&id={accession}&retmode=json"
    try:
        with urllib.request.urlopen(url, context=_ssl_ctx, timeout=15) as r:
            data = json.loads(r.read().decode())
        time.sleep(0.4)  # rate limit
        uids = list(data["result"].get("uids", []))
        if not uids:
            return None, None, None
        uid = uids[0]
        doc = data["result"].get(uid, {})
        title = doc.get("title", "N/A")
        organism = doc.get("organism", "N/A")
        taxid = doc.get("taxid", "N/A")
        return title, organism, taxid
    except Exception as e:
        return None, None, str(e)

# ICE marker genes from 108 kb region
ICE_WP_IDS = [
    ("WP_216539118.1",  "Recombinase",              "NKFIDM_01129"),
    ("WP_255222961.1",  "Mobilization protein",      "NKFIDM_01139"),
    ("WP_036623496.1",  "Relaxase",                  "NKFIDM_01140"),
    ("WP_036623490.1",  "VirD4 coupling protein",    "NKFIDM_01144"),
    ("WP_036623477.1",  "VirB4 ATPase",              "NKFIDM_01151"),
    ("WP_036628241.1",  "Conjugal transfer protein", "NKFIDM_01200"),
    ("WP_227872718.1",  "Conj. transposon recombinase", "NKFIDM_01201"),
]

# Second ICE (25 kb at 2.91 Mb) — get accessions from Bakta TSV
SECOND_ICE_WP_IDS = []  # will populate from TSV

# Read Bakta TSV for second ICE region (2,911,705–2,936,697)
M5_TSV = BASE / "data/sequence_data/TPHP3P_results/TPHP3P_1_M5__+_/annotation/TPHP3P_1_M5__+_.tsv"
with open(M5_TSV) as fh:
    for line in fh:
        if line.startswith("#"): continue
        p = line.strip().split("\t")
        if len(p) < 9: continue
        try:
            start, stop = int(p[2]), int(p[3])
        except ValueError: continue
        if 2911705 <= start <= 2936697:
            dbxrefs = p[8]
            product = p[7] if len(p) > 7 else ""
            locus   = p[5]
            for x in dbxrefs.split(","):
                x = x.strip()
                if x.startswith("RefSeq:WP_"):
                    wp = x.replace("RefSeq:", "")
                    SECOND_ICE_WP_IDS.append((wp, product[:50], locus))
                    break

print("=== ICE DONOR IDENTIFICATION VIA ENTREZ ===\n")
lines = [
    "# ICE Donor Identification via NCBI Entrez",
    "# Phase 4 — Session 4 (2026-03-23)",
    "# Method: NCBI esummary on WP_ protein accessions → source organism",
    "",
    "## 108 kb ICE (1,285,651–1,393,969) — key marker genes",
    "",
    f"  {'Locus':16}  {'WP Accession':20}  {'Organism':35}  Product",
    "  " + "-"*100,
]

for wp, product, locus in ICE_WP_IDS:
    title, organism, taxid = entrez_fetch_organism(wp)
    org_str = organism if organism else "N/A"
    print(f"  {locus:16}  {wp:20}  {org_str:35}  {product}")
    lines.append(f"  {locus:16}  {wp:20}  {org_str:35}  {product}")

lines += [
    "",
    "## 25 kb second ICE (2,911,705–2,936,697) — conjugation machinery",
    "",
    f"  {'Locus':16}  {'WP Accession':20}  {'Organism':35}  Product",
    "  " + "-"*100,
]

print("\n=== SECOND ICE (25 kb at 2.91 Mb) ===")
for wp, product, locus in SECOND_ICE_WP_IDS[:10]:
    title, organism, taxid = entrez_fetch_organism(wp)
    org_str = organism if organism else "N/A"
    print(f"  {locus:16}  {wp:20}  {org_str:35}  {product}")
    lines.append(f"  {locus:16}  {wp:20}  {org_str:35}  {product}")

lines += [
    "",
    "## Note on WP_ accessions",
    "  WP_ (non-redundant RefSeq protein) accessions represent conserved proteins",
    "  across multiple genomes. The 'organism' field shows one representative source,",
    "  but the protein may exist in many organisms. For donor identification, the",
    "  organism shown is the most taxonomically representative match in RefSeq.",
    "",
    "## Interpretation",
    "  [Filled from output above]",
]

out = OUTD / "phase4_ICE_donor_lookup.txt"
with open(out, "w") as fh:
    fh.write("\n".join(lines) + "\n")
print(f"\nReport written to: {out}")
