"""
Supplementary Figure S1: Hypothesis testing table — Paenibacillus spectinomycin resistance
Author: Melis Gencel
Date: 2026-06-18
Output: results/figures/genomics/figS1_hypothesis_table.png + .pdf
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import textwrap
import os

BASE = r"C:\Users\melis\Desktop\bioinformatics-portfolio\hgt-study"
OUT_DIR = os.path.join(BASE, "results", "figures", "genomics")

# ---------- data ----------
HYPOTHESES = [
    {
        'id': 'H1',
        'hypothesis': 'De novo rpsE mutation\nacquired during the\nexperiment',
        'approach': 'Direct comparison of rpsE\nsequence across M5, M6, M8\nisolates from 3 separate mice',
        'evidence': 'Identical 3-aa in-frame deletion\nΔ(A21-K22-V23) in all 3 isolates.\nProbability of 3 independent\nidentical events: ~10⁻²⁷',
        'verdict': 'FALSIFIED',
    },
    {
        'id': 'H2',
        'hypothesis': 'De novo 16S rRNA\nhelix 34 mutation\n(spectinomycin binding site)',
        'approach': 'Nucleotide comparison of\n16S h34 region (M5 vs I6);\nall 28 16S differences\nmapped to rRNA secondary\nstructure',
        'evidence': 'Helix 34 identical between\nM5 and I6. None of the 28\n16S differences fall in helix 34\nor the spectinomycin binding\nsite (E. coli pos 1063/1064/1066);\ndifferences lie outside h34',
        'verdict': 'FALSIFIED',
    },
    {
        'id': 'H3',
        'hypothesis': 'HGT of spectinomycin\nresistance gene (aadA /\nANT(9)) from E. coli',
        'approach': 'BLASTn (aadA query vs\nM5/M6/M8 genomes);\nAMRFinderPlus, RGI (CARD),\nResFinder on all 3 assemblies',
        'evidence': 'No aadA or ANT(9) detected\nin any isolate. Only amino-\nglycoside genes found (AAC(3),\nANT(6)) target other drugs\nand are ancestral to I6',
        'verdict': 'NOT SUPPORTED',
    },
    {
        'id': 'H4',
        'hypothesis': 'Intrinsic resistance via\npre-existing rpsE deletion\n(strain-level variant)',
        'approach': 'MAFFT MSA (4 species:\nM5/I6/E. coli K-12/\nB. subtilis 168);\nminimap2 + samtools\nread-level validation',
        'evidence': 'Δ(A21-K22-V23) removes\nuniversally conserved Lys\nfrom S5 β-hairpin contacting\n16S h34. Confirmed by 106/184\nreads with 9D261M CIGAR\nat 163× uniform coverage',
        'verdict': 'STRONGLY\nSUPPORTED',
    },
    {
        'id': 'H5',
        'hypothesis': 'rpsE gene duplication\n/ neofunctionalisation\n(one WT + one mutant copy)',
        'approach': 'Bakta annotation copy\nnumber (M5/M6/M8/I6);\nread coverage uniformity\nat rpsE locus',
        'evidence': 'Single rpsE locus in all\n4 genomes (Bakta-confirmed).\n163× coverage = average\ngenome depth; no tandem\nduplication detected',
        'verdict': 'FALSIFIED',
    },
    {
        'id': 'H6',
        'hypothesis': 'Efflux pump upregulation\nor novel efflux gene\nacquisition',
        'approach': 'Pan-genome efflux gene\ncensus; IS-element insertion\nscreening at efflux regulators;\nnovel region annotation',
        'evidence': 'Efflux genes present but not\nspectinomycin-specific (RND,\nMFS, MATE families). No IS\ninsertion into regulators.\nOnly 2 novel efflux genes;\nspectinomycin is a poor\nefflux substrate',
        'verdict': 'NOT SUPPORTED',
    },
    {
        'id': 'H7',
        'hypothesis': 'Novel resistance gene\nnot present in current\nAMR databases',
        'approach': 'Annotation of all 736 novel\ngene families (Panaroo);\nfull annotation of all 47\nnovel genomic regions vs I6\n(Bakta + domain screening)',
        'evidence': 'Zero resistance-related\ndomains in 736 novel gene\nfamilies. All 47 novel\ngenomic regions annotated —\nno resistance cargo in any',
        'verdict': 'NOT SUPPORTED',
    },
    {
        'id': 'H8',
        'hypothesis': 'Active HGT machinery\n(conjugation / transduction)\ntransferring resistance',
        'approach': 'nucmer whole-genome\nalignment; Bakta annotation\nof novel regions; NCBI\nBLAST of ICE donor origins',
        'evidence': '2 ICEs + 1 prophage (59 kb).\nClosest BLASTn matches:\nBrevibacillus agri (99.8%),\nLacrimispora saccharolytica\n(96.3%) — non-E. coli origin.\nFull conjugation machinery;\nzero resistance cargo in all',
        'verdict': 'CHARACTERISED\n(no cargo)',
    },
]

# ---------- colours ----------
VERDICT_STYLES = {
    'STRONGLY\nSUPPORTED': {'bg': '#1B5E20', 'fg': 'white', 'label': 'STRONGLY\nSUPPORTED'},
    'FALSIFIED':            {'bg': '#B71C1C', 'fg': 'white', 'label': 'FALSIFIED'},
    'NOT SUPPORTED':        {'bg': '#E65100', 'fg': 'white', 'label': 'NOT\nSUPPORTED'},
    'CHARACTERISED\n(no cargo)': {'bg': '#0D47A1', 'fg': 'white', 'label': 'CHARACTERISED\n(no cargo)'},
}

COL_HEADER_BG  = '#1A237E'   # dark navy
ROW_ALT        = '#F5F5F5'   # alternating row tint
ROW_H4         = '#E8F5E9'   # highlight for H4
BORDER_COLOR   = '#BDBDBD'

# ---------- layout ----------
N_ROWS  = len(HYPOTHESES)
N_COLS  = 5

# Column widths (relative, sum=1)
COL_W = [0.055, 0.165, 0.215, 0.315, 0.150]   # H#, Hypothesis, Approach, Evidence, Verdict
# will be converted to absolute below

FIG_W   = 22     # inches
FIG_H   = 18     # inches
MARGIN_L = 0.04
MARGIN_R = 0.02
MARGIN_T = 0.07
MARGIN_B = 0.04

HEADER_H = 0.065   # fraction of fig height
ROW_H    = (1 - MARGIN_T - MARGIN_B - HEADER_H) / N_ROWS

fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis('off')

# helpers
def col_x(ci):
    """Left x of column ci."""
    return MARGIN_L + sum(COL_W[:ci])

def row_y(ri):
    """Bottom y of data row ri (ri=0 is top data row)."""
    return 1 - MARGIN_T - HEADER_H - (ri + 1) * ROW_H

def add_cell(ax, x, y, w, h, text, bg, fg='black', fontsize=8.5,
             bold=False, va='center', ha='center', pad=0.012):
    rect = patches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle="square,pad=0",
        facecolor=bg, edgecolor=BORDER_COLOR, linewidth=0.5,
        transform=ax.transAxes, clip_on=False)
    ax.add_patch(rect)
    ax.text(x + w / 2, y + h / 2, text,
            ha=ha, va=va, fontsize=fontsize,
            fontweight='bold' if bold else 'normal',
            color=fg, transform=ax.transAxes,
            multialignment='center',
            wrap=False,
            linespacing=1.35)

def add_verdict_badge(ax, x, y, w, h, style):
    pad_x = w * 0.12
    pad_y = h * 0.18
    bx = x + pad_x
    by = y + pad_y
    bw = w - 2 * pad_x
    bh = h - 2 * pad_y
    rect = patches.FancyBboxPatch(
        (bx, by), bw, bh,
        boxstyle="round,pad=0.005",
        facecolor=style['bg'], edgecolor='none',
        transform=ax.transAxes, clip_on=False,
        zorder=3)
    ax.add_patch(rect)
    # cell background first
    bg_rect = patches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle="square,pad=0",
        facecolor='white', edgecolor=BORDER_COLOR, linewidth=0.5,
        transform=ax.transAxes, clip_on=False, zorder=2)
    ax.add_patch(bg_rect)
    # re-draw badge on top
    ax.add_patch(rect)
    ax.text(bx + bw / 2, by + bh / 2, style['label'],
            ha='center', va='center', fontsize=8,
            fontweight='bold', color=style['fg'],
            transform=ax.transAxes,
            multialignment='center', linespacing=1.3,
            zorder=4)

# ---------- title ----------
ax.text(0.5, 1 - MARGIN_T * 0.35,
        'Supplementary Table S1   |   Hypothesis testing: mechanism of spectinomycin resistance in Paenibacillus macerans (H1–H8)',
        ha='center', va='top', fontsize=12, fontweight='bold',
        color='#1A237E', transform=ax.transAxes)

# ---------- subtitle ----------
ax.text(0.5, 1 - MARGIN_T * 0.35 - 0.022,
        'Eight mechanistic hypotheses were evaluated by whole-genome sequencing and comparative genomics of isolates from mice m5–m8 (m7 excluded; Section 2.8).',
        ha='center', va='top', fontsize=8.5, color='#555555',
        fontstyle='italic', transform=ax.transAxes)

# ---------- column headers ----------
HEADERS = ['', 'Hypothesis', 'Analytical approach', 'Key evidence', 'Verdict']
header_y = 1 - MARGIN_T - HEADER_H
for ci, (hdr, cw) in enumerate(zip(HEADERS, COL_W)):
    cx = col_x(ci)
    rect = patches.FancyBboxPatch(
        (cx, header_y), cw, HEADER_H,
        boxstyle="square,pad=0",
        facecolor=COL_HEADER_BG, edgecolor='white', linewidth=1.0,
        transform=ax.transAxes, clip_on=False)
    ax.add_patch(rect)
    ax.text(cx + cw / 2, header_y + HEADER_H / 2, hdr,
            ha='center', va='center', fontsize=9.5,
            fontweight='bold', color='white',
            transform=ax.transAxes)

# ---------- data rows ----------
for ri, hyp in enumerate(HYPOTHESES):
    y = row_y(ri)
    is_h4 = hyp['id'] == 'H4'
    row_bg = ROW_H4 if is_h4 else (ROW_ALT if ri % 2 == 0 else 'white')

    # H# badge
    cx = col_x(0)
    cw = COL_W[0]
    rect = patches.FancyBboxPatch(
        (cx, y), cw, ROW_H,
        boxstyle="square,pad=0",
        facecolor=row_bg, edgecolor=BORDER_COLOR, linewidth=0.5,
        transform=ax.transAxes, clip_on=False)
    ax.add_patch(rect)
    # coloured dot for H4
    dot_color = '#1B5E20' if is_h4 else ('#B71C1C' if hyp['verdict'] == 'FALSIFIED'
                 else '#0D47A1' if 'CHARACTERISED' in hyp['verdict']
                 else '#E65100')
    circle = patches.Circle(
        (cx + cw / 2, y + ROW_H / 2), min(cw, ROW_H) * 0.33,
        facecolor=dot_color, transform=ax.transAxes, clip_on=False, zorder=3)
    ax.add_patch(circle)
    ax.text(cx + cw / 2, y + ROW_H / 2, hyp['id'],
            ha='center', va='center', fontsize=8,
            fontweight='bold', color='white',
            transform=ax.transAxes, zorder=4)

    # text columns: Hypothesis, Approach, Evidence
    for ci_offset, key, fs in [(1, 'hypothesis', 8.5), (2, 'approach', 8.0), (3, 'evidence', 8.0)]:
        cx = col_x(ci_offset)
        cw = COL_W[ci_offset]
        add_cell(ax, cx, y, cw, ROW_H,
                 hyp[key], bg=row_bg,
                 fg='#1A237E' if is_h4 else 'black',
                 fontsize=fs,
                 bold=(ci_offset == 1 and is_h4),
                 va='center', ha='center')

    # verdict badge
    vkey = hyp['verdict']
    style = VERDICT_STYLES.get(vkey, VERDICT_STYLES['NOT SUPPORTED'])
    add_verdict_badge(ax, col_x(4), y, COL_W[4], ROW_H, style)

# ---------- horizontal separator after H4 header line ----------
# Light rule under each row
for ri in range(N_ROWS):
    y = row_y(ri) + ROW_H
    if ri == 0:  # top of first row = bottom of header — skip (already has border)
        continue
    ax.plot([MARGIN_L, 1 - MARGIN_R], [y, y],
            color=BORDER_COLOR, linewidth=0.4,
            transform=ax.transAxes)

# Bold rule under H4
h4_bottom = row_y(3)
ax.plot([MARGIN_L, 1 - MARGIN_R], [h4_bottom, h4_bottom],
        color='#1B5E20', linewidth=1.5,
        transform=ax.transAxes)
ax.plot([MARGIN_L, 1 - MARGIN_R], [h4_bottom + ROW_H, h4_bottom + ROW_H],
        color='#1B5E20', linewidth=1.5,
        transform=ax.transAxes)

# ---------- outer border ----------
border = patches.FancyBboxPatch(
    (MARGIN_L, MARGIN_B), 1 - MARGIN_L - MARGIN_R, 1 - MARGIN_T - MARGIN_B,
    boxstyle="square,pad=0",
    facecolor='none', edgecolor='#1A237E', linewidth=1.5,
    transform=ax.transAxes, clip_on=False)
ax.add_patch(border)

# ---------- legend ----------
legend_items = [
    (VERDICT_STYLES['STRONGLY\nSUPPORTED']['bg'], 'Strongly supported'),
    (VERDICT_STYLES['FALSIFIED']['bg'],            'Falsified'),
    (VERDICT_STYLES['NOT SUPPORTED']['bg'],        'Not supported'),
    (VERDICT_STYLES['CHARACTERISED\n(no cargo)']['bg'], 'Characterised — no resistance cargo'),
]
lx = MARGIN_L + 0.01
ly = MARGIN_B - 0.001
for i, (col, lbl) in enumerate(legend_items):
    rx = lx + i * 0.18
    patch = patches.FancyBboxPatch(
        (rx, ly + 0.005), 0.018, 0.018,
        boxstyle="round,pad=0.002",
        facecolor=col, edgecolor='none',
        transform=ax.transAxes, clip_on=False)
    ax.add_patch(patch)
    ax.text(rx + 0.023, ly + 0.013, lbl,
            ha='left', va='center', fontsize=7.5, color='#444444',
            transform=ax.transAxes)

# highlighted row note
ax.text(0.99, ly + 0.013,
        '★  Green-shaded row = sole supported mechanism (H4)',
        ha='right', va='center', fontsize=7.5, color='#1B5E20',
        fontstyle='italic', transform=ax.transAxes)

# ---------- save ----------
out_png = os.path.join(OUT_DIR, 'figS1_hypothesis_table.png')
out_pdf = os.path.join(OUT_DIR, 'figS1_hypothesis_table.pdf')
plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(out_pdf, bbox_inches='tight', facecolor='white')
print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")
