"""Author `notebooks/reproduce_figures.ipynb` from readable cell sources.

The notebook is generated rather than hand-edited: an executed notebook is
megabytes of JSON with base64 images in it, which is not a thing to edit by
hand or to read a diff of. This file is the editable source. Rebuild with:

    python tools/build_notebook.py            # write the notebook
    python tools/build_notebook.py --execute  # write it and run it end to end

Executing needs `jupyter` (`pip install jupyter nbconvert`) and takes a few
minutes, most of it Figure 1C-D reading 10.4 million rows.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "notebooks" / "reproduce_figures.ipynb"

C: list[dict] = []


def md(text: str) -> None:
    C.append({"cell_type": "markdown", "id": f"md{len(C):03d}", "metadata": {},
              "source": text.strip("\n").splitlines(keepends=True)})


def code(text: str) -> None:
    C.append({"cell_type": "code", "id": f"cd{len(C):03d}", "execution_count": None,
              "metadata": {}, "outputs": [],
              "source": text.strip("\n").splitlines(keepends=True)})


# ═══════════════════════════════════════════════════════════════════════════
# Front matter
# ═══════════════════════════════════════════════════════════════════════════
md("""
# Reproducing the figures, one panel at a time

This notebook rebuilds every data-derived panel of the manuscript from the
tables in `data/`, **one panel per cell**, with the reasoning next to the
picture.

## How to read it

Every panel is one cell that runs `run_panel(...)`, and every one prints the
same three steps:

| | | |
|---|---|---|
| **Step 1** | *what runs* | the script and function about to be called, the shell command that runs the same thing, the table it reads, and the pipeline script that made that table |
| **Step 2** | *the code* | the source of that function, read out of `figures/*.py` at the moment the cell runs |
| **Step 3** | *run it* | the function called for real — then the figure it returned, and a list of the files it wrote while the cell was running |

**No picture in this notebook is loaded from disk.** What you see under Step 3
is the matplotlib figure object the function just returned, rendered inline. The
file list under it is built by comparing modification times in `figures/_out/`
before and after the call, so it names the PNG and EPS that *this cell* wrote,
not files that happened to be lying there.

The notebook keeps no copy of the plotting code — it prints the real source and
then calls it — so what you read is exactly what ran, and exactly what
`make figures` runs.

Analysis and drawing are kept apart: everything that computes a number lives in
the `succession` package (`io`, `jacobian`, `diversity`, `stats`, `anchors`,
`timeaxis`) and `figures/*.py` only draws. Before each figure, the analysis
functions its panels rest on are printed too, so the chain

> raw reads → `pipeline/` (R) → `data/*.csv` → `src/succession/` → `figures/` → panel

is visible end to end rather than asserted. Where the paper claims a number, the
cell recomputes it and prints it.

**Two things to know before the figures**, because they decide numbers rather
than appearance: the two experimental groups are on different clocks, and the
sliding-window estimator has a warm-up rule that sets the sample size. Both are
worked through first.

**Runtime.** Most cells are fast. Figure 1C–D reads a 10.4-million-row table and
takes a minute or two; it is placed at the end of the Figure 1 section so you
can skip it. Figure 2B re-estimates every off-diagonal element at every
evaluation time and takes about a minute.
""")

md("""
---

## Setup

One cell. It puts `src/` and `figures/` on the path, applies the shared
matplotlib style, and defines the two helpers the rest of the notebook uses:
`show_code`, which prints a function's real source, and `run_panel`, which does
the three steps above. Both are short enough to read.
""")

code('''
%matplotlib inline

import inspect
import sys
import time
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
from IPython.display import Markdown, display

# Works whether the kernel starts in notebooks/ or at the repo root.
ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
for sub in ("src", "figures"):
    sys.path.insert(0, str(ROOT / sub))

from succession import anchors, diversity, io, jacobian, stats, style
from succession.config import COLONISED, CONTROLS, COHORT_1, COHORT_2, WINDOW
from succession.timeaxis import to_days

style.apply()

# The panels save at savefig.dpi = 200. Inline rendering uses figure.dpi, so it
# is raised here to keep what you read close to what the file looks like. It
# changes resolution only -- matplotlib lays out in inches, so nothing moves.
mpl.rcParams["figure.dpi"] = 140

OUT = ROOT / "figures" / "_out"


def show_code(fn, why=None):
    """Print a function's source, read from the file it lives in, right now."""
    path = Path(inspect.getmodule(fn).__file__).relative_to(ROOT).as_posix()
    head = f"**`{path}` &rarr; `{fn.__name__}()`**"
    display(Markdown(f"{head}<br><small>{why}</small>" if why else head))
    print(inspect.getsource(fn))


def _arg(a):
    """Short label for a call argument -- a whole DataFrame is not a label."""
    r = repr(a)
    return r if len(r) <= 30 else f"<{type(a).__name__}>"


def _written_since(t0):
    """PNG/EPS in figures/_out/ whose mtime says this call wrote them."""
    return sorted(p for p in OUT.iterdir()
                  if p.suffix in (".png", ".eps") and p.stat().st_mtime >= t0)


def run_panel(fn, *args, reads=None, made_by=None, **kwargs):
    """Name the script, print its code, run it, show what it returned.

    The figure displayed is the object `fn` returns, not a file read back --
    `style.save` closes the figure but hands it back intact, and a closed
    figure still renders. The file list is whatever appeared in figures/_out/
    while the call was running, so it is evidence rather than a claim.
    """
    path = Path(inspect.getmodule(fn).__file__).relative_to(ROOT).as_posix()
    call = f"{fn.__name__}({', '.join(_arg(a) for a in args)})"

    step1 = ["**Step 1 — what runs**",
             f"- **`{path}` &rarr; `{call}`**",
             f"- same thing from a shell: `python {path}`"]
    if reads:
        step1.append(f"- reads: {reads}")
    if made_by:
        step1.append(f"- that table was made by: {made_by}")
    display(Markdown("\\n".join(step1)))

    display(Markdown(f"**Step 2 — the code, read out of `{path}` just now**"))
    print(inspect.getsource(fn))

    display(Markdown("**Step 3 — run it**"))
    t0 = time.time()
    result = fn(*args, **kwargs)
    elapsed = time.time() - t0

    figs = (list(result.values()) if isinstance(result, dict)
            else list(result) if isinstance(result, (list, tuple))
            else [result])
    for fig in figs:
        display(fig)

    print(f"{call} returned {len(figs)} figure(s) in {elapsed:.1f} s.")
    print("files this cell wrote:")
    for p in _written_since(t0):
        print(f"    {p.relative_to(ROOT).as_posix():<46}"
              f"{p.stat().st_size / 1024:>8.0f} KB")
    return result


print("repo root:", ROOT)
print("colonised:", COLONISED, "| controls:", CONTROLS, "| window:", WINDOW)
''')

md("""
### What is where

`data/` is shipped and is all a figure needs. `pipeline/` is the R and genomics
code that produced it, kept for provenance and **not run by this notebook** —
it reads raw sequencing data that is not in the repository.
""")

code('''
for d in ("data", "src/succession", "figures", "pipeline", "tests"):
    files = sorted(p for p in (ROOT / d).rglob("*")
                   if p.is_file() and "__pycache__" not in p.parts)
    size = sum(p.stat().st_size for p in files) / 1e6
    print(f"{d:<18} {len(files):>4} files   {size:>7.1f} MB")

print("\\nshipped tables:")
for p in sorted((ROOT / "data").rglob("*.csv*")) + sorted((ROOT / "data").rglob("*.fasta")):
    print(f"    {p.relative_to(ROOT / 'data').as_posix():<44}"
          f"{p.stat().st_size / 1024:>9.0f} KB")
''')

# ═══════════════════════════════════════════════════════════════════════════
# 1. The two clocks
# ═══════════════════════════════════════════════════════════════════════════
md("""
---

## 1. The two clocks

Every table is indexed by a sampling index, and the index does **not** mean the
same thing in the two groups.

Colonised mice were sampled three times on the first day and daily afterwards,
so their index is a *sampling slot*: 1 = 3 h, 2 = 6 h, 3 = 12 h, 4 = day 1, and
slot *k* = day *k* − 3 thereafter. Control mice were sampled once a day from day
1, so their index **is** the day.

`to_days` therefore has no default for `group`, and no function in the module
guesses. The table below is why.
""")

code('''
show_code(to_days, "the conversion every time axis in the paper depends on")
''')

code('''
idx = np.arange(1, 11)
pd.DataFrame({
    "sampling index":      idx,
    "colonised → days":    to_days(idx, group="colonised"),
    "control → days":      to_days(idx, group="control"),
}).set_index("sampling index")
''')

md("""
Reading control indices with the colonised converter compresses days 1–10 into
0.125–7 d, misplacing every control point on every axis built from it. It also
silently resized the Figure 4C test, because that test's comparison window had
been defined relative to the control span — section 6 works through what that
did and did not change.

Asking for a clock that does not exist is an error rather than a fallback:
""")

code('''
try:
    to_days(4, group="colonized")     # American spelling — not a valid group
except ValueError as e:
    print("ValueError:", e)
''')

# ═══════════════════════════════════════════════════════════════════════════
# 2. The estimator
# ═══════════════════════════════════════════════════════════════════════════
md("""
---

## 2. The estimator, and the rule that sets the sample size

The community Jacobian is estimated on a dense grid (0.1 sampling slots) from
log₁₀ abundances *z*:

$$J[i \\leftarrow j](t) \\;=\\; \\mathrm{cov}\\!\\left(\\frac{dz_i}{dt},\\; z_j\\right)
\\quad\\text{over the half-open window } (t-w,\\; t]$$

*i* is the target whose rate responds, *j* the driver whose abundance acts.

This is a covariance, **not** a normalised Jacobian — the textbook
`/ var(z_j)` denominator is dropped to match the published analysis. That is
defensible for numerical stability, but not because the denominator is close to
one: during collapse the resident drivers are near-constant, `var(z_j)`
approaches zero, and the normalised slope diverges.

This is the analysis behind Figures 2, 3 and 4, so it is printed in full rather
than described. Three functions: build the state matrix, choose the evaluation
times, compute the coefficients.
""")

code('''
show_code(jacobian.build_state,
          "clone and 16S trajectories onto one dense grid, in log10 space")
''')

md("""
An evaluation time owns the window `(t − w, t]`. Near the start of a series that
window runs off the beginning of the data, making it an *expanding* window
rather than a sliding one. The warm-up rule discards those, and it is what sets
the published sample size.
""")

code('''
show_code(jacobian.evaluation_times, "the warm-up rule, and nothing else")
show_code(jacobian.offdiagonal, "every off-diagonal J[i ← j] at each evaluation time")
''')

code('''
states = {m: jacobian.build_state(m) for m in COLONISED}

rows = []
for w in (1, 2, 3, 4, 5, 10):
    n = sum(len(jacobian.evaluation_times(m, states[m], window=w))
            for m in COLONISED)
    rows.append({"window": w, "evaluations kept (n)": n})
pd.DataFrame(rows).set_index("window")
''')

md("""
Each extra step of window width costs exactly one evaluation per mouse. The
published main figures use **window 5**, which leaves **n = 113** — the sample
size quoted in the Figure 2C caption.

Turning the rule off entirely takes the pooled Figure 2C correlation from
ρ = 0.73 to 0.60 and the number of individually significant mice from 7 of 8 to
5 of 8. Sweeping the width across 1–10 shows window 5 is the maximum of that
sweep, which is worth knowing when reading the headline coefficient.
""")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 1
# ═══════════════════════════════════════════════════════════════════════════
md("""
---

# Figure 1 — colonisation of the barcoded *E. coli*

Figure 1 is the one figure with no counterpart in the manuscript's Python code,
which contains no barcode-trajectory or co-clustering script. Its panels are
drawn here from the tables the R pipeline produced.

**Panel A** is the experimental schematic and has no code.

Everything in this section comes out of `figures/figure1.py`, so a shell can
rebuild all five panels at once with `python figures/figure1.py`. The cells
below call the panel functions individually, which is the same code on the same
data.
""")

code('''
import figure1
''')

md("""
### 1B — *E. coli* load

Colony-forming units per gram of faeces, all eight colonised mice. Establishes
that the colonisation took and that the load is broadly comparable across mice,
so later per-mouse differences are not simply differences in inoculum success.

Note this is the one table indexed in **hours**, not sampling slots.
""")

code('''
run_panel(figure1.panel_b,
          reads="`data/cfu/cfu_m1-m4.csv`, `data/cfu/cfu_m5-m8.csv`",
          made_by="plate counts — measured, not derived from any script")
''')

md("""
### 1E–F — dominant clonal clusters

Barcodes are grouped into clonal clusters by correlation on their log₁₀
frequency trajectories; the line drawn is the LOESS consensus of each cluster.
This is the resolution at which lineages enter the interaction model — the
Jacobian's *E. coli* dimension is these clusters, not individual barcodes.

The loader applies the m2/m8 C1↔C2 relabel: clusters are meant to be numbered by
mean frequency at the end of the experiment, and in those two mice the stored
numbering has them the wrong way round. That relabel is a loader rule, not a
drawing rule, so it is printed here too.
""")

code('''
show_code(io.load_clone_loess, "and the relabel it applies")
''')

code('''
run_panel(figure1.panels_ef,
          reads="`data/barcodes/clusters/<mouse>_loess_clusters.csv`",
          made_by="`pipeline/03_clustering/01_hclust.R` → `02_select.R`")
''')

md("""
### 1G–H — 16S community composition

Family-level relative abundance, cohort 1 then cohort 2. This is where the
central observation of the paper is visible: the resident community collapses
under antibiotic, and *Paenibacillaceae* rises to dominate before the community
re-establishes.

Two display rules, both reconstructed from the published panel: zeros get a 1e-6
pseudocount so a vanishing family drops to the axis floor instead of leaving a
gap, and a *run* of zeros becomes NaN so the line breaks rather than being
bridged by a straight segment. Only families whose own mean in that mouse
exceeds 1e-3 are drawn — applied per mouse, because Akkermansiaceae averages
2e-5 in m1 but 7.9e-3 in m5.

`panels_gh` returns one figure per cohort, so Step 3 shows two.
""")

code('''
run_panel(figure1.panels_gh,
          reads="`data/16s/family/<mouse>_family.csv`",
          made_by="`pipeline/01_16s/02_process.R` → `03_filter.R`")
''')

md("""
The twelve families in the legend are the published set. Worth being explicit:
**no abundance threshold reproduces exactly that set.** Lachnospiraceae peaks at
22.2% and is excluded, while Marinifilaceae peaks at 2.0% and is included.
Applying the per-mouse gate to all 49 families detected in the colonised mice
yields these twelve plus Lachnospiraceae and Rikenellaceae — so the published
legend is curated, not thresholded. It is reproduced here as published.
""")

md("""
### 1I — co-clustering of lineages with families

UPGMA on shape-based distances between every *E. coli* clonal trajectory and
every bacterial family trajectory within a mouse. Clones and families that land
in the same clade rise and fall together, which is the qualitative motivation
for estimating a joint interaction model over both.
""")

code('''
run_panel(figure1.panel_i,
          reads="`data/coclustering/<mouse>_sbd_distance.csv`",
          made_by="`pipeline/04_coclustering/01_coclustering.R` → `02_all_sbd.R`")
''')

md("""
### 1C–D — every barcode

**Slow cell** — reads 10.4 million rows, roughly a minute or two.

Every barcode in every mouse. The caption colours the 1000 most abundant; the
stored colour assignment is broader (~4000 per mouse), which floods the panel
and hides the grey background of the published figure, so barcodes are ranked by
peak frequency and the top 1000 kept. Each keeps its stored colour, so a lineage
looks the same in every mouse it appears in.
""")

code('''
run_panel(figure1.panels_cd,
          reads="`data/barcodes/trajectories.csv.gz` (10.4 M rows, 36 MB)",
          made_by="`pipeline/02_barcode/02_process.R`, twelve per-sample tables "
                  "concatenated to five columns")
''')

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 2
# ═══════════════════════════════════════════════════════════════════════════
md("""
---

# Figure 2 — the shift tracks community diversity

From here on, every panel is a port of the manuscript's own Python code.

Two pieces of analysis sit under this figure: Hill ¹D, and the phase anchor that
marks when *Paenibacillaceae* takes over. Both are printed before the panels
that use them.
""")

code('''
import figure2

show_code(diversity.hill_q1, "Hill ¹D = exp(Shannon)")
show_code(diversity.hill_q1_from_taxa,
          "computed from the curated taxon set that feeds the Jacobian — "
          "not the pipeline's own diversity table; Figure 2C needs this one")
show_code(anchors.paeni_onset, "the dashed line in 2A: first slot at 10% of the mouse's own maximum")
''')

md("""
### 2A — diversity and the *Paenibacillaceae* transition

Per mouse: 16S Hill ¹D, barcode lineage Hill ¹D (log₁₀), and the relative
abundance of *Paenibacillaceae* and Enterobacteriaceae on the right axis. The
dashed vertical line is the *Paenibacillaceae* onset.

The 3 h barcode point is dropped for m1, m2, m4 and m5. Those samples exist but
are built from a handful of reads: Hill ¹D of 5.7–13.0 against 10⁴–10⁵ six hours
later. That is a sampling artefact, not low diversity, and plotting it would put
a point near 1 on a log₁₀ axis whose real signal starts near 4. They are *kept*
in the barcode LOESS feeding the Jacobian, because dropping them would start
cohort 1's dense grid a slot later and change every downstream evaluation count.
""")

code('''
run_panel(figure2.panel_a,
          reads="`data/16s/diversity.csv`, `data/barcodes/diversity.csv`, "
                "`data/16s/family/`",
          made_by="`pipeline/01_16s/05_diversity.R`, "
                  "`pipeline/02_barcode/04_diversity.R`")
''')

code('''
anchors.table().set_index("mouse").round(3)
''')

md("""
### 2C — inhibition against diversity

The central quantitative claim. For each evaluation time, the mean of the
**negative** off-diagonal Jacobian elements — inhibitory interaction strength —
against 16S Hill ¹D at the nearest sampling slot. Less diverse communities
interact more inhibitorily.

Two details decide whether the published coefficient reproduces:

- Hill ¹D comes from the **curated taxon set** that feeds the Jacobian
  (`16s/for_jacobian/`), not from the pipeline's own diversity table. The two
  differ by a mean of 0.010 at every timepoint, and the wrong source misses the
  per-mouse values.
- `mean_negative` averages only the negative elements. It is not `mean`, which
  averages inhibitory and facilitative together and sits near zero.

The manuscript reports ρ = 0.73, n = 113, significant in 7 of 8 mice.
""")

code('''
show_code(jacobian.summarise, "where mean_negative comes from")
show_code(figure2.figure2c_points, "the points the panel plots and the test uses")
''')

code('''
run_panel(figure2.panel_c,
          reads="`data/jacobian/<mouse>_jacobian_timeseries.csv`, "
                "`data/16s/for_jacobian/<mouse>_16S_taxa.csv`",
          made_by="`pipeline/05_jacobian/01_jacobian.R`, `pipeline/01_16s/03_filter.R`")
''')

md("""
And the published statistics, recomputed from those same points:
""")

code('''
pts = figure2.figure2c_points()
rho, p, n = stats.spearman(pts["q1"], pts["mean_negative"])
print(f"pooled: rho = {rho:.3f}   P = {p:.2e}   n = {n}")

per = pd.DataFrame([
    {"mouse": m, "rho": r, "P": pv, "n": nn}
    for m, g in pts.groupby("mouse")
    for r, pv, nn in [stats.spearman(g["q1"], g["mean_negative"])]
]).set_index("mouse")
per["significant"] = per["P"] < 0.05
print(f"significant in {int(per['significant'].sum())} of {len(per)} mice")
per.round(4)
''')

md("""
The asymptotic p-values above are anticonservative — these are autocorrelated
timecourses, and the manuscript quotes them anyway. `stats.circular_shift_p`
gives a null that preserves the serial structure by circularly shifting one
series, and should be preferred when a claim rests on significance:
""")

code('''
show_code(stats.circular_shift_p)
''')

code('''
pd.DataFrame([
    {"mouse": m,
     "asymptotic P": stats.spearman(g["q1"], g["mean_negative"])[1],
     "circular-shift P": stats.circular_shift_p(g["q1"], g["mean_negative"])}
    for m, g in pts.groupby("mouse")
]).set_index("mouse").round(4)
''')

md("""
### 2B — the distribution behind the mean

**Slow cell** — re-estimates every off-diagonal element at every evaluation time.

Panel C collapses each evaluation to one number. This panel shows the whole
distribution of $J_{ij}$ as a ridgeline, one row per evaluation, coloured by
whether it falls before or after the *Paenibacillaceae* onset. It is the
evidence that the shift is a change in the shape of the interaction
distribution, not only in its mean.

It uses a **different estimator** from panel C — windowed by slot count, and
expanding rather than warm-up-discarded — which is a deliberate choice, not an
inconsistency. The function says why, so it is printed first.
""")

code('''
show_code(jacobian.offdiagonal_by_slot,
          "the Figure 2B estimator, and why it is not the panel-C one")
''')

md("""
**Every evaluation is drawn here, not only the ones panel C uses.** The warm-up
rule exists so that a *correlation* is computed over a consistent window width;
it is not a statement that the early windows are wrong. Withholding them from a
panel whose whole purpose is to show the collapse would hide the phase of
interest, so B shows 17/16/15 rows per mouse where C uses 15/14/13. Pass
`warmup_tolerance=WARMUP_TOLERANCE` to reproduce the panel-C subset.

**The scaling is per mouse, not per ridge**, and this is the choice that decides
what the panel says. One scale — `RIDGE_SCALE / median(peak densities)` — is
shared by every ridge in a panel, so a concentrated distribution draws tall and
a diffuse one draws flat. Normalising each ridge to its own maximum instead
would force every row to the same height and delete exactly the signal the
panel exists to show. (Both variants exist in the manuscript's code; this
follows `12_ridgeline_jacobian.py`, the script behind the submitted panel.)

Ridges sit at their **real day** on a continuous axis rather than on evenly
spaced rows, so the sub-day points crowd near the origin and the daily samples
spread out. That also aligns the eight panels with no padding — a mouse with no
evaluation at a time simply has no ridge there:

| | first evaluation | last evaluation |
|---|---|---|
| m1–m4 | 6 h | 15.9 d |
| m5 | 6 h | 14.9 d |
| m6–m8 | **12 h** | 14.9 d |

The x-axis is clipped to the submitted panel's view, (−0.42, 0.22). The full
range of the data reaches −1.68, but 99% of the mass lies within ±0.5, so
drawing it uncut compresses every ridge into a spike; `panel_b(xlim=None)`
shows the full range.
""")

code('''
run_panel(figure2.panel_b,
          reads="`data/jacobian/`, `data/barcodes/clusters/`, `data/16s/for_jacobian/`",
          made_by="`pipeline/05_jacobian/01_jacobian.R`")
''')

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 3
# ═══════════════════════════════════════════════════════════════════════════
md("""
---

# Figure 3 — stability, and who acts on whom

Panels A and B are eigenvalues of the stored per-window Jacobians; C and D are
single directed coefficients. Both pieces of analysis first.
""")

code('''
import figure3

show_code(jacobian.eigenvalues, "A and B")
show_code(jacobian.directed_group, "C and D — mean over a mouse's clones, then across mice")
''')

code('''
eig = figure3.all_eigenvalues()
print(f"{len(eig)} eigenvalues from {eig['mouse'].nunique()} mice")
eig.head()
''')

md("""
### 3A — the eigenvalue spectrum of one mouse

Eigenvalues of the per-window Jacobian for m1 in the complex plane, coloured by
time. The community is locally stable where every Re(λ) < 0, so the stable
region is the shaded left half-plane.

**A correction is applied here.** The submitted panel labels x as Im and y as
Re. The code that produced it computes the opposite — x = Re, y = Im — which is
also what makes the stability boundary the vertical line it is drawn as. The
panel is drawn here as the code computes; only the labels change.
""")

code('''
run_panel(figure3.panel_a, eig,
          reads="`data/jacobian/<mouse>_jacobian_matrices_by_time.csv` (m1)",
          made_by="`pipeline/05_jacobian/01_jacobian.R` → `03_stability.R`")
''')

md("""
### 3B — Re(λ) against time, every eigenvalue of every mouse

The stability trend across the whole study. The manuscript reports ρ = −0.43.
The view is clipped to the bulk — a handful of eigenvalues reach −9 — but the
correlation is computed on all of them.
""")

code('''
run_panel(figure3.panel_b, eig,
          reads="`data/jacobian/<mouse>_jacobian_matrices_by_time.csv` (all mice)",
          made_by="`pipeline/05_jacobian/01_jacobian.R` → `03_stability.R`")
''')

code('''
rho, p, n = stats.spearman(eig["day"], eig["re"])
print(f"Re(lambda) vs time:  rho = {rho:.3f}   P = {p:.2e}   n = {n}")
''')

md("""
### 3C–D — directed interactions between colonizer and resident

The direction convention matters and is easy to invert:
$J[i \\leftarrow j] = \\mathrm{cov}(dz_i/dt,\\, z_j)$ — **i responds, j acts**.

- **C** is resident → colonizer: $J[\\text{clone} \\leftarrow \\textit{Paeni}]$
- **D** is colonizer → resident: $J[\\textit{Paeni} \\leftarrow \\text{clone}]$

Within a mouse the value is the mean over its clonal clusters; across mice it is
mean ± s.d. Enterobacteriaceae is the single 16S taxon standing for *E. coli* in
the community model, drawn dashed as a consistency check on the clone-based
estimate.
""")

code('''
run_panel(figure3.panels_cd,
          reads="`data/jacobian/`, `data/16s/family/`",
          made_by="`pipeline/05_jacobian/01_jacobian.R` → `04_interactions.R`")
''')

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 4
# ═══════════════════════════════════════════════════════════════════════════
md("""
---

# Figure 4 — the antibiotic-only controls

The controls received spectinomycin but no *E. coli*. They separate what the
antibiotic does from what the colonizer does.

**This is the figure where the clock matters.** Controls carry no barcode
dimension, so their Jacobian is estimated over 16S taxa only — the two groups
are therefore not built from the same number of variables, which is stated on
the panel itself.
""")

code('''
import figure4
''')

md("""
### 4A — control community composition

Family composition of the four control mice, on the control clock (index = day).
*Paenibacillaceae* is drawn heavier. It does **not** take over in the absence of
the colonizer, which is what makes the bloom in the colonised mice interesting
rather than a plain antibiotic effect.
""")

code('''
run_panel(figure4.panel_a,
          reads="`data/16s/family/c_m<n>_family.csv`",
          made_by="`pipeline/01_16s/02_process.R` → `03_filter.R`")
''')

md("""
### 4B — control diversity against the colonised band

16S Hill ¹D for the controls over the colonised mean ± s.d.

The gaps are real and are drawn as breaks: c_m1 has usable ¹D on days 3, 6, 7, 8
and 9 only. Reindexing onto the full day range inserts NaN so the line breaks
rather than being drawn straight across the missing days.

The `clock` argument is explicit and has no default — that is the point of
section 1.
""")

code('''
run_panel(figure4.panel_b, "corrected",
          reads="`data/16s/diversity.csv`",
          made_by="`pipeline/01_16s/05_diversity.R`")
''')

md("""
### 4C — the comparison, and the clock

Mean pairwise $J_{ij}$ over time for controls and colonised, with a one-sided
Mann–Whitney testing whether control interactions are *less negative* than
colonised over the days the controls span.

The published version put control indices through the colonised converter. Both
are computed below from the same data — only the clock differs.
""")

code('''
show_code(figure4.comparison, "the test, with the window as a named argument")
show_code(stats.mann_whitney)
''')

code('''
rows = []
for clock in ("corrected", "published"):
    c = figure4.comparison(clock)
    rows.append({
        "clock": clock,
        "controls span (days)": f"{c['min_ctl']:g} – {c['max_ctl']:g}",
        "n control": len(c["a"]),
        "n colonised": len(c["b"]),
        "one-sided P": c["p"],
        "significant at 0.05": c["p"] < 0.05,
    })
pd.DataFrame(rows).set_index("clock")
''')

md("""
Same measurements, same test, different reading of the index — and the span the
controls appear to cover changes. Under the colonised converter the controls
appear to span 0.125–7 d rather than 1–10 d.

**This repository defaults to the corrected clock.** The submitted version is
still one command away: `python figures/figure4.py --clock published`, or
`make figures-published`.
""")

code('''
run_panel(figure4.panel_c, "corrected",
          reads="`data/jacobian/c_m<n>_jacobian_timeseries.csv` and the colonised set",
          made_by="`pipeline/05_jacobian/01_jacobian.R`")
''')

md("""
The two clocks give the **same** P. That is the point worth taking away: the
clock never changed a measured value, it changed which colonised evaluations
fell inside a window that had been defined relative to the control span. Fixing
the window at 1–7 d makes the test the same test under either clock.

What the result actually depends on is the window:
""")

code('''
pd.DataFrame([
    {"colonised window (days)": f"{w[0]:g}–{w[1]:g}",
     "n colonised": len(c["b"]), "n control": len(c["a"]),
     "one-sided P": round(c["p"], 4)}
    for w in [(1, 7), (1, 9.9), (1, 16)]
    for c in [figure4.comparison("corrected", colonised_window=w)]
]).set_index("colonised window (days)")
''')

md("""
So the control comparison is significant over the collapse-and-onset period the
controls actually cover, and not significant once the colonised recovery phase
is folded in. The honest statement of the result has to name its window; an
open-ended one let a labelling error look like a finding.
""")

# ═══════════════════════════════════════════════════════════════════════════
# FIGURE 5
# ═══════════════════════════════════════════════════════════════════════════
md("""
---

# Figure 5 — the resistance genotype

Why *Paenibacillaceae* survives spectinomycin at all.

**Panels B** (the deletion mapped onto the 70S ribosome) and **C** (colony
morphology) are a structural render and a photograph. Neither has code.

Both shipped genomics tables were literals inside the manuscript's own figure
scripts rather than files; they were extracted unchanged so that
`figures/figure5.py` only draws. `pipeline/genomics/` holds the assembly,
alignment and annotation code that established them upstream — see
`pipeline/README.md`.
""")

code('''
import figure5
''')

md("""
### 5A — the RpsE deletion

Ribosomal protein S5 of the gut isolates is three residues shorter (162 aa) than
the *Paenibacillus macerans* I6 reference (165 aa), with no other difference.
S5 contacts 16S h34 at the spectinomycin binding site, so a deletion there is a
target-alteration mechanism — the resistance is chromosomal and intrinsic, not
acquired.

Which three residues are "deleted" is ambiguous. The reference reads
`…I-N-R-V-A-K-V-V-K-G-G…`, and removing either `VAK` at 20–22 or `AKV` at 21–23
gives the identical isolate sequence, because the region is a short repeat. The
caption's "ΔVAK at 20–22" and the genomics record's "Δ(A21-K22-V23)" are the
same event described from opposite ends of that ambiguity, not a contradiction.
""")

code('''
seqs = io.load_rpse_alignment()
for name, s in seqs.items():
    print(f"{len(s):>4} aa   {name}")
    print(f"          {s[:40]}...")
''')

code('''
run_panel(figure5.panel_a,
          reads="`data/genomics/rpsE_alignment.fasta`",
          made_by="`pipeline/genomics/build_rpsE_msa_input.py`, "
                  "`interpret_mafft_rpsE.py` (MAFFT, Nanopore assemblies)")
''')

md("""
### 5D — the resistance gene survey

Antimicrobial resistance gene survey of the *P. macerans* I6 genome against
CARD (RGI, Strict threshold), following the manuscript's
`32_figure_card_rgi.py`. **A** counts hits per drug class; **B** lists them.

Detected genes fall within the glycopeptide, lincosamide, fosfomycin and
disinfecting-agent classes. The point of the panel is the row with no bar:
**no spectinomycin resistance determinant was identified** — no *aadA* / ANT(9)
anywhere in the genome. That is what confirms resistance in the gut isolates is
attributable solely to the chromosomal ΔVAK deletion in *rpsE*.
""")

code('''
hits = io.load_card_hits()
print(f"{len(hits)} strict hits;", hits["drug_class"].nunique(), "drug classes")
print("aminoglycoside hits:",
      int(hits["drug_class"].str.contains("minoglycoside", case=False).sum()))
hits
''')

code('''
run_panel(figure5.panel_d,
          reads="`data/genomics/card_rgi_hits.csv`",
          made_by="CARD/RGI on the I6 assembly; the 13 strict hits as extracted "
                  "from the manuscript's `32_figure_card_rgi.py`")
''')

# ═══════════════════════════════════════════════════════════════════════════
# Closing
# ═══════════════════════════════════════════════════════════════════════════
md("""
---

## Every panel, as PNG and EPS

`style.save` writes each panel twice: a **PNG** to read and an **EPS** to edit
or send to a journal. EPS has no alpha channel, so semi-transparent fills — the
s.d. bands, the ridgeline fills — flatten against white; the PNG is the faithful
one, the EPS is the editable one. The densest panels rasterise their barcode
layer inside the EPS so the file stays openable.

The cell below checks that both formats exist for every panel in
`figures/_out/`, and says so if one is missing.
""")

code('''
show_code(style.save)
''')

code('''
names = sorted({p.stem for p in OUT.glob("*.png")} | {p.stem for p in OUT.glob("*.eps")})
rows = []
for name in names:
    png, eps = OUT / f"{name}.png", OUT / f"{name}.eps"
    rows.append({"panel": name,
                 "PNG (KB)": round(png.stat().st_size / 1024) if png.exists() else None,
                 "EPS (KB)": round(eps.stat().st_size / 1024) if eps.exists() else None})
table = pd.DataFrame(rows).set_index("panel")

missing = table[table.isna().any(axis=1)]
print(f"{len(table)} panels in figures/_out/")
print("missing a format:", "none" if missing.empty else list(missing.index))
table
''')

md("""
The supplementary panels S1–S6 are drawn by `figures/supplementary/figS*.py`
into the same directory and are not walked through here; `make supplementary`
(or `python run_all.py supplementary`) rebuilds them.

## Checking the whole thing

Everything above is also asserted as a test. The suite fails if a change to the
estimator moves a published number:

```bash
make verify              # or: python run_all.py test
```

`tests/test_reproduces_published.py` pins Figure 2C (ρ = 0.73, n = 113, 7 of 8),
Figure 3B (ρ = −0.43, n = 1545), the window sweep, and the fact that the two
clocks are distinct. `tests/test_timeaxis.py` covers the index conversions on
their own.
""")


# ═══════════════════════════════════════════════════════════════════════════
def build() -> Path:
    nb = {
        "cells": C,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python",
                           "name": "python3"},
            "language_info": {"name": "python", "pygments_lexer": "ipython3"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    NOTEBOOK.parent.mkdir(exist_ok=True)
    NOTEBOOK.write_text(json.dumps(nb, indent=1, ensure_ascii=False), encoding="utf-8")
    print(f"wrote {NOTEBOOK.relative_to(ROOT)}  ({len(C)} cells)")
    return NOTEBOOK


def execute(path: Path) -> None:
    """Run the notebook in place, so it ships with every output rendered."""
    cmd = [sys.executable, "-m", "jupyter", "nbconvert", "--to", "notebook",
           "--execute", "--inplace", "--ExecutePreprocessor.timeout=1800",
           str(path)]
    print("$", " ".join(cmd[2:]))
    subprocess.run(cmd, cwd=ROOT, check=True)
    print(f"executed {path.relative_to(ROOT)}  "
          f"({path.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--execute", action="store_true",
                    help="run the notebook after writing it")
    args = ap.parse_args()
    nb_path = build()
    if args.execute:
        execute(nb_path)
