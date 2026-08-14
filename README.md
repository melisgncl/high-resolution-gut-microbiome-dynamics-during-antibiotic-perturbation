# succession-paper

Reproduction code and derived data for

> ***Inhibitory-to-Facilitative Interaction Shift During Ecological Succession
> in Gut Microbiome Driven by Antimicrobial Resistance***
> Gencel, Matta, Hui, Marrero Cofino, Ramanathan, Menendez, Bershtein &
> Serohijos.

---

## What this repository is

Antibiotic treatment empties the gut of most of its resident bacteria. What
happens next — which organisms return, in what order, and how they act on one
another while the community reassembles — is a succession process, but one that
is difficult to watch directly, because the coloniser and the recovering
residents are usually measured on different instruments and cannot be placed in
a single community model.

The study behind this repository makes that possible by colonising
antibiotic-treated mice with a **barcoded *Escherichia coli* library** under
continuous antibiotic selection, so that the invader's population is resolved
lineage by lineage, while **16S rRNA amplicon sequencing of the same faecal
pellets** follows the resident community around it. Placing both layers on one
time axis allows a **time-varying community Jacobian** to be estimated — a
matrix of who responds to whom — and its behaviour over the course of the
succession to be measured rather than assumed.

This repository holds the derived data and the code that turn it into the
manuscript's figures and numbers. It exists so that a reader can do three things
without contacting the authors:

1. **See** the main figures without installing anything — `figures/_out/`
   ships every panel of Figures 1–5 as PNG and EPS.
2. **Rebuild** all 24 panels, main and supplementary, from the shipped tables
   with a single command, and confirm that the main figures come back
   byte-identical to what is committed.
3. **Interrogate** the result — check a published statistic against an
   assertion, change an analysis parameter and see which conclusions move, and
   read, in one place, every point at which the rebuilt figures differ from the
   ones submitted.

**No raw sequencing data is required.** Raw reads and assemblies are deposited
separately; see the manuscript's data availability statement. What ships here is
the derived layer: 41 MB of tables, documented column by column in
[`data/README.md`](data/README.md).

`24 panels` · `41 MB derived data` · `35 tests pinning the published numbers` ·
`Python 3.13`

---

## The experiment

Mice were pretreated for four weeks with an antibiotic cocktail (metronidazole,
neomycin, ampicillin, vancomycin) to deplete the resident microbiota. After
three days' recovery, eight animals were gavaged with ~10⁸ CFU of a barcoded
*E. coli* MG1655 library carrying ~10⁶ unique heritable chromosomal barcodes,
under continuous spectinomycin selection. Four control animals received the same
antibiotics and spectinomycin but no *E. coli*, separating the effect of the
drug regimen from the effect of colonisation.

| group | mice | sampling |
|---|---|---|
| colonised | `m1`–`m4` (cohort 1), `m5`–`m8` (cohort 2), run at two facilities | 3 h, 6 h, 12 h, then daily to day 16 / day 15 |
| control | `c_m1`–`c_m4` | daily, days 1–10 |

Four measurement layers were collected, and each has a directory under `data/`:

| layer | what it measures | `data/` |
|---|---|---|
| CFU | *E. coli* load per gram of faeces | `cfu/` |
| barcodes | ~10.4 M lineage-frequency observations | `barcodes/` |
| 16S | resident community at family level, plus Hill diversity | `16s/` |
| genomics | Nanopore assemblies of the recovered *Paenibacillus macerans* | `genomics/` |

---

## The analysis

Six steps. Each writes tables into `data/`, which the figure scripts then read.
The full script order is in [`pipeline/README.md`](pipeline/README.md).

**1 · 16S community composition** → `16s/family/`, `16s/diversity.csv`
Rarefaction thresholds set per cohort by sequencing depth; relative abundance at
family level; Hill numbers q0/q1/q∞ per sample.
*Figures 1G–H, 2A, 3C–D, 4A, 4B.*

**2 · Barcode lineage dynamics** → `barcodes/trajectories.csv.gz`, `barcodes/diversity.csv`
Every barcode at every timepoint, plus Hill diversity over lineages. The shipped
36 MB cache replaces ~1.1 GB of per-barcode processing state.
*Figures 1C–D, 2A, S1.*

**3 · Clonal clustering** → `barcodes/clusters/`
Barcodes with similar trajectories are grouped into **clonal clusters**
(C1, C2, …), the unit every later analysis uses. Distance is Pearson correlation
on log10 frequency with average linkage; the cutoff is chosen where the
scaled-distance and cluster-count curves cross. A LOESS consensus is fitted per
cluster on a dense grid. *Figures 1E–F, S2.*

**4 · Co-clustering of clones with families** → `coclustering/`
Shape-Based Distance on z-normalised trajectories, average linkage, with clonal
clusters and bacterial families placed in a single tree. This asks which
*E. coli* lineages track which residents. *Figures 1I, S3, S6.*

**5 · Time-varying community Jacobian** → `jacobian/`
Clones and families are treated as one community, and interaction coefficients
are estimated in a sliding window over log10 abundances on a dense grid:

```
J[i ← j](t) = cov( dz_i/dt , z_j )       over (t − w, t],  w = 5 slots
                                          i responds, j acts
```

The eigenvalues of each J give the local stability of the community at that time
(stable when Re λ < 0). Note that this is a **covariance**, not a normalised
Jacobian: the textbook `/ var(z_j)` denominator is dropped to match the
published analysis. That is defensible, but not on the usual grounds that the
denominator is near unity — during community collapse `var(z_j)` approaches
zero. *Figures 2B–C, 3A–D, 4C, S4, S5.*

**6 · Genomics of the resistant isolate** → `genomics/`
Nanopore assembly and QC, average nucleotide identity against the *P. macerans*
I6 reference, *rpsE* alignment with read-level confirmation, a CARD/RGI
resistance survey, pan-genome analysis and a survey of novel genomic regions.
*Figure 5.*

---

## Principal findings

- ***E. coli* colonises the emptied gut, and its lineages are rapidly selected
  down** — a handful of clonal clusters come to carry the population.
- **The resident community recovers, and Paenibacillaceae leads that recovery** —
  peaking at 53–67% of the 16S community in every colonised mouse, against
  ≤1.4% in the antibiotic-only controls. The bloom is therefore a feature of
  colonisation, not of the drug regimen alone.
- **Interaction strength tracks community diversity.** Mean negative *J* against
  16S Hill ¹D gives a pooled Spearman **ρ = 0.732, n = 113**, and is
  individually significant in **7 of 8** colonised mice.
- **Stability increases over the succession.** Re(λ) against time gives
  **ρ = −0.430** over 1545 eigenvalues.
- **The returning *Paenibacillus* is spectinomycin-resistant, and the mechanism
  is chromosomal rather than acquired**: a three-residue deletion in *rpsE*, the
  gene for 30S ribosomal protein S5, at the residues that contact the drug's
  binding site. The deleted lysine is conserved in *E. coli*, *Bacillus
  subtilis* and the *P. macerans* type strain, and is absent from every
  sequenced isolate. The CARD/RGI survey returns 13 strict hits and **no
  aminoglycoside-modifying determinant — no *aadA* or ANT(9) anywhere in the
  genome**.

---

## Reproducing the figures

```bash
git clone <repo-url> && cd succession-paper
conda env create -f environment.yml && conda activate succession-paper
make all
```

`pip install -r requirements.txt` is an equivalent route. If `make` is not
available — typically on Windows — every target has a twin:
`python run_all.py all`.

| command | what it does |
|---|---|
| `make all` | all 24 panels → `figures/_out/` |
| `make figures` / `make supplementary` | main figures 1–5 / supplementary S1–S6 only |
| `make notebook` | opens the panel-by-panel walkthrough |
| `make verify` | re-derives the published statistics and asserts them |
| `make figures-published` | regenerates the *submitted* Fig 4B/4C for side-by-side comparison |

**A rebuild should change nothing.** The 17 panels of Figures 1–5 are
committed, so `make figures` overwrites files already in the tree and
`git status` comes back clean. That is the reproducibility check being made
here — not "the code ran" but "the code produced the same bytes". The
supplementary panels are drawn by `make supplementary` into the same directory
and are deliberately left untracked, so they do not enter that comparison.

For a reader who would rather watch the analysis than run it,
`notebooks/reproduce_figures.ipynb` ships executed and is readable directly on
GitHub. Each panel cell prints the function that is about to run, then that
function's source read live from `figures/*.py`, then the call itself, and
renders the figure object returned. **No image in the notebook is loaded from
disk**, so it cannot drift from what the figure scripts produce.

---

## Verification against the published numbers

```bash
make verify
```

| statistic | manuscript | this repository |
|---|---|---|
| Fig 2C pooled Spearman ρ | 0.73 | **0.732** |
| Fig 2C sample size | 113 | **113** |
| Fig 2C individually significant mice | 7 of 8 | **7 of 8** (m8 the exception) |
| Fig 3B Re(λ) against time, ρ | −0.43 | **−0.430**, P = 1.5 × 10⁻⁷⁰, n = 1545 |
| Fig S4 window sweep, n range | 73–129 | **73–129** |

The figure code is a *port of the code that produced the submitted figures*, not
an independent reimplementation whose agreement would itself have to be argued.
**Figure 1 is the exception**: it has no counterpart in the submitted code, so
its panels are drawn from the pipeline's tables.

---

## Known differences from the submitted figures

These are stated here rather than left to be discovered.

**Fig 4B/4C — the control clock.** The submitted panels read control indices
through the colonised converter. This repository defaults to the corrected
clock, and `make figures-published` regenerates the submitted versions alongside
it.

The Fig 4C *test* turns out not to depend on the clock at all. It depends on the
comparison window, which the submitted analysis defined relative to the control
span — so a mislabelled axis quietly resized the test rather than changing any
measured value. With the window stated explicitly (`figure4.COLONISED_WINDOW`,
default 1–7 d), both clocks agree exactly:

| colonised window | corrected clock | published clock |
|---|---|---|
| **1–7 d** (default) | P = 0.019 | P = 0.019 |
| 1–9.9 d | P = 0.118 | P = 0.118 |
| full series | P = 0.266 | P = 0.266 |

The controls do differ from the colonised animals over the collapse-and-onset
period they actually cover, and do not once the colonised recovery phase is
folded in. Either statement is defensible; only a stated window makes the choice
checkable.

**Fig 3A — axis labels.** The submitted panel labels x = Im, y = Re; the code
computes the opposite. It is drawn here as computed, with the stable half-plane
(Re < 0) shaded.

**S6 is an addition, not a reproduction** — the submitted set runs S1–S5.

---

## Reading the data correctly

**The two groups are on different clocks.** The colonised index is a *sampling
slot* (1 = 3 h, 2 = 6 h, 3 = 12 h, 4 = day 1, then slot *k* = day *k* − 3); the
control index is *the day itself*. `timeaxis.to_days` therefore requires `group`
and deliberately has no default. Passing control indices through the colonised
converter compresses days 1–10 into 0.125–7 d and misplaces every control point
on every time axis.

**Every choice that moves a result is a named argument** rather than a buried
literal (`src/succession/jacobian.py`):

| parameter | default | effect |
|---|---|---|
| `window` | 5 | window width in sampling slots |
| `pseudocount` | 1e-4 | added before the log; the 16S tables are mostly zeros, so this is not cosmetic |
| `warmup_tolerance` | 2.0 | keeps an evaluation only if `t − window ≥ grid_start − tol` |
| `min_points` | 4 | minimum dense points per window |

`warmup_tolerance` is what sets n — 129 at widths 1–3, **113 at 5**, 73 at 10.
Disabling it takes the pooled Fig 2C ρ from 0.73 to 0.60, and the number of
individually significant mice from 7 of 8 to 5 of 8. **Figure S4 sweeps the
window width so that this choice is visible rather than asserted.**

**Three traps are documented in `data/README.md`**, and each has produced a
wrong number at some point: two Hill q1 series exist for the same mice and are
not interchangeable (Fig 2C needs the recomputed one); clusters C1 and C2 are
stored the wrong way round in `m2` and `m8`; and the 3 h barcode samples in four
mice are near-empty rather than genuinely low-diversity.

---

## Repository layout

```
data/            derived tables, ~41 MB — 16s/ barcodes/ cfu/ coclustering/ jacobian/ genomics/
src/succession/  the package: loading, analysis, the estimator. No plotting except style.py
figures/         one script per figure. Drawing only, no analysis
figures/_out/    Figures 1–5, every panel as PNG + EPS — committed
                 (S1–S6 regenerate here, untracked)
pipeline/        provenance: the code that produced data/. Cannot run from a clone
notebooks/       the walkthrough, shipped executed
tools/           builds that notebook — it is generated, never hand-edited
tests/           the published numbers, as assertions
assets/          the figures as submitted, for comparison
```

Two conventions keep the boundary sharp: **no analysis in `figures/`** and **no
plotting in `src/`** (except `style.py`). The consequence is that "where does
this number come from?" is always answerable without opening a figure script.

| I want to… | Go to |
|---|---|
| see the figures without installing anything | `figures/_out/` — Figures 1–5 as PNG + EPS, committed |
| watch the analysis run, panel by panel | `notebooks/reproduce_figures.ipynb` |
| know what a table's columns mean | **`data/README.md`** — the data dictionary |
| find how a number was computed | `src/succession/` — never in a figure script |
| find how a panel was drawn | `figures/figure{1..5}.py`, `figures/supplementary/` |
| check a published statistic | `tests/test_reproduces_published.py` |
| see how the derived tables were made | `pipeline/` + `pipeline/README.md` |
| compare against the submitted version | `assets/figure_N_as_submitted.png` |
| change a parameter and see what moves | `src/succession/config.py` — every constant documented with its effect |

### Which figure reads which data

| Fig | Panels | Reads | Notes |
|---|---|---|---|
| 1 | A | — | schematic, no code |
| | B | `cfu/` | |
| | C–D | `barcodes/trajectories.csv.gz` | top 1000 by peak frequency coloured, rest grey |
| | E–F | `barcodes/clusters/` | m2 and m8 carry a C1↔C2 relabel |
| | G–H | `16s/family/` | per-mouse mean > 1e-3, 12 legend families |
| | I | `coclustering/` | UPGMA on shape-based distances |
| 2 | A | `16s/`, `barcodes/diversity.csv` | 3 h barcode point dropped for m1, m2, m4, m5 |
| | B | `jacobian/` | ridgelines, shared time grid across mice |
| | C | `jacobian/`, `16s/for_jacobian/` | warm-up rule; Hill ¹D from the curated taxon set |
| 3 | A | `jacobian/…matrices_by_time` | m1 spectrum; x = Re, y = Im |
| | B | `jacobian/…matrices_by_time` | every eigenvalue, all mice |
| | C–D | `jacobian/`, `16s/family/` | mean ± s.d. across mice |
| 4 | A | `16s/family/` (controls) | control clock |
| | B | `16s/diversity.csv` | controls against the colonised band |
| | C | `jacobian/` (controls) | one-sided Mann–Whitney |
| 5 | A | `genomics/rpsE_alignment.fasta` | the three-residue deletion |
| | B, C | — | ribosome render, plate photograph — no code |
| | D | `genomics/card_rgi_hits.csv` | CARD/RGI strict hits |

Panels **1A, 5B and 5C** have no code — a drawing, a structural render and a
photograph. They exist only inside the composites in `assets/`.

---

## Notes on the repository itself

- **The main-figure panels are committed** so that they are downloadable
  without running anything. EPS has no alpha channel, so semi-transparent fills
  flatten against white — **the PNG is the faithful rendering**. The densest
  panels (1C–D, S1, S2B) rasterise their barcode layer inside the EPS to avoid
  millions of line segments.
- **The supplementary panels ship as code, not as images.** `figures/supplementary/`
  holds one script per panel and `make supplementary` draws S1–S6; the rendered
  files are not tracked.
- **`pipeline/` cannot run from a clone.** It reads raw inputs that are not part
  of this repository. It is included so that the derived tables are auditable
  rather than taken on trust. One discrepancy is worth knowing: the pipeline
  computes an *expanding*-window Jacobian, whereas the published figures use the
  *sliding*-window estimator that `src/succession/jacobian.py` re-derives from
  the stored matrices.
- **Not shipped:** raw reads and assemblies (~3.7 GB), intermediate genomics
  output (~3.3 GB), and the per-barcode processing state.
- **Docstrings carry the caveats.** Where a choice is contestable, the module
  says so *and quantifies the alternative*.

## Licence

Code is released under the MIT licence (`LICENSE`); derived data and assets
under CC BY 4.0 (`LICENSE-data`).
