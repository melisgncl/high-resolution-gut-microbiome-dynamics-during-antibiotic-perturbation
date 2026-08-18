# high-resolution-gut-microbiome-dynamics-during-antibiotic-perturbation

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

## Licence

Code is released under the MIT licence (`LICENSE`); derived data and assets
under CC BY 4.0 (`LICENSE-data`).
