# Shipped data

Every table here is *derived*. Raw reads, assemblies and the per-barcode
processing state are not part of this repository — see `../pipeline/` for the
code that produced these files from them.

Total ≈ 41 MB, of which 36 MB is the single barcode trajectory cache.

## The two clocks

Every table is indexed by a sampling index, and **the index does not mean the
same thing in the two groups**:

| group | mice | index means | range |
|---|---|---|---|
| colonised | `m1`–`m8` | sampling **slot** — 1 = 3 h, 2 = 6 h, 3 = 12 h, 4 = day 1, then slot *k* = day *k* − 3 | 1–19 (cohort 1), 1–18 (cohort 2) |
| control | `c_m1`–`c_m4` | the **day** itself | 1–10 |

Convert with `succession.timeaxis.to_days(index, group=...)`, which has no
default for `group` on purpose. Passing control indices through the colonised
converter misplaces every control point on every time axis, and — because the
Figure 4C comparison window was defined relative to the control span — silently
resizes that test as well.

---

## `16s/`

| file | contents | consumed by |
|---|---|---|
| `family/<mouse>_family.csv` | Family-level composition, long form. `Family`, `Time` (index), `Abundance` (reads), `Abundance.family` (relative abundance, 0–1), `Sample`. 49 families across the eight colonised mice, 134 across the four controls, 160 distinct in total; 20 in `m1`. Loaded by `io.load_family`, which keeps `Family`, `Time` and the relative abundance renamed to `abundance`. | 1G–H, 2A, 3C–D, 4A |
| `diversity.csv` | Hill numbers per sample from the **full** family table: `Time`, `Sample`, `q_0` (richness), `q_1` (exp Shannon), `q_inf`, `Cohort`. 166 rows. Cohort-1 mice have usable `q_1` at only 15 of 19 slots — the gaps are real and figures break the line rather than bridge them. | 2A, 4B |
| `for_jacobian/<mouse>_16S_taxa.csv` (colonised)<br>`for_jacobian/<mouse>_16S_family.csv` (controls) | Wide matrix, `Time` plus one column per taxon, relative abundance. Only the taxa that survived filtering into the interaction model — 5 for `m1`. Present at **every** slot, unlike `diversity.csv`. | 2C, and the state matrix behind every Jacobian |

`diversity.csv` and Hill q1 recomputed from `for_jacobian/` are **not** the same
series (mean difference 0.010, max 0.032 for `m1`), because they are computed
over different taxon sets. Figure 2C uses the recomputed one, via
`diversity.hill_q1_from_taxa`. Using the wrong source is enough to miss the
published per-mouse coefficients.

## `barcodes/`

| file | contents | consumed by |
|---|---|---|
| `trajectories.csv.gz` | **36 MB**, ~10.4 M rows: every barcode at every timepoint. `Sample`, `Time`, `ID`, `Freq` (frequency, 0–1), `hex_line` (stored colour; `#cccccc` = unassigned). | 1C–D, S1 |
| `diversity.csv` | Hill numbers over barcode lineages, same columns as the 16S table. 127 rows. | 2A |
| `clusters/<mouse>_loess_clusters.csv` | LOESS consensus per clonal cluster on the dense grid (0.1 slots). `cluster`, `time`, `loess_value` = **log10** frequency. | 1E–F, S2B, and the clone rows of every colonised state matrix |
| `clusters/<mouse>_clustered_series.csv` | The individual barcodes behind each cluster: `ID`, `time`, `frequency` (linear), `cluster`, `cluster_orig`. | S2B |
| `clusters/<mouse>_threshold_selection.csv` | Cluster count against correlation cutoff: `cutoff`, `n_clusters`, `dist_spread`. | S2A |

**Cluster numbering.** Clusters are meant to be numbered by mean frequency at
the end of the experiment. In `m2` and `m8` C1 and C2 are the wrong way round in
the stored files, and `io.load_clone_loess` swaps them (`config.CLONE_RELABEL`).
The swap is applied for display but *not* when building the Jacobian state, where
labels are only identifiers.

**Near-empty 3 h samples.** `m1`, `m2`, `m4`, `m5` have a 3 h barcode sample
built from a handful of reads (Hill q1 of 5.7–13.0 against 10⁴–10⁵ six hours
later). Figure 2A drops them; the Jacobian keeps them, because removing them
would start the dense grid a slot later and change every downstream evaluation
count. See `config.BARCODE_3H_EMPTY`.

## `cfu/`

`cfu_m1-m4.csv`, `cfu_m5-m8.csv` — wide, `Time` in **hours** (not slots, the one
exception in the repo) plus one column per mouse, *E. coli* CFU per gram. Zeros
mean below detection and are dropped at load. Consumed by 1B.

## `jacobian/`

| file | contents |
|---|---|
| `<mouse>_jacobian_timeseries.csv` | Long form: `time` (evaluation index), `window`, `effector_j` (driver), `target_i` (target), `strength`. |
| `<mouse>_jacobian_matrices_by_time.csv` | The same values as square matrices: `time`, `window`, `target_i`, then one column per driver. Rows = targets, columns = drivers. |

`strength(i, j) = cov(dz_i/dt, z_j) = J[i ← j]` — **i responds, j acts**. `z` is
log10 abundance. This is a covariance, not a normalised Jacobian; the textbook
`/ var(z_j)` denominator is dropped to match the published analysis.

The evaluation grid is not the list of sampling slots: 17 values for `m1`–`m4`
(2–17, then 18.9), 16 for `m5`, 15 for `m6`–`m8`. Slot 1 has no derivative, and
the last two slots are represented by one dense-grid endpoint. A **warm-up rule**
then discards evaluations whose window runs off the start of the series, which is
what sets the published n = 113 at window 5. Consumed by 2B, 2C, 3A–D, 4C, S4, S5.

## `coclustering/`

| file | contents | consumed by |
|---|---|---|
| `<mouse>_sbd_distance.csv` | Pairwise shape-based distance between every clone and family trajectory in one mouse: `series1`, `series2`, `dist`, `pair_id`. | 1I, S3 |
| `<mouse>_sbd_clusters.csv` | Cluster assignment per series: `series`, `sbd_cluster`. | S3 |
| `all_sbd_dist.csv` | The full 95 × 95 distance matrix pooled across mice, series labelled `<mouse>.<series>`. | S6 |
| `similarity_all_overlap.csv` | 72 × 72 barcode-overlap similarity between clonal clusters across mice. | S6 |

## `genomics/`

| file | contents | consumed by |
|---|---|---|
| `rpsE_alignment.fasta` | Two protein sequences: the *Paenibacillus macerans* I6 reference RpsE (165 aa) and the gut-isolate RpsE (162 aa). | 5A |
| `card_rgi_hits.csv` | 13 CARD/RGI strict hits: `gene`, `drug_class`, `mechanism`, `identity_pct`. | 5D |

The isolate protein is three residues shorter. Which three are "deleted" is
ambiguous — the reference reads `…I-N-R-V-A-K-V-V-K-G-G…` and removing either
`VAK` at 20–22 or `AKV` at 21–23 gives the identical isolate sequence, because
the region is a short repeat. The caption's "ΔVAK at 20–22" and the genomics
record's "Δ(A21-K22-V23)" describe the same event from opposite ends of that
ambiguity.

No aminoglycoside-modifying determinant appears among the 13 hits — no *aadA* /
ANT(9) anywhere in the genome — which is what rules out an acquired
spectinomycin cassette and leaves the chromosomal *rpsE* deletion.

> **Caveat on the FASTA header.** The isolate record is named
> `gut_isolates_M5_M6_M7_M8`, but the assembly QC behind this project excluded
> **M7** (9× coverage, 7 contigs, contaminated) and confirmed the identical
> *rpsE* in M5, M6 and M8 only. The header overstates what was verified; the
> sequence itself is the M5/M6/M8 consensus. Nothing downstream reads the header
> except the Figure 5A title, which prints "gut isolates".
