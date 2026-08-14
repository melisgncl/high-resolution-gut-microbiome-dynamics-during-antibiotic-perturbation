# pipeline/ — provenance only

**Nothing in this directory is needed to draw a figure.** The figures read the
derived tables in `../data/`, which are shipped. This directory is here so a
reader can see *how those tables were made*, and so the analysis is auditable
rather than taken on trust.

These scripts **cannot be run from a clone**. They read raw inputs that are not
in this repository:

- 16S and barcode amplicon read tables (`data/raw/…`)
- ~3.7 GB of Nanopore reads and assemblies
- ~3.3 GB of intermediate genomics output (Bakta, Panaroo, CARD/RGI, nucmer)

Raw sequencing data is deposited separately; see the manuscript's data
availability statement.

---

## Execution order

The modules below produced every table in `../data/`. They ship as-is: rerunning
them is not part of the release, and reworking them would risk silently changing
published numbers for no benefit. `../tests/` pins the figure layer to those
values.

The modules are sequential; each reads what the previous wrote. Every script
carries a header block naming its exact inputs and outputs.

| module | scripts | produces in `../data/` |
|---|---|---|
| `01_16s/` | `00_qc` → `01_import` → `02_process` → `03_filter` → `04_plot_community` → `05_diversity` | `16s/family/*.csv`, `16s/diversity.csv` |
| `02_barcode/` | `00_qc` → `01_import` → `02_process` → `03_dynamics` → `04_diversity` → `05_intersect` → `06_fitness` → `07_popgen` → `08_freq_spectrum` → `09_crosscohort` | `barcodes/trajectories.csv.gz` (see note), `barcodes/diversity.csv` |
| `03_clustering/` | `00_filter` → `01_hclust` → `02_select` → `03_cross_sample` → `04_cluster_similarity` | `barcodes/clusters/*.csv`, `coclustering/similarity_all_overlap.csv` |
| `04_coclustering/` | `01_coclustering` → `02_all_sbd` | `coclustering/*_sbd_*.csv`, `coclustering/all_sbd_dist.csv` |
| `05_jacobian/` | `01_jacobian` → `02_lag_time` → `03_stability` → `04_interactions` → `05_outliers` | `jacobian/*_jacobian_timeseries.csv`, `jacobian/*_jacobian_matrices_by_time.csv` |
| `genomics/` | see below | `genomics/rpsE_alignment.fasta`, `genomics/card_rgi_hits.csv` |

`barcodes/trajectories.csv.gz` is not written directly by any script: it is the
twelve per-sample `<nm>_processed.csv` files from `02_barcode/02_process.R`
concatenated and gzipped, keeping only the five columns Figure 1C–D and S1 need.
It replaces ~1.1 GB of per-barcode processing state with 36 MB.

`metadata/` holds the configuration every module reads via `here::here()`:
rarefaction thresholds, interpolation targets, and the barcode colour
assignments (`all_top_max2.csv`, `top_colors3.csv`) that give a lineage the same
colour in every mouse.

Key method choices, as recorded in the script headers:

- **16S**: rarefaction thresholds set per cohort by sequencing depth; relative
  abundance at family level.
- **Clustering**: distance = Pearson correlation on log10 frequency, average
  linkage; cutoff chosen where the scaled-distance and cluster-count curves
  cross.
- **Co-clustering**: Shape-Based Distance on z-normalised trajectories, average
  linkage.
- **Jacobian**: the pipeline computes an **expanding**-window Jacobian plus
  eigenvalues, KPCA and changepoints. The published figures use the
  **sliding**-window estimator, which is re-derived in
  `../src/succession/jacobian.py` from the stored per-window matrices. The two
  are not the same estimator — see that module's docstring.

## `genomics/`

Scripts for the *Paenibacillus macerans* isolate analysis: assembly
comparison, *rpsE* alignment and read-level confirmation, ANI, pan-genome, and
the survey of novel regions.

Roughly grouped:

- **rpsE** — `build_rpsE_msa_input.py`, `interpret_mafft_rpsE.py`,
  `rpsE_ecoli_alignment.py`, `fig1*_rpsE_msa*.py`, `map_rpsE_vs_I6.sh`,
  `extract_rpsE_depth.sh`, `fig5*_rpsE_pileup*.py`
- **genome comparison** — `run_fastani.sh`, `fig3_genome_comparison.py`,
  `phase4_annotate_novel_regions.py`, `phase4_blast_novel_regions.py`,
  `phase4_lookup_ICE_donors.py`
- **pan-genome and gene content** — `fig4_panaroo_heatmap.py`,
  `fig6_pangenome_functional.py`, `multicopy_gene_comparison.py`,
  `gene_classification.py`
- **targeted resistance mechanisms** — `phase33_16S_h34_fixed.py`,
  `phase33_34_rsmI_16S_h34.py`, `phase35_ribosomal_proteins.py`
- **16S context** — `00_paenibacillus_16S_dynamics.R`,
  `fig2_paenibacillaceae_RA.R`

### How the two shipped genomics tables were made

Both were **literals inside the manuscript's own figure scripts**, not files.
`data/genomics/rpsE_alignment.fasta` and `data/genomics/card_rgi_hits.csv` are
faithful extractions of the sequences and the 13 RGI hits hard-coded in the
published `31_figure5_wgs_rpse.py` and `32_figure_card_rgi.py`, moved out so
`../figures/figure5.py` only draws. The values are unchanged.

The scripts in this directory are what established those results upstream.

### Known discrepancy

The published Figure 5 describes the *rpsE* deletion as shared by isolates
**M5/M6/M7/M8**, and the shipped FASTA header preserves that wording. The
assembly QC in this project excluded **M7** (9× coverage, 7 contigs,
contaminated) and verified the identical 162 aa *rpsE* in **M5, M6 and M8**
only. The sequence is the M5/M6/M8 consensus; the header overstates what was
verified. Nothing downstream reads the header.
