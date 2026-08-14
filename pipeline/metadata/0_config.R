# =============================================================================
# Title:   Global Configuration — HGT Study
# Author:  Melis Gencel
# Date:    2026-03-19
# Purpose: Single source of truth for all analysis modules (01_16S through
#          05_DCM). Source this file at the top of every analysis script via:
#            source(here::here("metadata/0_config.R"))
#          Requires the 'here' package to be loaded before sourcing.
# =============================================================================


# =============================================================================
# HOW TO ADAPT THIS PIPELINE FOR NEW SAMPLES
# =============================================================================
#
# This is the ONLY file you need to edit to add new samples.
# The clustering scripts (00_filter.R, 01_hclust.R, 02_select.R) will
# automatically pick up any samples listed in the cohort vectors below.
#
# STEP 1 — Add your sample ID to one of the cohort vectors, e.g.:
#   cohort_colonized_2 <- c("m5", "m6", "m7", "m8", "m10")
#
# STEP 2 — Add your raw barcode file to data/raw/barcodes/
#   Expected format: tab-separated, columns ID | <timepoint> | <timepoint> ...
#   File must be named consistently (see bc_file_map below).
#
# STEP 3 — Add the file → sample ID mapping in bc_file_map (Section 3).
#
# STEP 4 — Run the scripts in order:
#   Rscript analysis/03_clustering/00_filter.R
#   Rscript analysis/03_clustering/01_hclust.R
#   Rscript analysis/03_clustering/02_select.R
#
# Optional overrides (leave as NA to use automatic detection):
#   - time_threshold_override in 00_filter.R  (min timepoints with Freq > 0)
#   - cutoff_override         in 01_hclust.R  (HC distance cutoff)
#   Both fall back to sensible defaults for any sample not listed.
#
# =============================================================================


# =============================================================================
# SECTION 1 — SAMPLE MANIFEST
# =============================================================================

# Mapping from original file names to standardised sample IDs used throughout
# this project. m1-m4 = colonised cohort 1 (files M5-M8),
# m5-m8 = colonized cohort 2, separate facility (files NM5-NM8).
sample_map <- c(
  "M5"  = "m1", "M6"  = "m2", "M7"  = "m3", "M8"  = "m4",
  "NM5" = "m5", "NM6" = "m6", "NM7" = "m7", "NM8" = "m8"
)

cohort_colonized   <- c("m1", "m2", "m3", "m4")   # Colonized cohort 1 (E. coli + intact microbiome)
cohort_colonized_2 <- c("m5", "m6", "m7", "m8")   # Colonized cohort 2 (E. coli, separate facility)
cohort_controls    <- c("c_m1", "c_m2", "c_m3", "c_m4")  # SPN antibiotic controls
all_samples        <- c(cohort_colonized, cohort_colonized_2, cohort_controls)


# =============================================================================
# SECTION 1b — CLONAL CLUSTER RELABELLING
# =============================================================================

# Doblin numbers clonal clusters by its own internal ordering, which does not
# always match the convention used in the manuscript: C1 is the cluster with the
# highest mean frequency at the end of the experiment, C2 the next, and so on.
#
# In m2 and m8 the stored C1 and C2 are the wrong way round by that convention
# (end-of-experiment frequency ratios C1/C2 = 0.58 and 0.09 respectively), so
# they are swapped here. m4 also has C2 ending above C1 (ratio 0.55) but its
# labelling is correct and must be left alone.
#
# The swap is applied when tables are read, not by rewriting them. Stored
# outputs under results/tables/ keep the original Doblin numbering; anything
# reported or plotted goes through the helpers below. Use the relabelled
# clusters everywhere clonal clusters are referenced.
clone_relabel <- list(
  m2 = c("1" = 2L, "2" = 1L),
  m8 = c("1" = 2L, "2" = 1L)
)

# Remap integer cluster numbers. `cluster` and `mouse` are parallel vectors,
# so this works on a table spanning several mice.
relabel_clone <- function(cluster, mouse) {
  out <- as.integer(cluster)
  for (nm in intersect(unique(mouse), names(clone_relabel))) {
    map <- clone_relabel[[nm]]
    idx <- which(mouse == nm)
    key <- as.character(out[idx])
    hit <- key %in% names(map)
    out[idx[hit]] <- unname(map[key[hit]])
  }
  out
}

# Remap "C1"/"C2" style labels, and the "m2.C1" form used by the pooled
# co-clustering tables. Labels that are not clusters (family names) pass through
# untouched. `mouse` may be omitted when the label carries the mouse prefix.
relabel_clone_label <- function(label, mouse = NULL) {
  label <- as.character(label)
  prefixed <- grepl("^(m[0-9]+)\\.C[0-9]+$", label)
  bare     <- grepl("^C[0-9]+$", label)

  if (any(prefixed)) {
    nm  <- sub("\\..*$", "", label[prefixed])
    num <- as.integer(sub("^.*\\.C", "", label[prefixed]))
    label[prefixed] <- paste0(nm, ".C", relabel_clone(num, nm))
  }
  if (any(bare)) {
    if (is.null(mouse)) stop("relabel_clone_label(): bare 'C<n>' labels need `mouse`.")
    nm  <- if (length(mouse) == 1) rep(mouse, sum(bare)) else mouse[bare]
    num <- as.integer(sub("^C", "", label[bare]))
    label[bare] <- paste0("C", relabel_clone(num, nm))
  }
  label
}


# =============================================================================
# SECTION 2 — DATA & OUTPUT PATHS
# =============================================================================

path_16S_raw      <- here::here("data/raw/16S")
path_16S_control  <- here::here("data/raw/16S/control_SPN.tsv")
path_barcodes     <- here::here("data/raw/barcodes")
path_metadata     <- here::here("metadata")
path_results_16S  <- here::here("results/figures/16S")
path_results_bc   <- here::here("results/figures/barcodes")
path_results_clust <- here::here("results/figures/clustering")
path_tables       <- here::here("results/tables")
path_tables_clust <- here::here("results/tables/clustering")

# Figure output settings
fig_formats      <- c("pdf", "png")
fig_dpi          <- 300
fig_width_single <- 8.25   # single panel, inches
fig_width_panel4 <- 25     # 4-panel row, inches
fig_height       <- 6      # standard height, inches


# =============================================================================
# SECTION 3 — RAREFACTION THRESHOLDS
# =============================================================================

# Threshold set per cohort based on minimum sample sequencing depth.
# Do not reduce higher-depth cohorts to match lower-depth ones.
rarefaction_thresholds <- c(
  m1   = 90000, m2   = 90000, m3   = 90000, m4   = 90000,
  m5   = 60000, m6   = 60000, m7   = 60000, m8   = 60000,
  c_m1 = 60000, c_m2 = 60000, c_m3 = 60000, c_m4 = 60000
)


# =============================================================================
# SECTION 4 — PER-MOUSE INTERPOLATION TARGETS
# =============================================================================

# BIOLOGICAL DECISION — do not change without checking raw data.
# Missing timepoints occurred during the family expansion window.
# Families listed here are interpolated linearly across missing timepoints.
# Empty string means no interpolation applied for that sample.
interpolation_targets <- list(
  m1   = "Paenibacillaceae",                          # expansion window gap
  m2   = "Paenibacillaceae",                          # expansion window gap
  m3   = "Ruminococcaceae",                           # expansion window gap
  m4   = "",
  m5   = "",
  m6   = "Lactobacillaceae",                          # expansion window gap
  m7   = c("Lactobacillaceae", "Paenibacillaceae"),   # expansion window gap
  m8   = "Lactobacillaceae",                          # expansion window gap
  c_m1 = "",
  c_m2 = "",
  c_m3 = "",
  c_m4 = ""
)


# =============================================================================
# SECTION 5 — TIME AXIS SETTINGS
# =============================================================================

# Named lists replace the original hash::hash() objects.
# Access is identical: breaks.hash[["m"]], labels.hash[["c2"]], etc.
# Keys: "m"  = colonised cohort 1 (m1-m4)
#        "c2" = colonized cohort 2 (m5-m8)
#        "c"  = controls (c_m1-c_m4)

breaks.hash <- list(
  m  = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19),
  c2 = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18),
  c  = c(1,2,3,4,5,6,7,8,9,10)
)

labels.hash <- list(
  m  = c("3h","","12h","","2d","","","","6d","","","","10d","","","","14d","","16d"),
  c2 = c("3h","","12h","","2d","","","","6d","","","","10d","","","","14d",""),
  c  = c("1d","","3d","","5d","","7d","","9d","")
)

limits.hash <- list(
  m  = c(1, 19),
  c2 = c(1, 18),
  c  = c(1, 10)
)


# =============================================================================
# SECTION 6 — PLOT THEMES
# =============================================================================

theme_Publication <- function(base_size = 32, base_family = "Helvetica",
                               aspect.ratio = 0.75) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size = base_size, base_family = base_family)
    + theme(
      plot.title   = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      text         = element_text(),
      panel.background = element_rect(colour = NA, fill = "#FCFCFC"),
      plot.background  = element_rect(colour = NA),
      panel.border = element_rect(fill = NA, colour = "black", size = 1.5),
      axis.title   = element_text(size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text    = element_text(),
      axis.line    = element_line(colour = "black"),
      axis.ticks   = element_line(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin  = unit(c(10, 5, 5, 5), "mm"),
      aspect.ratio = aspect.ratio
    ))
}

theme_Publication_noYaxis <- function(base_size = 32, base_family = "Helvetica",
                                       aspect.ratio = 0.75) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size = base_size, base_family = base_family)
    + theme(
      plot.title   = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      text         = element_text(),
      panel.background = element_rect(colour = NA, fill = "#FCFCFC"),
      plot.background  = element_rect(colour = NA),
      panel.border = element_rect(fill = NA, colour = "black", size = 1.5),
      axis.title   = element_text(size = rel(1)),
      axis.title.x = element_text(vjust = -0.2),
      axis.title.y = element_blank(),
      axis.text    = element_text(),
      axis.line    = element_line(colour = "black"),
      axis.ticks   = element_line(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin  = unit(c(10, 5, 5, 5), "mm"),
      aspect.ratio = aspect.ratio
    ))
}

theme_Publication_bottomLegend <- function(base_size = 32, base_family = "Helvetica",
                                            aspect.ratio = 0.75) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size = base_size, base_family = base_family)
    + theme(
      plot.title   = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      text         = element_text(),
      panel.background = element_rect(colour = NA, fill = "#FCFCFC"),
      plot.background  = element_rect(colour = NA),
      panel.border = element_rect(fill = NA, colour = "black", size = 1.5),
      axis.title   = element_text(size = rel(1)),
      axis.title.x = element_text(vjust = -0.2),
      axis.title.y = element_blank(),
      axis.text    = element_text(),
      axis.line    = element_line(colour = "black"),
      axis.ticks   = element_line(),
      legend.position   = c(0.35, 0.15),
      legend.direction  = "horizontal",
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      plot.margin  = unit(c(10, 5, 5, 5), "mm"),
      aspect.ratio = aspect.ratio
    ))
}

theme_Matrix <- function(base_size = 32, base_family = "Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size = base_size, base_family = base_family)
    + theme(
      plot.title   = element_text(size = rel(1.2), hjust = 0.5),
      text         = element_text(),
      panel.background  = element_rect(colour = NA, fill = "#FCFCFC"),
      plot.background   = element_rect(colour = NA),
      panel.border = element_rect(fill = NA, colour = "black", size = 1.5),
      axis.title   = element_text(size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text.x  = element_text(angle = 45, hjust = 1),
      axis.text    = element_text(),
      axis.line    = element_line(colour = "black"),
      axis.ticks   = element_line(),
      legend.title       = element_text(angle = -90),
      legend.title.align = 0.5,
      panel.grid.major   = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.margin  = unit(c(10, 5, 5, 5), "mm")
    ))
}


# =============================================================================
# SECTION 7 — COHORT COLOR PALETTES
# =============================================================================

pal1 <- c("#045a8d", "#2b8cbe", "#74a9cf", "#bdc9e1")  # colonized_1 m1-m4 (blues)
pal2 <- c("#7a0177", "#c51b8a", "#9138A7", "#552586")  # colonized_2 m5-m8 (purples)
pal3 <- c("#E8A20C", "#FFAF19", "#E2522F", "#CB6015")  # controls  c_m1-c_m4 (oranges)

c1 <- "#045a8d"   # cohort 1 single colour
c2 <- "#552586"   # cohort 2 single colour
c3 <- "#CB6015"   # control single colour


# =============================================================================
# SECTION 8 — TAXONOMIC COLOR PALETTES
# =============================================================================

# --- Genus ---
# Note: empty-string entry removed from original (was a silent bug).
genus.colors <- c(
  "Escherichia/Shigella"           = "#003D18",
  "Desulfovibrio"                  = "#00A341",
  "Paenibacillus"                  = "#561AA3",
  "Lachnospiraceae_NK4A136_group"  = "#063D74",
  "Shuttleworthia"                 = "#060874",
  "Acetatifactor"                  = "#0A66C2",
  "Lachnospiraceae_UCG-001"        = "#1685F3",
  "Clostridium_sensu_stricto_1"    = "#51A3F6",
  "Ruminiclostridium_9"            = "#08519C",
  "Romboutsia"                     = "#77B8F8",
  "Bacteroides"                    = "#BD0026",
  "Prevotellaceae_UCG-001"         = "#FD8D3C",
  "Alloprevotella"                 = "#FEB24C",
  "Alistipes"                      = "#FED976",
  "Methanobacterium"               = "#a6a600",
  "other"                          = "black",
  "Parasutterella"                 = "#571150",
  "Ruminococcaceae_UCG-005"        = "#ab1e96",
  "Lachnoclostridium"              = "#8bbada",
  "Ruminiclostridium_6"            = "#48cb49",
  "Lactobacillus"                  = "#2ec7f4",
  "Akkermansia"                    = "#a37110",
  "Anaeroplasma"                   = "#fee24b",
  "Veillonella"                    = "#2a01ae",
  "Oscillibacter"                  = "#e17d68"
)

# --- Family ---
# Canonical version from 0_config/0_config.R. Do not change colors inline
# in analysis scripts — edit here only.
family.colors <- c(
  "Enterobacteriaceae"                      = "#003D18",
  "Moraxellaceae"                           = "#006D2C",
  "Desulfovibrionaceae"                     = "#00A341",
  "Lachnospiraceae"                         = "#063D74",
  "Ruminococcaceae"                         = "#08519C",
  "Peptostreptococcaceae"                   = "#0A66C2",
  "Clostridiaceae_1"                        = "#77B8F8",
  "Clostridiales_vadinBB60_group"           = "#51A3F6",
  "Peptococcaceae"                          = "#1685F3",
  "Paenibacillaceae"                        = "#561AA3",
  "Lactobacillaceae"                        = "#7E1CFF",
  "Erysipelotrichaceae"                     = "#a6007a",
  "Muribaculaceae"                          = "#9D0211",
  "Bacteroidaceae"                          = "#BD0026",
  "Rikenellaceae"                           = "#F03B20",
  "Prevotellaceae"                          = "#FD8D3C",
  "Marinifilaceae"                          = "#FEB24C",
  "Tannerellaceae"                          = "#FED976",
  "Anaeroplasmataceae"                      = "#633228",
  "Deferribacteraceae"                      = "#ff028d",
  "Saccharimonadaceae"                      = "#A18E66",
  "Xiphinematobacteraceae"                  = "#CEA2FD",
  "Methanobacteriaceae"                     = "#a6a600",
  "other"                                   = "black",
  "Sutterellaceae"                          = "#5F021F",
  "Acholeplasmataceae"                      = "#b21561",
  "Oscillospiraceae"                        = "#e17d68",
  "Akkermansiaceae"                         = "#30bc9e",
  "[Eubacterium] coprostanoligenes group"   = "#aa862f",
  "Enterococcaceae"                         = "#a9cea8",
  "Bacillaceae"                             = "#8f7e83",
  "Corynebacteriaceae"                      = "#1f194d",
  "Microbacteriaceae"                       = "#8bbada",
  "Erwiniaceae"                             = "#d4e125",
  "Veillonellaceae"                         = "#2a01ae",
  "Comamonadaceae"                          = "#b087b7",
  "Xanthomonadaceae"                        = "#640b66",
  "Leptotrichiaceae"                        = "#be6a5b",
  "Pseudomonadaceae"                        = "#29e2b2",
  "Sphingomonadaceae"                       = "#fffc78",
  "Weeksellaceae"                           = "#527d54"
)

# Append rare-family colors loaded from CSV
bottom_genera  <- readr::read_csv(here::here("metadata/bottom_genera_color.csv"),
                                  show_col_types = FALSE)
family_bottom  <- as.character(bottom_genera$hex)
names(family_bottom) <- as.character(bottom_genera$all.main.families)
family.colors  <- c(family.colors, family_bottom)

# --- Phylum ---
phylum.colors <- c(
  "Bacteroidetes"  = "#BD0026",
  "Firmicutes"     = "#0A66C2",
  "Proteobacteria" = "#00A341",
  "Deferribacteres"= "#ff028d",
  "Patescibacteria"= "#A18E66",
  "Tenericutes"    = "#633228"
)


# =============================================================================
# SECTION 9 — ORDERED TAXONOMIC VECTORS
# =============================================================================

# Use these to set factor levels so taxa appear in consistent order across plots.
major.family <- c(
  "Enterobacteriaceae", "Moraxellaceae", "Desulfovibrionaceae",
  "Lachnospiraceae", "Ruminococcaceae", "Peptostreptococcaceae",
  "Clostridiaceae_1", "Clostridiales_vadinBB60_group", "Peptococcaceae",
  "Paenibacillaceae", "Lactobacillaceae", "Erysipelotrichaceae",
  "Muribaculaceae", "Bacteroidaceae", "Rikenellaceae", "Prevotellaceae",
  "Marinifilaceae", "Tannerellaceae", "Anaeroplasmataceae",
  "Deferribacteraceae", "Saccharimonadaceae", "Xiphinematobacteraceae",
  "Methanobacteriaceae", "other", "Acholeplasmataceae", "Oscillospiraceae",
  "Akkermansiaceae", "[Eubacterium] coprostanoligenes group",
  "Enterococcaceae", "Bacillaceae", "Corynebacteriaceae",
  "Microbacteriaceae", "Erwiniaceae", "Veillonellaceae", "Sutterellaceae",
  "Comamonadaceae", "Xanthomonadaceae", "Leptotrichiaceae",
  "Pseudomonadaceae", "Sphingomonadaceae", "Weeksellaceae"
)

major.genus <- c(
  "Escherichia/Shigella", "Desulfovibrio", "Paenibacillus",
  "Lachnospiraceae_NK4A136_group", "Acetatifactor", "Lachnospiraceae_UCG-001",
  "Clostridium_sensu_stricto_1", "Bacteroides", "Ruminiclostridium_9",
  "Shuttleworthia", "Prevotellaceae_UCG-001", "Alloprevotella",
  "Alistipes", "Romboutsia", "Methanobacterium", "other",
  "Parasutterella", "Ruminiclostridium_6", "Lachnoclostridium",
  "Ruminococcaceae_UCG-005", "Lactobacillus", "Akkermansia",
  "Anaeroplasma", "Veillonella", "Oscillibacter"
)

major.phylum <- c(
  "Proteobacteria", "Firmicutes", "Bacteroidetes",
  "Patescibacteria", "Tenericutes", "Deferribacteres"
)


# =============================================================================
# SECTION 10 — BARCODE COLORS
# =============================================================================

all.top.max <- readr::read_csv(here::here("metadata/all_top_max2.csv"),
                               show_col_types = FALSE)

top_colors2 <- readr::read_csv(here::here("metadata/top_colors3.csv"),
                                show_col_types = FALSE)

# Long shuffled color list for low-frequency barcodes.
# Repetition is intentional — we just need every barcode colored.
# set.seed ensures the rare-barcode color assignment is reproducible.
long.color.list        <- rep(top_colors2$hex, 50)
set.seed(42)
long.color.list.random <- sample(long.color.list)


# =============================================================================
# SECTION 10b — BARCODE FILE MAP & NOISE FILTERS
# =============================================================================

# Maps each project sample ID to its raw file name and cohort subfolder.
# Colonized 1 (m1-m4): data/raw/barcodes/processed_sample/M/M{5-8}/
# Colonized 2 (m5-m8): data/raw/barcodes/processed_sample/NM/M{5-8}/
bc_file_map <- list(
  m1 = list(cohort = "M",  file = "M5"),
  m2 = list(cohort = "M",  file = "M6"),
  m3 = list(cohort = "M",  file = "M7"),
  m4 = list(cohort = "M",  file = "M8"),
  m5 = list(cohort = "NM", file = "M5"),
  m6 = list(cohort = "NM", file = "M6"),
  m7 = list(cohort = "NM", file = "M7"),
  m8 = list(cohort = "NM", file = "M8")
)

# Barcode IDs confirmed as technical noise (sequencing artefacts).
# Source: visual inspection of raw clustering data; IDs cluster at unexpected
# high frequency with no biological plausible explanation.
# NULL means no ID filter for that sample.
bc_noise_ids <- list(
  m1 = NULL, m2 = NULL, m3 = NULL, m4 = NULL,
  m5 = c(68808L, 68805L, 68804L, 66523L, 68793L, 68794L, 68795L,
         68796L, 68797L, 68798L, 68812L),
  m6 = NULL,
  m7 = c(55725L),
  m8 = NULL
)

# Timepoints removed due to poor coverage / technical failure.
# NULL means no timepoint filter for that sample.
bc_noise_times <- list(
  m1 = NULL, m2 = NULL, m3 = NULL, m4 = NULL,
  m5 = NULL,
  m6 = 14L,
  m7 = 14L,
  m8 = NULL
)

# Barcode-specific output paths
path_bc_m   <- here::here("data/raw/barcodes/processed_sample/M")
path_bc_nm  <- here::here("data/raw/barcodes/processed_sample/NM")


# =============================================================================
# SECTION 11 — OTHER PALETTES & CLUSTER LABELS
# =============================================================================

# Colors for DCM eigenvalue plots
eigen.colors <- c(
  '#9b2335', '#0F4C81', '#e15d44', '#6667AB', '#aebf9a', '#c5a2e8',
  '#e1bc44', '#884c5e', '#218EAE', '#009B77', '#fabebe', '#008080',
  '#e6beff', '#9a6324', '#C96115', '#800000', '#808000', '#ffd8b1',
  '#00756d', "#B0BA97", "#E940E2", "#69A1AD", "#E447AD", "#84F046",
  "#ADEAA5", "#E2ECBC", "#BC31EB", "#E29B3F", "#904D7E", "#C756D7",
  "#B5DFE9", "#D8537A", "#DE898C", "#60ABE4", "#9CE2CA", "#8A764C",
  "#E582CF", "#E7BBAC", "#ABE36B"
)

# Colors for barcode cluster plots
cluster.colors <- c(
  "#3cb44b", "#4363d8", "#e6194B", "#e8ca00", "#911eb4", "#f58231",
  "#22766dff", "#42d4f4", "#f032e6", "#9A6324", "#2F4C39", "#1D3F6E",
  "#94170f", "#665679", "#F17829", "#97A69C", "#606EA9", "#A9606E",
  "#A99060", "#F8F1AE", "#bcf60c", "#fabebe", "#008080", "#e6beff",
  "#9a6324", "#fffac8", "#003D18", "#003D18", "#82a479", "#c74b0e",
  "#77b5fe", "#ccff00"
)

# Cluster label placeholders used in barcode/co-clustering modules
cluster_labels <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10")


# =============================================================================
# SECTION 12 — UTILITY FUNCTIONS
# =============================================================================

# Captures a base-graphics plot as a grob for compositing with ggplot.
# Used in co-clustering and DCM heatmap figures.
grab_grob <- function() {
  grid.echo()
  grid.grab()
}

# Save a ggplot to PDF and PNG (both formats) in path_results_16S.
# w/h default to fig_width_panel4 / fig_height from Section 2.
# cairo_pdf is used for PDF to support Unicode and custom fonts.
save_fig <- function(p, filename_stem, w = fig_width_panel4, h = fig_height) {
  for (fmt in fig_formats) {
    path <- file.path(path_results_16S, paste0(filename_stem, ".", fmt))
    ggplot2::ggsave(path, plot = p, width = w, height = h,
                    dpi    = fig_dpi,
                    device = if (fmt == "pdf") cairo_pdf else NULL)
    message("  Saved: ", path)
  }
}

# Save a clustering figure to the per-mouse subfolder under path_results_clust.
# Creates the directory if it does not exist.
save_fig_clust <- function(p, nm, filename_stem, w = 10, h = 8) {
  out_dir <- file.path(path_results_clust, nm)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  for (fmt in fig_formats) {
    path <- file.path(out_dir, paste0(filename_stem, ".", fmt))
    ggplot2::ggsave(path, plot = p, width = w, height = h,
                    dpi       = fig_dpi,
                    limitsize = FALSE,
                    device    = if (fmt == "pdf") cairo_pdf else NULL)
    message("  Saved: ", path)
  }
}

# Identical to save_fig() but writes to path_results_bc (barcode module).
save_fig_bc <- function(p, filename_stem, w = fig_width_panel4, h = fig_height) {
  if (!dir.exists(path_results_bc)) dir.create(path_results_bc, recursive = TRUE)
  for (fmt in fig_formats) {
    path <- file.path(path_results_bc, paste0(filename_stem, ".", fmt))
    ggplot2::ggsave(path, plot = p, width = w, height = h,
                    dpi    = fig_dpi,
                    device = if (fmt == "pdf") cairo_pdf else NULL)
    message("  Saved: ", path)
  }
}

# Save a DCM figure to results/figures/DCM/ (or a per-mouse subfolder).
# nm = NULL → saves to flat DCM dir; nm = "<mouse>" → saves to DCM/<nm>/ subdir.
save_fig_dcm <- function(p, filename_stem, nm = NULL, w = 10, h = 7) {
  base <- here::here("results/figures/DCM")
  out_dir <- if (is.null(nm)) base else file.path(base, nm)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  for (fmt in fig_formats) {
    path <- file.path(out_dir, paste0(filename_stem, ".", fmt))
    ggplot2::ggsave(path, plot = p, width = w, height = h,
                    dpi       = fig_dpi,
                    limitsize = FALSE,
                    device    = if (fmt == "pdf") cairo_pdf else NULL)
    message("  Saved: ", path)
  }
}


# =============================================================================
# SECTION 13 — DCM PARAMETERS
# =============================================================================

# Expanding-window Jacobian step size (dense-grid units).
# With a ~170-point dense grid, window=10 yields ~17 Jacobian matrices/mouse.
dcm_window      <- 10L

# LOESS span for smoothing 16S family series to the barcode dense grid.
dcm_loess_span  <- 0.2

# Maximum zero timepoints a 16S family may have before being dropped.
# Applied to colonised mice.  Controls use dcm_mean_filter instead.
dcm_zero_filter <- 8L

# Mean-abundance floor for control 16S families (archive: filter_columns_by_mean).
dcm_mean_filter <- 5e-3

# Minimum number of mice in which a Jacobian pair must be an outlier
# (in both raw-strength AND delta-strength) to be classified as "common".
dcm_min_mice    <- 4L

# Focal taxon for stability contribution and lag-time analyses.
dcm_hgt_target  <- "Paenibacillaceae"

# Density-plot color palette for Jacobian interaction strength (before/after HGT).
# Exact archive palette — do not change.
dcm_density_palette_anchors <- c(
  "#11695a", "#1b9e77",   # Teal
  "#a74701", "#d95f02",   # Orange
  "#524d94", "#7570b3",   # Violet
  "#b01f6b", "#e7298a",   # Magenta
  "#447713", "#66a61e",   # Olive
  "#b18301", "#e6ab02"    # Yellow-gold
)

# Maps each sample to the time-axis key used in breaks/labels/limits hashes.
# "m"  → cohort 1 (m1-m4),  "c2" → cohort 2 (m5-m8),  "c" → controls
dcm_time_type <- c(
  setNames(rep("m",  length(cohort_colonized)),   cohort_colonized),
  setNames(rep("c2", length(cohort_colonized_2)), cohort_colonized_2),
  setNames(rep("c",  length(cohort_controls)),    cohort_controls)
)
