# =============================================================================
# Title:   01_coclustering.R — Per-mouse SBD co-clustering of barcode + 16S
# Author:  Melis Gencel
# Date:    2026-03-20
# Input:   results/tables/clustering/<nm>/<nm>_clustered_loess.csv
#          results/tables/16S/family/<nm>_family.csv
#          metadata/0_config.R
# Output:  results/tables/coclustering/<nm>/<nm>_sbd_distance.csv
#          results/tables/coclustering/<nm>/<nm>_sbd_clusters.csv
#          results/figures/coclustering/<nm>/<nm>_sbd_dendrogram.png
#          results/figures/coclustering/<nm>/<nm>_sbd_pvclust.png
# Notes:   Exact archive methodology (Gencel 2020):
#          - Distance: Shape-Based Distance (SBD), z-normalized
#          - HC: average linkage
#          - 16S preprocessing: drop families with >= 8 zero timepoints,
#            log10(x + 1e-6), reinterpolate to barcode LOESS grid length
#          - Dendrogram: barcode cluster labels coloured by cluster.colors[k],
#            16S family labels in black (same as archive)
#          - pvclust: 1000 bootstrap replicates, SBD distance
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(dtwclust)    # registers SBD method in proxy
library(proxy)
library(pvclust)
library(dendextend)

source(here::here("metadata/0_config.R"))


# =============================================================================
# DATA LOADING
# =============================================================================

# Barcode LOESS → wide matrix: rows = dense grid, cols = C1, C2, ...
load_loess_wide <- function(nm) {
  df <- read_csv(
    here::here("results/tables/clustering", nm,
               paste0(nm, "_clustered_loess.csv")),
    show_col_types = FALSE
  ) %>%
    filter(!is.na(loess_value)) %>%
    pivot_wider(id_cols = time, names_from = cluster,
                values_from = loess_value) %>%
    arrange(time)

  mat           <- as.matrix(df[, -1])
  colnames(mat) <- paste0("C", colnames(mat))
  mat
}

# 16S family abundances → wide matrix: rows = timepoints, cols = families
# Uses the curated per-mouse taxa files from the archive 6_coclustering workflow
# (data/16S/taxa_long_format/<nm>_16S_taxa.tsv) so co-clustering reproduces the
# original results. These are already wide (Time + one column per family).
load_16S_wide <- function(nm) {
  df <- read_tsv(
    here::here("data/16S/taxa_long_format", paste0(nm, "_16S_taxa.tsv")),
    show_col_types = FALSE
  ) %>%
    arrange(Time)

  as.matrix(df[, setdiff(names(df), "Time")])
}


# =============================================================================
# 16S PREPROCESSING  (archive: transformTaxa)
# =============================================================================

# Filter: keep families with < 8 zero timepoints (present at most timepoints)
# Transform: log10(x + 1e-6)
# Reinterpolate: stretch to target_len to match barcode dense grid
transform_taxa <- function(taxa_mat, target_len) {
  keep     <- apply(taxa_mat, 2, function(col) sum(col == 0) < 8)
  taxa_mat <- taxa_mat[, keep, drop = FALSE]
  ts_list  <- dtwclust::tslist(t(log10(taxa_mat + 1e-6)))
  dtwclust::reinterpolate(ts_list, new.length = target_len)
}


# =============================================================================
# DENDROGRAM COLOURING  (exact archive scheme)
# =============================================================================
# Barcode cluster labels (^C[0-9]+$ or \.C[0-9]+$) → cluster.colors[k]
# 16S family labels                                 → black

color_dend <- function(dend) {
  leaf_lab <- labels(dend)
  col_vec  <- setNames(rep("black", length(leaf_lab)), leaf_lab)

  bc_idx <- grep("(\\.C|^C)[0-9]+$", leaf_lab)
  for (i in bc_idx) {
    k <- as.integer(regmatches(leaf_lab[i],
                               regexpr("[0-9]+$", leaf_lab[i])))
    if (!is.na(k) && k >= 1L && k <= length(cluster.colors))
      col_vec[i] <- cluster.colors[k]
  }
  labels_colors(dend) <- col_vec
  dend
}


# =============================================================================
# PER-MOUSE CO-CLUSTERING
# =============================================================================

run_coclustering <- function(nm) {
  message("  ", nm, ": loading series ...")

  clust_mat  <- load_loess_wide(nm)
  taxa_mat   <- load_16S_wide(nm)
  target_len <- nrow(clust_mat)

  taxa_series  <- transform_taxa(taxa_mat, target_len)
  clust_series <- dtwclust::tslist(t(clust_mat))
  series       <- append(clust_series, taxa_series)

  message("    ", length(clust_series), " clusters  +  ",
          length(taxa_series), " families  =  ",
          length(series), " series total")

  # ---- SBD distance + average-linkage HC ----
  message("  ", nm, ": computing SBD distance ...")
  sbd_dist <- proxy::dist(series, method = "SBD", znorm = TRUE)
  attr(sbd_dist, "Labels") <- names(series)

  hc   <- hclust(sbd_dist, method = "average")
  dend <- color_dend(as.dendrogram(hc))

  # ---- Output directories ----
  out_tbl <- here::here("results/tables/coclustering", nm)
  out_fig <- here::here("results/figures/coclustering", nm)
  dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

  # ---- Save cophenetic distance (lower triangle) ----
  cop_mat <- as.matrix(cophenetic(hc))
  n_cop   <- nrow(cop_mat)
  idx     <- which(lower.tri(matrix(0, n_cop, n_cop)), arr.ind = TRUE)
  write_csv(
    data.frame(
      series1 = rownames(cop_mat)[idx[, 1]],
      series2 = colnames(cop_mat)[idx[, 2]],
      dist    = cop_mat[idx]
    ) %>% mutate(pair_id = paste(series1, series2, sep = "_")),
    file.path(out_tbl, paste0(nm, "_sbd_distance.csv"))
  )

  # ---- Save cluster assignments ----
  # Cut at the archive's per-mouse k (auto-selected by dtwclust in 6_coclustering)
  # so assignments reproduce the original; fall back to min(8, n-1) otherwise.
  coclust_k <- c(m1 = 4L, m2 = 3L, m3 = 9L, m4 = 4L,
                 m5 = 2L, m6 = 2L, m7 = 3L, m8 = 2L)
  k_cut <- if (!is.na(coclust_k[nm])) coclust_k[[nm]] else min(8L, length(series) - 1L)
  write_csv(
    data.frame(series      = hc$labels,
               sbd_cluster = cutree(hc, k = k_cut)),
    file.path(out_tbl, paste0(nm, "_sbd_clusters.csv"))
  )
  message("    Saved tables.")

  # ---- Coloured horizontal dendrogram (base R, matches archive) ----
  png(file.path(out_fig, paste0(nm, "_sbd_dendrogram.png")),
      res = 300, width = 8, height = 8, units = "in")
  par(mar = c(4, 2, 2, 8))
  plot(dend, horiz = TRUE,
       main = paste0(nm, " \u2014 SBD co-clustering"),
       xlab = "Height")
  dev.off()
  message("    Saved: ", nm, "_sbd_dendrogram.png")

  # ---- pvclust (1000 bootstrap replicates) ----
  message("  ", nm, ": pvclust (1000 boots) ...")

  series_mat           <- do.call(cbind, series)
  colnames(series_mat) <- names(series)

  sbd_fn <- function(x)
    proxy::dist(dtwclust::tslist(t(x)), method = "SBD", znorm = TRUE)

  set.seed(456L)
  pv <- pvclust::pvclust(series_mat,
                         method.hclust = "average",
                         method.dist   = sbd_fn,
                         nboot         = 1000,
                         quiet         = TRUE)

  png(file.path(out_fig, paste0(nm, "_sbd_pvclust.png")),
      res = 300, width = 12, height = 6, units = "in")
  par(mar = c(10, 2, 2, 0))
  plot(pv, hang = -1, cex = 0.6,
       main = paste0(nm, " \u2014 pvclust (SBD, 1000 boots)"))
  pvclust::pvrect(pv, alpha = 0.95)
  dev.off()
  message("    Saved: ", nm, "_sbd_pvclust.png")

  message("  ", nm, ": done.")
  invisible(list(hc = hc, series = series))
}


# =============================================================================
# RUN ALL MICE
# =============================================================================

message("--- 01_coclustering: per-mouse SBD co-clustering ---")

all_samples_bc <- c(cohort_colonized, cohort_colonized_2)

for (nm in all_samples_bc) {
  tryCatch(
    run_coclustering(nm),
    error = function(e)
      message("  ERROR in ", nm, ": ", conditionMessage(e))
  )
}

message("--- 01_coclustering: done ---")
