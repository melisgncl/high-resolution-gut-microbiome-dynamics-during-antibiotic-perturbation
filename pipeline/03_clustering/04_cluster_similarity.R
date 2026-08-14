# =============================================================================
# Title:   04_cluster_similarity.R — Barcode overlap and cluster correlation
# Author:  Melis Gencel
# Date:    2026-03-20
# Input:   results/tables/clustering/<nm>/<nm>_clustered_series.csv
#          results/tables/clustering/cross_<group>_dist.csv
#          data/Gavage_index_overlap.csv
#          metadata/0_config.R
# Output:  results/figures/clustering/similarity_<group>.pdf/.png
#          results/tables/clustering/similarity_<group>_overlap.csv
#          results/tables/clustering/similarity_<group>_z.csv
#          results/tables/clustering/similarity_<group>_scatter.csv
# Notes:   Exact archive methodology (Gauthier & Gencel 2020):
#          - Overlap metric: Szymkiewicz-Simpson coefficient
#            overlapCoef(s1, s2) = |s1 ∩ s2| / min(|s1|, |s2|)
#          - Null distribution: 1000 permutations per cluster pair,
#            sampling from each mouse's full barcode library
#          - Z-score: (observed - mean(null)) / sd(null)
#          - P-value: two-sided, -log10(2 * pnorm(-|z|))
#          - Scatter: overlap coefficient vs Pearson correlation (from step 3)
#          Three runs: cohort 1, cohort 2, all mice.
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(patchwork)

source(here::here("metadata/0_config.R"))


# =============================================================================
# OVERLAP AND Z-SCORE FUNCTIONS  (exact archive methodology)
# =============================================================================

# Szymkiewicz-Simpson overlap coefficient
overlapCoef <- function(s1, s2) {
  i <- length(intersect(s1, s2))
  u <- min(length(s1), length(s2))
  if (u == 0) return(0)
  i / u
}

# Single permutation: draw random samples of the observed sizes from libraries
sampleSimul <- function(s1_size, s2_size, lib_x, lib_y) {
  s1_size <- min(s1_size, length(lib_x))   # guard against edge cases
  s2_size <- min(s2_size, length(lib_y))
  xs <- sample(lib_x, size = s1_size, replace = FALSE)
  ys <- sample(lib_y, size = s2_size, replace = FALSE)
  m  <- min(length(xs), length(ys))
  if (m == 0) return(0)
  length(intersect(xs, ys)) / m
}

# Empirical Z-score: observed overlap vs null distribution of 1000 permutations
empiricalZ <- function(s1, s2, lib_x, lib_y, n_trials = 1000) {
  obs       <- overlapCoef(s1, s2)
  null_dist <- replicate(n_trials,
                         sampleSimul(length(s1), length(s2), lib_x, lib_y))
  mu    <- mean(null_dist)
  sigma <- sd(null_dist)
  if (sigma == 0) return(0)
  (obs - mu) / sigma
}


# =============================================================================
# DATA LOADING
# =============================================================================

message("--- 04_cluster_similarity: loading data ---")

gavage         <- read_csv(here::here("data/Gavage_index_overlap.csv"),
                           show_col_types = FALSE)
all_samples_bc <- c(cohort_colonized, cohort_colonized_2)

series_list <- setNames(
  lapply(all_samples_bc, function(nm) {
    read_csv(
      here::here("results/tables/clustering", nm,
                 paste0(nm, "_clustered_series.csv")),
      show_col_types = FALSE
    )
  }),
  all_samples_bc
)


# =============================================================================
# BUILD PER-MOUSE STRUCTURES
# =============================================================================

# Join clustered_series IDs to universal Center barcodes via Gavage index.
# Returns data.frame with columns: cluster, Center
merge_format <- function(series, nm) {
  uniq   <- unique(series[, c("ID", "cluster")])
  joined <- merge(uniq, gavage, by.x = "ID", by.y = paste0(nm, ".ID"))
  joined[, c("cluster", "Center")]
}

# Returns named list: cluster_number -> vector of Center IDs
split_clusters <- function(center_df) {
  lst <- split(center_df$Center, f = center_df$cluster)
  lapply(lst, function(x) x[!is.na(x)])
}

# Library: all unique Center IDs for a mouse across all clusters
build_library <- function(series, nm) {
  all_ids <- unique(series$ID)
  joined  <- merge(data.frame(id = all_ids), gavage,
                   by.x = "id", by.y = paste0(nm, ".ID"))
  joined$Center[!is.na(joined$Center)]
}

mouse_centers   <- setNames(
  lapply(all_samples_bc, function(nm) merge_format(series_list[[nm]], nm)),
  all_samples_bc
)
mouse_hash      <- setNames(
  lapply(all_samples_bc, function(nm) split_clusters(mouse_centers[[nm]])),
  all_samples_bc
)
mouse_libraries <- setNames(
  lapply(all_samples_bc, function(nm) build_library(series_list[[nm]], nm)),
  all_samples_bc
)


# =============================================================================
# MATRIX COMPUTATION
# =============================================================================

# Ordered flat list of (mouse, cluster_key) entries for a group
build_entries <- function(sample_ids, hash) {
  entries <- list()
  for (nm in sample_ids) {
    for (k in names(hash[[nm]])) {
      entries <- c(entries, list(list(nm = nm, k = k)))
    }
  }
  entries
}

compute_overlap_matrix <- function(entries, hash) {
  n      <- length(entries)
  labels <- sapply(entries, function(e) paste0(e$nm, ".C", e$k))
  mat    <- matrix(NA_real_, n, n, dimnames = list(labels, labels))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      s1        <- hash[[entries[[i]]$nm]][[entries[[i]]$k]]
      s2        <- hash[[entries[[j]]$nm]][[entries[[j]]$k]]
      mat[i, j] <- overlapCoef(s1, s2)
    }
  }
  mat
}

compute_z_matrix <- function(entries, hash, libraries, n_trials = 1000) {
  n      <- length(entries)
  labels <- sapply(entries, function(e) paste0(e$nm, ".C", e$k))
  mat    <- matrix(NA_real_, n, n, dimnames = list(labels, labels))
  total  <- n * n
  done   <- 0L
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      s1        <- hash[[entries[[i]]$nm]][[entries[[i]]$k]]
      s2        <- hash[[entries[[j]]$nm]][[entries[[j]]$k]]
      lib_x     <- libraries[[entries[[i]]$nm]]
      lib_y     <- libraries[[entries[[j]]$nm]]
      mat[i, j] <- tryCatch(
        empiricalZ(s1, s2, lib_x, lib_y, n_trials),
        error = function(e) NA_real_
      )
      done <- done + 1L
      if (done %% 200L == 0L)
        message("    Z-matrix: ", done, " / ", total, " pairs")
    }
  }
  mat
}


# =============================================================================
# OVERLAP HEATMAP + SCATTER PLOT
# =============================================================================

plot_similarity <- function(sample_ids, group_name, group_label,
                            hash, libraries, dist_csv_path) {

  message("  ", group_name, ": building cluster entries ...")
  entries <- build_entries(sample_ids, hash)
  n_pairs <- length(entries)
  message("    clusters: ", n_pairs)

  message("  ", group_name, ": computing overlap matrix ...")
  overlap_mat <- compute_overlap_matrix(entries, hash)

  message("  ", group_name, ": computing Z-score matrix (1000 permutations) ...")
  z_mat    <- compute_z_matrix(entries, hash, libraries, n_trials = 1000)
  pval_mat <- -log10(2 * pnorm(-abs(z_mat)))
  pval_mat[!is.finite(pval_mat)] <- 0

  # --- Load Pearson distance from step 3 → convert to correlation ---
  dist_df   <- read_csv(dist_csv_path, show_col_types = FALSE)
  dist_mat  <- as.matrix(dist_df[, -1])
  rownames(dist_mat) <- dist_df$cluster
  colnames(dist_mat) <- dist_df$cluster
  cor_mat   <- (dist_mat - 1) * -1   # r = 1 - dist

  # Align to common cluster labels (same order as dist_mat)
  common <- intersect(colnames(dist_mat), colnames(overlap_mat))
  if (length(common) < 2) {
    warning("  ", group_name, ": too few common labels — skipping.")
    return(invisible(NULL))
  }
  overlap_mat <- overlap_mat[common, common]
  z_mat       <- z_mat[common, common]
  pval_mat    <- pval_mat[common, common]
  cor_mat     <- cor_mat[common, common]

  # --- Save matrices ---
  out_tbl <- here::here("results/tables/clustering")
  write_csv(
    tibble(cluster = rownames(overlap_mat), as.data.frame(overlap_mat)),
    file.path(out_tbl, paste0("similarity_", group_name, "_overlap.csv"))
  )
  write_csv(
    tibble(cluster = rownames(z_mat), as.data.frame(z_mat)),
    file.path(out_tbl, paste0("similarity_", group_name, "_z.csv"))
  )
  message("    Saved overlap and Z matrices.")

  # --- Overlap heatmap (ordered by HC from Pearson distance, step 3) ---
  clust_ord     <- hclust(as.dist(dist_mat[common, common]), method = "average")
  ordered_names <- common[clust_ord$order]

  overlap_long <- as.data.frame(overlap_mat) %>%
    mutate(row = rownames(overlap_mat)) %>%
    pivot_longer(-row, names_to = "col", values_to = "overlap") %>%
    mutate(
      row = factor(row, levels = rev(ordered_names)),
      col = factor(col, levels = ordered_names)
    )

  p_overlap <- ggplot(overlap_long, aes(x = col, y = row, fill = overlap)) +
    geom_tile(colour = "white", linewidth = 0.1) +
    scale_fill_gradient(low = "white", high = "blue",
                        limits = c(0, 1), name = "Overlap\ncoef.") +
    labs(
      title    = paste0("Barcode overlap \u2014 ", group_label),
      subtitle = "Szymkiewicz\u2013Simpson overlap coefficient",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 7) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
      axis.text.y = element_text(size = 5),
      legend.position = "right"
    )

  # --- Scatter: overlap coefficient vs Pearson correlation ---
  nc      <- length(common)
  tri_idx <- which(lower.tri(matrix(0, nc, nc)), arr.ind = TRUE)

  scatter_df <- data.frame(
    label1      = common[tri_idx[, 1]],
    label2      = common[tri_idx[, 2]],
    overlap     = overlap_mat[tri_idx],
    z_score     = z_mat[tri_idx],
    neg_log_p   = pval_mat[tri_idx],
    correlation = cor_mat[tri_idx]
  ) %>%
    mutate(significant = neg_log_p >= -log10(0.05))

  write_csv(scatter_df,
            file.path(out_tbl,
                      paste0("similarity_", group_name, "_scatter.csv")))

  sig_df    <- filter(scatter_df, significant)
  nonsig_df <- filter(scatter_df, !significant)

  p_scatter <- ggplot() +
    geom_point(data = nonsig_df,
               aes(x = overlap, y = correlation),
               colour = "#6a6a6a", shape = 1, size = 3, stroke = 1) +
    geom_point(data = sig_df,
               aes(x = overlap, y = correlation, size = neg_log_p),
               colour = "blue", shape = 1, stroke = 1.2) +
    scale_size_continuous(range = c(3, 8), name = "-log10(p)") +
    labs(
      title    = paste0("Cluster similarity \u2014 ", group_label),
      subtitle = "Overlap coefficient vs LOESS trajectory correlation",
      x        = "Overlap coefficient",
      y        = "Cluster correlation (Pearson r)"
    ) +
    xlim(0, 1) +
    theme_Publication(base_size = 12, aspect.ratio = 1)

  # --- Compose: overlap heatmap | scatter ---
  p_full <- p_overlap + p_scatter +
    plot_layout(widths = c(2, 1)) +
    plot_annotation(
      theme = theme(plot.background = element_rect(fill = "white", colour = NA))
    )

  out_fig <- here::here("results/figures/clustering")
  for (fmt in fig_formats) {
    out_path <- file.path(out_fig,
                          paste0("similarity_", group_name, ".", fmt))
    ggplot2::ggsave(out_path, plot = p_full,
                    width = 20, height = 10, dpi = fig_dpi,
                    limitsize = FALSE,
                    device    = if (fmt == "pdf") cairo_pdf else NULL)
    message("  Saved: ", out_path)
  }
}


# =============================================================================
# RUN THREE GROUPS
# =============================================================================

message("--- 04_cluster_similarity: running analyses ---")

plot_similarity(
  cohort_colonized, "colonized", "Cohort 1 (m1\u2013m4)",
  mouse_hash, mouse_libraries,
  here::here("results/tables/clustering/cross_colonized_dist.csv")
)

plot_similarity(
  cohort_colonized_2, "colonized_2", "Cohort 2 (m5\u2013m8)",
  mouse_hash, mouse_libraries,
  here::here("results/tables/clustering/cross_colonized_2_dist.csv")
)

plot_similarity(
  c(cohort_colonized, cohort_colonized_2), "all", "All mice (m1\u2013m8)",
  mouse_hash, mouse_libraries,
  here::here("results/tables/clustering/cross_all_dist.csv")
)

message("--- 04_cluster_similarity: done ---")
