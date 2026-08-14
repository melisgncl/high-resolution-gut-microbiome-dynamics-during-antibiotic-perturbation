# =============================================================================
# Title:   01_hclust.R — Hierarchical clustering + automatic cutoff detection
# Author:  Melis Gencel
# Date:    2026-03-20
# Input:   results/tables/clustering/<nm>/<nm>_filtered.csv (from 00_filter.R)
#          metadata/0_config.R
# Output:  results/tables/clustering/<nm>/<nm>_hclust_all_cutoffs.csv
#          results/tables/clustering/<nm>/<nm>_threshold_selection.csv
#          results/tables/clustering/cutoffs_summary.csv
#          results/figures/clustering/<nm>/<nm>_threshold_diagnostic.pdf/.png
#          results/figures/clustering/<nm>/<nm>_dist_heatmap.pdf/.png
# Notes:   Distance metric: Pearson correlation on log10(freq)
#          HC method: average linkage
#          Optimal cutoff: crossing of scaled distance and cluster-count curves.
#          If multiple crossings, picks the one with highest cluster count.
#          Override cutoffs via cutoff_override vector below (NA = auto).
#          Saves cutoffs_summary.csv for review across all 8 samples.
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

source(here::here("metadata/0_config.R"))


# =============================================================================
# LOESS HELPER (same as 02_select.R — needed for threshold selection curve)
# =============================================================================

adjust_span <- function(x, y, span = 0.2) {
  if (any(is.nan(x)) || any(is.infinite(x)) ||
      any(is.nan(y)) || any(is.infinite(y)))
    stop("NaN/Inf in x or y (adjust_span).")
  fit <- suppressWarnings(loess(log10(y + 1e-7) ~ x, span = span))
  comps <- c("enp", "s", "one.delta", "two.delta",
             "trace.hat", "divisor", "robust", "pars", "kd")
  if (any(is.na(fit[comps]))) {
    if (span < 1) return(adjust_span(x, y, span + 0.1))
    stop("adjust_span: failed to fit LOESS even at span = 1.")
  }
  fit
}


# =============================================================================
# CUTOFF OVERRIDES
# NA = auto-detect crossing point
# Set numeric to force a specific cutoff (e.g. after inspecting the diagnostic plot)
# =============================================================================

# Old manual cutoffs from archive 4_clustering_doblin/2_select_barcodes.R
# (process_clustering calls). Pinned to reproduce the original clustering exactly.
cutoff_override <- c(
  m1 = 0.59, m2 = 0.69, m3 = 0.51, m4 = 0.35,
  m5 = 0.33, m6 = 0.35, m7 = 0.24, m8 = 0.32
)


# =============================================================================
# CUTOFF DETECTION FUNCTIONS
# =============================================================================

find_optimal_cutoff <- function(df) {
  # df: threshold_selection tibble with cols: cutoff, dist_spread, n_clusters
  # Threshold axis is reversed (high = coarse, low = fine).
  # Returns the cutoff at which dist_small and scaled n_clusters cross,
  # choosing the crossing with the highest cluster count when multiple exist.

  df <- df %>%
    arrange(desc(cutoff)) %>%
    mutate(
      cluster_scaled = n_clusters * (max(dist_spread) / max(n_clusters)),
      diff           = dist_spread - cluster_scaled
    )

  crossings <- which(diff(sign(df$diff)) != 0)

  if (length(crossings) == 0) {
    warning("No crossing found — using minimum |diff|.")
    return(df$cutoff[which.min(abs(df$diff))])
  }

  # Pick the first downward crossing (dist_spread drops below cluster_scaled)
  # coming from coarse → fine, matching archive visual selection logic.
  # Downward = sign goes from +1 to -1, i.e. diff(sign) < 0
  sign_changes   <- diff(sign(df$diff))
  downward       <- which(sign_changes < 0)
  best           <- if (length(downward) > 0) downward[1] else crossings[1]

  # Linear interpolation for precise crossing value
  x1 <- df$cutoff[best];     y1 <- df$diff[best]
  x2 <- df$cutoff[best + 1]; y2 <- df$diff[best + 1]
  cutoff_opt <- x1 - y1 * (x2 - x1) / (y2 - y1)

  round(cutoff_opt, 3)
}

resolve_cutoff <- function(nm, threshold_df) {
  ov <- cutoff_override[[nm]]
  if (!is.na(ov)) {
    message("    cutoff: OVERRIDE = ", ov)
    return(list(value = ov, is_override = TRUE))
  }
  opt <- find_optimal_cutoff(threshold_df)
  message("    cutoff: auto = ", opt)
  list(value = opt, is_override = FALSE)
}


# =============================================================================
# THRESHOLD SELECTION CURVE
# Matches archive (plotHCQuantification + applyLOESS):
#   - Clusters with < min_members barcodes are dropped, UNLESS mean freq >= 1e-3
#   - Per surviving cluster: fit LOESS on pooled (time, freq) barcodes,
#     predict at 10x dense time grid (log10 scale)
#   - min pairwise Euclidean between cluster LOESS trajectories
#   - n_clusters = surviving cluster count (post-filter)
# =============================================================================

compute_threshold_selection <- function(filtered_wide, clust,
                                        min_members = 8,
                                        dominant_freq_threshold = 1e-3) {
  time_cols <- setdiff(names(filtered_wide), c("ID", "mean", "points"))

  # Long format: rows = barcode × timepoint
  long_df <- filtered_wide %>%
    select(ID, all_of(time_cols)) %>%
    pivot_longer(-ID, names_to = "time", values_to = "freq") %>%
    mutate(time = as.integer(time))

  # Dense prediction grid (10x resolution, matching archive)
  min_t <- min(long_df$time)
  max_t <- max(long_df$time)
  xx    <- seq(min_t, max_t, length.out = (max_t - min_t) * 10 + 1)

  cutoffs <- seq(from = 0.1, to = max(clust$height), by = 0.01)

  results <- lapply(cutoffs, function(h) {
    assignments <- cutree(clust, h = h)

    # Attach cluster assignment to long data
    assign_df <- data.frame(ID = filtered_wide$ID,
                            cluster = assignments,
                            stringsAsFactors = FALSE)
    ld <- left_join(long_df, assign_df, by = "ID")

    # --- min_members filter (archive: filterHC) ---
    # Keep clusters with >= min_members barcodes OR mean_freq >= dominant_freq_threshold
    cluster_stats <- ld %>%
      group_by(cluster) %>%
      summarise(
        n_bc      = n_distinct(ID),
        mean_freq = mean(freq, na.rm = TRUE),
        .groups   = "drop"
      )

    keep <- cluster_stats %>%
      filter(n_bc >= min_members | mean_freq >= dominant_freq_threshold) %>%
      pull(cluster)

    ld_filt <- ld %>% filter(cluster %in% keep)
    n_surviving <- length(keep)

    if (n_surviving < 2) {
      return(data.frame(cutoff = round(h, 3), n_clusters = n_surviving,
                        dist_spread = 0))
    }

    # --- LOESS per surviving cluster (archive: applyLOESS) ---
    loess_mat <- tryCatch({
      do.call(rbind, lapply(keep, function(k) {
        sub <- ld_filt %>% filter(cluster == k)
        fit <- tryCatch(
          adjust_span(sub$time, sub$freq, span = 0.2),
          error = function(e) NULL
        )
        if (is.null(fit)) rep(NA_real_, length(xx))
        else               predict(fit, newdata = data.frame(x = xx))
      }))
    }, error = function(e) NULL)

    dist_spread <- if (is.null(loess_mat) || nrow(loess_mat) < 2) {
      0
    } else {
      min(as.numeric(dist(loess_mat)), na.rm = TRUE)
    }

    data.frame(cutoff = round(h, 3), n_clusters = n_surviving,
               dist_spread = dist_spread)
  })

  bind_rows(results)
}


# =============================================================================
# PER-SAMPLE HIERARCHICAL CLUSTERING
# =============================================================================

run_hclust <- function(nm) {
  message("  ", nm, ":")

  tbl_dir <- here::here("results/tables/clustering", nm)
  fig_dir <- here::here("results/figures/clustering", nm)

  filtered <- read_csv(file.path(tbl_dir, paste0(nm, "_filtered.csv")),
                       show_col_types = FALSE)

  time_cols <- setdiff(names(filtered), c("ID", "mean", "points"))
  message("    ", nrow(filtered), " barcodes × ", length(time_cols), " timepoints")

  # Build log10 frequency matrix (zeros → NA)
  mat <- as.matrix(filtered[, time_cols])
  mat[mat == 0] <- NA
  mat_log <- log10(mat)

  # Pearson distance: 1 - r(log10 freq_i, log10 freq_j)
  distmat <- as.matrix(1 - cor(t(mat_log),
                                use    = "pairwise.complete.obs",
                                method = "pearson"))

  # Guard against non-finite values (e.g. barcodes with all-NA after log)
  if (any(!is.finite(distmat))) {
    n_bad <- sum(!is.finite(distmat))
    warning(nm, ": ", n_bad, " non-finite distance values replaced with 2")
    distmat[!is.finite(distmat)] <- 2
  }

  # Hierarchical clustering — average linkage
  clust <- hclust(as.dist(distmat), method = "average")

  # ---- Save all-cutoffs assignment table ----
  cutoff_seq <- seq(from = 0.1, to = max(clust$height), by = 0.01)
  cut_list   <- lapply(cutoff_seq, function(h) cutree(clust, h = h))
  cut_df     <- as.data.frame(do.call(cbind, cut_list))
  colnames(cut_df) <- as.character(round(cutoff_seq, 3))
  cut_df$ID        <- filtered$ID
  cut_df           <- cut_df[, c("ID", setdiff(names(cut_df), "ID"))]

  write_csv(cut_df, file.path(tbl_dir, paste0(nm, "_hclust_all_cutoffs.csv")))
  message("    Saved: ", file.path(tbl_dir, paste0(nm, "_hclust_all_cutoffs.csv")))

  # ---- Threshold selection curve ----
  message("    computing threshold selection curve ...")
  thresh_df <- compute_threshold_selection(filtered, clust)

  write_csv(thresh_df, file.path(tbl_dir, paste0(nm, "_threshold_selection.csv")))
  message("    Saved: ", file.path(tbl_dir, paste0(nm, "_threshold_selection.csv")))

  # ---- Resolve cutoff ----
  cutoff_res  <- resolve_cutoff(nm, thresh_df)
  final_cut   <- cutoff_res$value
  is_override <- cutoff_res$is_override
  n_at_cut    <- length(unique(cutree(clust, h = final_cut)))
  message("    n_clusters at cutoff: ", n_at_cut)

  # ---- Plot 1: Threshold diagnostic ----
  scale_fac <- max(thresh_df$dist_spread) / max(thresh_df$n_clusters)

  p_thresh <- ggplot(thresh_df, aes(x = cutoff)) +
    geom_line(aes(y = dist_spread),
              colour = "black", linewidth = 1.2) +
    geom_line(aes(y = n_clusters * scale_fac),
              colour = "#2b8cbe", linewidth = 1) +
    geom_vline(xintercept = final_cut,
               linetype = "dashed", colour = "darkred", linewidth = 0.8) +
    scale_x_reverse() +
    scale_y_continuous(
      name     = "Min distance between cluster means (log10)",
      sec.axis = sec_axis(~ . / scale_fac, name = "Number of clusters")
    ) +
    annotate("text",
             x = final_cut, y = max(thresh_df$dist_spread) * 0.95,
             label  = paste0("cutoff = ", final_cut,
                             if (is_override) " [override]" else " [auto]"),
             hjust  = 1.1, size = 3.5, colour = "darkred") +
    labs(title = paste0(nm, " — threshold selection"),
         x     = "Threshold (reversed: coarse → fine)") +
    theme_minimal(base_size = 12)

  save_fig_clust(p_thresh, nm, paste0(nm, "_threshold_diagnostic"), w = 8, h = 6)

  # ---- Plot 2: Distance heatmap (ggplot tile, top 200 barcodes) ----
  n_show  <- min(nrow(distmat), 200)
  ord_idx <- clust$order[seq_len(n_show)]
  dist_sub <- distmat[ord_idx, ord_idx]

  bc_labels <- as.character(seq_len(n_show))
  rownames(dist_sub) <- bc_labels
  colnames(dist_sub) <- bc_labels

  dist_long <- as.data.frame(dist_sub) %>%
    mutate(row = bc_labels) %>%
    pivot_longer(-row, names_to = "col", values_to = "dist") %>%
    mutate(row = factor(row, levels = bc_labels),
           col = factor(col, levels = bc_labels))

  p_heat <- ggplot(dist_long, aes(x = col, y = row, fill = dist)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 1, name = "1 - r") +
    labs(title = paste0(nm, " — Pearson distance (top ", n_show, " barcodes)"),
         x = NULL, y = NULL) +
    theme_minimal(base_size = 8) +
    theme(axis.text = element_blank(), axis.ticks = element_blank())

  save_fig_clust(p_heat, nm, paste0(nm, "_dist_heatmap"), w = 8, h = 7)

  list(
    auto_cutoff  = if (is_override) NA_real_ else final_cut,
    override     = if (is_override) final_cut else NA_real_,
    final_cutoff = final_cut,
    n_clusters   = n_at_cut
  )
}


# =============================================================================
# RUN ALL SAMPLES + SAVE CUTOFF SUMMARY
# =============================================================================

message("--- 01_hclust (clustering): running for all samples ---")

all_samples_bc <- c(cohort_colonized, cohort_colonized_2)

cutoff_rows <- lapply(all_samples_bc, function(nm) {
  res <- run_hclust(nm)
  tibble(
    sample       = nm,
    auto_cutoff  = res$auto_cutoff,
    override     = res$override,
    final_cutoff = res$final_cutoff,
    n_clusters   = res$n_clusters
  )
})

cutoff_summary <- bind_rows(cutoff_rows)

message("\n--- Cutoff summary ---")
print(cutoff_summary, n = Inf)

if (!dir.exists(path_tables_clust))
  dir.create(path_tables_clust, recursive = TRUE)

write_csv(cutoff_summary,
          file.path(path_tables_clust, "cutoffs_summary.csv"))
message("Saved: ", file.path(path_tables_clust, "cutoffs_summary.csv"))

message("--- 01_hclust (clustering): done ---")
