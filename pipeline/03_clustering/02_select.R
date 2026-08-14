# =============================================================================
# Title:   02_select.R — Apply cutoff, order clusters, fit LOESS, plot
# Author:  Melis Gencel
# Date:    2026-03-20
# Input:   results/tables/clustering/<nm>/<nm>_filtered.csv
#          results/tables/clustering/<nm>/<nm>_hclust_all_cutoffs.csv
#          results/tables/clustering/cutoffs_summary.csv
#          metadata/0_config.R
# Output:  results/tables/clustering/<nm>/<nm>_clustered_series.csv
#          results/tables/clustering/<nm>/<nm>_clustered_loess.csv
#          results/figures/clustering/<nm>/<nm>_cluster<k>.pdf/.png
#          results/figures/clustering/<nm>/<nm>_loess_overview.pdf/.png
# Notes:   Cluster ordering:
#            1) Clusters with non-zero mean at last timepoint → ranked by mean freq desc
#            2) Extinct clusters → ranked by time of disappearance desc
#          LOESS smoothing applied to log10(freq) per cluster (span auto-adjusted).
#          Per-cluster trajectory plots: individual barcode lines + LOESS overlay.
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

source(here::here("metadata/0_config.R"))


# =============================================================================
# HELPERS
# =============================================================================

adjust_span <- function(x, y, span = 0.2) {
  # Fit LOESS on log10(y + epsilon); widen span until stable or give up.
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

apply_loess_to_clusters <- function(long_df) {
  # long_df: tibble with cols time (int), frequency (num), cluster (int), ID
  # Returns: tibble with cols cluster, time (dense grid), loess_value (log10 scale)
  max_t <- max(long_df$time)
  min_t <- min(long_df$time)
  xx    <- seq(min_t, max_t, length.out = (max_t - min_t) * 10 + 1)

  long_df %>%
    group_by(cluster) %>%
    group_modify(~ {
      fit <- tryCatch(
        adjust_span(.x$time, .x$frequency, span = 0.2),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        tibble(time = xx, loess_value = NA_real_)
      } else {
        tibble(time = xx, loess_value = predict(fit, newdata = data.frame(x = xx)))
      }
    }) %>%
    ungroup()
}

order_clusters <- function(long_df) {
  # Returns a vector of cluster numbers in display order:
  # persistent clusters (mean freq > 0 at last TP) sorted desc by mean freq,
  # then extinct clusters sorted desc by time of last non-zero cumulative freq.
  last_t <- max(long_df$time)

  means_at_last <- long_df %>%
    filter(time == last_t) %>%
    group_by(cluster) %>%
    summarise(avg = mean(frequency, na.rm = TRUE), .groups = "drop")

  persistent   <- means_at_last %>% filter(avg >  0) %>% arrange(desc(avg))
  extinct_raw  <- means_at_last %>% filter(avg == 0)

  # For extinct clusters: last timepoint where cumulative freq > 0
  last_nonzero <- long_df %>%
    group_by(cluster, time) %>%
    summarise(cum = sum(frequency, na.rm = TRUE), .groups = "drop") %>%
    filter(cum > 0) %>%
    group_by(cluster) %>%
    summarise(last_tp = max(time), .groups = "drop")

  extinct <- extinct_raw %>%
    left_join(last_nonzero, by = "cluster") %>%
    arrange(desc(last_tp))

  c(persistent$cluster, extinct$cluster)
}


# =============================================================================
# PER-SAMPLE PROCESSING
# =============================================================================

process_sample <- function(nm, cutoff_summary) {
  message("  ", nm, ":")

  tbl_dir <- here::here("results/tables/clustering", nm)
  fig_dir <- here::here("results/figures/clustering", nm)

  final_cutoff <- cutoff_summary$final_cutoff[cutoff_summary$sample == nm]
  n_expected   <- cutoff_summary$n_clusters[cutoff_summary$sample == nm]
  message("    cutoff = ", final_cutoff, " | expected clusters = ", n_expected)

  # ---- Load data ----
  filtered <- read_csv(file.path(tbl_dir, paste0(nm, "_filtered.csv")),
                       show_col_types = FALSE)
  all_cuts <- read_csv(file.path(tbl_dir, paste0(nm, "_hclust_all_cutoffs.csv")),
                       show_col_types = FALSE)

  # Find nearest available cutoff column
  cut_cols  <- setdiff(names(all_cuts), "ID")
  cut_nums  <- as.numeric(cut_cols)
  nearest   <- cut_cols[which.min(abs(cut_nums - final_cutoff))]
  message("    using cutoff column: ", nearest)

  cluster_assignment <- all_cuts %>%
    select(ID, cluster = all_of(nearest))

  # ---- Build long series ----
  time_cols <- setdiff(names(filtered), c("ID", "mean", "points"))

  long_df <- filtered %>%
    select(ID, all_of(time_cols)) %>%
    pivot_longer(-ID, names_to = "time", values_to = "frequency") %>%
    mutate(time = as.integer(time)) %>%
    left_join(cluster_assignment, by = "ID") %>%
    filter(!is.na(cluster))

  # ---- min_members filter (matches threshold selection curve) ----
  # Keep clusters with >= 8 barcodes OR mean frequency >= 1e-3 (dominant small clusters)
  cluster_stats <- long_df %>%
    group_by(cluster) %>%
    summarise(
      n_bc      = n_distinct(ID),
      mean_freq = mean(frequency, na.rm = TRUE),
      .groups   = "drop"
    )
  keep_clusters <- cluster_stats %>%
    filter(n_bc >= 8 | mean_freq >= 1e-3) %>%
    pull(cluster)

  n_dropped <- n_distinct(long_df$cluster) - length(keep_clusters)
  if (n_dropped > 0)
    message("    dropped ", n_dropped, " clusters (< 8 members, mean freq < 1e-3)")

  long_df <- long_df %>% filter(cluster %in% keep_clusters)

  # Order clusters
  cluster_order <- order_clusters(long_df)
  long_df <- long_df %>%
    mutate(
      cluster_orig = cluster,
      cluster      = match(cluster, cluster_order)   # renumber 1..K in display order
    ) %>%
    filter(!is.na(cluster))

  n_clusters <- length(unique(long_df$cluster))
  message("    n_clusters (final): ", n_clusters)

  # ---- Save clustered series ----
  write_csv(long_df, file.path(tbl_dir, paste0(nm, "_clustered_series.csv")))
  message("    Saved: ", file.path(tbl_dir, paste0(nm, "_clustered_series.csv")))

  # ---- LOESS smoothing ----
  loess_df <- apply_loess_to_clusters(long_df)
  write_csv(loess_df, file.path(tbl_dir, paste0(nm, "_clustered_loess.csv")))
  message("    Saved: ", file.path(tbl_dir, paste0(nm, "_clustered_loess.csv")))

  # ---- Per-cluster trajectory plots ----
  eff_breaks <- sort(unique(long_df$time))
  eff_limits <- range(eff_breaks)

  cluster_plots <- vector("list", n_clusters)

  for (k in seq_len(n_clusters)) {
    sub_k  <- long_df  %>% filter(cluster == k)
    loes_k <- loess_df %>% filter(cluster == k)
    col_k  <- cluster.colors[k]

    y_min <- max(min(sub_k$frequency[sub_k$frequency > 0], na.rm = TRUE) * 0.5,
                 1e-7)

    p_k <- ggplot(sub_k, aes(x = time, y = frequency)) +
      geom_line(aes(group = ID), colour = col_k, alpha = 0.5, linewidth = 0.5) +
      geom_line(data = loes_k,
                aes(x = time, y = 10^loess_value),
                colour = "black", linewidth = 1.2, na.rm = TRUE) +
      scale_y_log10(limits = c(y_min, 1)) +
      scale_x_continuous(breaks = eff_breaks, limits = eff_limits) +
      labs(title  = paste0("C", k, "  (n=", length(unique(sub_k$ID)), ")"),
           x = "Time", y = "Frequency") +
      theme_Publication(base_size = 9, aspect.ratio = 0.75)

    save_fig_clust(p_k, nm, paste0(nm, "_cluster", k), w = 8, h = 6)
    cluster_plots[[k]] <- p_k
  }

  # ---- Combined all-clusters strip (one row per mouse) ----
  p_strip <- patchwork::wrap_plots(cluster_plots, nrow = 1) +
    patchwork::plot_annotation(
      title = paste0(nm, " — all clusters"),
      theme = theme(plot.title = element_text(size = 12, face = "bold"))
    )
  strip_w <- n_clusters * 4
  strip_h <- 4
  save_fig_clust(p_strip, nm, paste0(nm, "_clusters_strip"), w = strip_w, h = strip_h)

  # ---- Combined LOESS overview ----
  loess_plot_df <- loess_df %>%
    mutate(cluster = factor(paste0("C", cluster),
                            levels = paste0("C", seq_len(n_clusters))))

  eff_col <- cluster.colors[seq_len(n_clusters)]
  names(eff_col) <- paste0("C", seq_len(n_clusters))

  p_overview <- ggplot(loess_plot_df,
                       aes(x = time, y = 10^loess_value,
                           group = cluster, colour = cluster)) +
    geom_line(linewidth = 1, na.rm = TRUE) +
    scale_colour_manual(values = eff_col, name = "Cluster") +
    scale_y_log10() +
    scale_x_continuous(breaks = eff_breaks, limits = eff_limits) +
    labs(title = paste0(nm, " — cluster LOESS overview"),
         x = "Time index", y = "Mean frequency (log10)") +
    theme_Publication(base_size = 12, aspect.ratio = 0.75) +
    theme(legend.position = "right")

  save_fig_clust(p_overview, nm, paste0(nm, "_loess_overview"), w = 10, h = 6)

  invisible(NULL)
}


# =============================================================================
# RUN ALL SAMPLES
# =============================================================================

message("--- 02_select (clustering): applying cutoffs and plotting ---")

cutoff_summary <- read_csv(
  here::here("results/tables/clustering/cutoffs_summary.csv"),
  show_col_types = FALSE
)

all_samples_bc <- c(cohort_colonized, cohort_colonized_2)

for (nm in all_samples_bc) {
  process_sample(nm, cutoff_summary)
}

message("--- 02_select (clustering): done ---")
