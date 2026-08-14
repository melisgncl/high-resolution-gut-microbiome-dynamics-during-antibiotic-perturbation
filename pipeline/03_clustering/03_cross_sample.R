# =============================================================================
# Title:   03_cross_sample.R — Cluster-of-clusters: cross-sample HC
# Author:  Melis Gencel
# Date:    2026-03-20
# Input:   results/tables/clustering/<nm>/<nm>_clustered_loess.csv (from 02_select.R)
#          results/tables/clustering/cutoffs_summary.csv
#          metadata/0_config.R
# Output:  results/figures/clustering/cross_colonized.pdf/.png
#          results/figures/clustering/cross_colonized_2.pdf/.png
#          results/figures/clustering/cross_all.pdf/.png
#          results/tables/clustering/cross_<group>_dist.csv
# Notes:   Pearson correlation distance between cluster LOESS trajectories
#          (log10 scale, pairwise.complete.obs), average-linkage HC.
#          Heatmap: ggplot2 tiles + ggdendro dendrogram + mouse colour strip.
#          Three runs: cohort 1 only, cohort 2 only, all 8 mice.
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(patchwork)
library(ggdendro)

source(here::here("metadata/0_config.R"))


# =============================================================================
# LOAD + COLLATE LOESS TRAJECTORIES
# =============================================================================

message("--- 03_cross_sample (clustering): loading LOESS trajectories ---")

load_loess_wide <- function(nm) {
  path <- here::here("results/tables/clustering", nm,
                     paste0(nm, "_clustered_loess.csv"))

  read_csv(path, show_col_types = FALSE) %>%
    filter(!is.na(loess_value)) %>%
    mutate(label = paste0(nm, ".C", cluster)) %>%
    select(time, label, loess_value) %>%
    pivot_wider(names_from = label, values_from = loess_value)
}

loess_list <- setNames(
  lapply(c(cohort_colonized, cohort_colonized_2), load_loess_wide),
  c(cohort_colonized, cohort_colonized_2)
)

# Full-join all samples in a group on time, drop time column
collate_group <- function(sample_ids) {
  dfs <- loess_list[sample_ids]
  combined <- Reduce(function(a, b) full_join(a, b, by = "time"), dfs)
  as.data.frame(combined %>% arrange(time) %>% select(-time))
}


# =============================================================================
# CROSS-SAMPLE HC + HEATMAP
# =============================================================================

# Colour strip: one colour per mouse, mapped from project palettes
mouse_col_map <- c(
  setNames(pal1, cohort_colonized),
  setNames(pal2, cohort_colonized_2)
)

plot_cross_heatmap <- function(mat, group_name, group_label) {
  message("  ", group_name, ": ", ncol(mat), " clusters")

  cluster_names <- colnames(mat)
  n             <- length(cluster_names)

  # --- Pearson distance + HC (exact archive methodology) ---
  distmat <- as.matrix(
    1 - cor(mat, use = "pairwise.complete.obs", method = "pearson")
  )
  distmat[!is.finite(distmat)] <- 2

  clust         <- hclust(as.dist(distmat), method = "average")
  ordered_names <- cluster_names[clust$order]

  # --- Save distance matrix CSV ---
  write_csv(
    tibble(cluster = rownames(distmat), as.data.frame(distmat)),
    here::here("results/tables/clustering",
               paste0("cross_", group_name, "_dist.csv"))
  )
  message("    Saved: cross_", group_name, "_dist.csv")

  # ------------------------------------------------------------------ #
  #  Panel 1: dendrogram (column ordering)                              #
  # ------------------------------------------------------------------ #
  dend_data <- ggdendro::dendro_data(clust, type = "rectangle")

  p_dend <- ggplot() +
    geom_segment(data = dend_data$segments,
                 aes(x = x, y = y, xend = xend, yend = yend),
                 linewidth = 0.4, colour = "grey30") +
    scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    ggdendro::theme_dendro() +
    theme(plot.margin = margin(4, 5, 0, 5))

  # ------------------------------------------------------------------ #
  #  Panel 2: mouse colour annotation strip                             #
  # ------------------------------------------------------------------ #
  annot_df <- data.frame(
    name  = factor(ordered_names, levels = ordered_names),
    mouse = sub("\\..*", "", ordered_names)
  )

  p_annot <- ggplot(annot_df, aes(x = name, y = 1, fill = mouse)) +
    geom_tile() +
    scale_fill_manual(values = mouse_col_map, guide = "none") +
    scale_x_discrete(expand = c(0, 0)) +
    theme_void() +
    theme(plot.margin = margin(0, 5, 0, 5))

  # ------------------------------------------------------------------ #
  #  Panel 3: heatmap tiles                                             #
  # ------------------------------------------------------------------ #
  dist_long <- as.data.frame(distmat) %>%
    mutate(row = rownames(distmat)) %>%
    pivot_longer(-row, names_to = "col", values_to = "dist") %>%
    mutate(
      row = factor(row, levels = rev(ordered_names)),
      col = factor(col, levels = ordered_names)
    )

  p_heat <- ggplot(dist_long, aes(x = col, y = row, fill = dist)) +
    geom_tile(colour = "white", linewidth = 0.1) +
    scale_fill_gradient2(
      low      = "blue", mid = "white", high = "red",
      midpoint = 1, limits = c(0, 2), name = "1 \u2212 r"
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 7) +
    theme(
      axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      axis.text.y  = element_text(size = 6),
      legend.position = "right",
      plot.margin  = margin(0, 5, 5, 5)
    )

  # ------------------------------------------------------------------ #
  #  Compose: dendrogram / colour strip / heatmap                       #
  # ------------------------------------------------------------------ #
  p_full <- (p_dend / p_annot / p_heat) +
    plot_layout(heights = c(1.5, 0.15, 5)) +
    plot_annotation(
      title    = paste0("Cluster-of-clusters \u2014 ", group_label),
      subtitle = "Pearson correlation between LOESS trajectories (log10)",
      theme    = theme(
        plot.title    = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8)
      )
    )

  out_dir_f <- here::here("results/figures/clustering")
  if (!dir.exists(out_dir_f)) dir.create(out_dir_f, recursive = TRUE)

  for (fmt in fig_formats) {
    out_path <- file.path(out_dir_f, paste0("cross_", group_name, ".", fmt))
    ggplot2::ggsave(out_path, plot = p_full,
                    width = 14, height = 13, dpi = fig_dpi,
                    limitsize = FALSE,
                    device    = if (fmt == "pdf") cairo_pdf else NULL)
    message("  Saved: ", out_path)
  }
}


# =============================================================================
# RUN THREE COMPARISONS
# =============================================================================

message("--- 03_cross_sample (clustering): computing cross-sample HC ---")

mat_col1 <- collate_group(cohort_colonized)
mat_col2 <- collate_group(cohort_colonized_2)
mat_all  <- collate_group(c(cohort_colonized, cohort_colonized_2))

plot_cross_heatmap(mat_col1, "colonized",   "Colonized cohort 1 (m1\u2013m4)")
plot_cross_heatmap(mat_col2, "colonized_2", "Colonized cohort 2 (m5\u2013m8)")
plot_cross_heatmap(mat_all,  "all",         "All mice (m1\u2013m8)")

message("--- 03_cross_sample (clustering): done ---")
