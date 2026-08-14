# =============================================================================
# Title:   05_outliers.R — Cross-mouse common outlier Jacobian pairs
# Author:  Melis Gencel
# Date:    2026-03-20
# Input:   results/tables/DCM/<nm>/<nm>_interactions.csv   (from 04)
#          results/tables/DCM/<nm>/<nm>_changepoints_geo.csv (from 01)
#          metadata/0_config.R
# Output:  results/tables/DCM/<nm>/<nm>_outliers_raw.csv
#          results/tables/DCM/<nm>/<nm>_outliers_delta.csv
#          results/tables/DCM/common_outlier_pairs.csv
#          results/figures/DCM/<nm>/<nm>_outliers.png/pdf
#          results/figures/DCM/common_outlier_pairs.png/pdf
# Notes:   Exact archive methodology (Gencel 2020):
#          - Outlier criteria (both applied per time window):
#              1. Raw strength: boxplot.stats(strength)$out (Tukey 1.5×IQR fence)
#              2. Delta strength: boxplot.stats(delta_strength)$out
#          - Common pair: flagged in both raw AND delta in >= dcm_min_mice mice.
#          - Tile plot + ridge plot per mouse with common pairs highlighted.
#          - Bar plot: cross-mouse pair count faceted by outlier type.
#          - Pair colors: deterministic random hex (set.seed(123)).
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggridges)
library(patchwork)

source(here::here("metadata/0_config.R"))


# =============================================================================
# OUTLIER DETECTION
# =============================================================================

# Tukey 1.5×IQR fence applied per time window (window column).
# Returns rows flagged as outliers.
find_raw_outliers <- function(df) {
  df %>%
    group_by(window) %>%
    filter(strength %in% boxplot.stats(strength)$out) %>%
    ungroup()
}

# Compute time-to-time delta in strength for each directional pair,
# then apply Tukey fence per window on the deltas.
find_delta_outliers <- function(df) {
  df %>%
    arrange(pair_directional, window) %>%
    group_by(pair_directional) %>%
    mutate(
      lag_strength   = lag(strength),
      delta_strength = strength - lag_strength
    ) %>%
    ungroup() %>%
    group_by(window) %>%
    filter(delta_strength %in% boxplot.stats(delta_strength)$out) %>%
    ungroup()
}

# Run both outlier methods for one mouse.
# Returns list: out_raw, out_delta, common_pairs
analyze_outliers <- function(df, nm) {
  df <- df %>%
    mutate(pair_directional = paste0(species_i, "\u2192", species_j),
           mouse = nm)

  out_raw   <- find_raw_outliers(df)
  out_delta <- find_delta_outliers(df)

  common_pairs <- intersect(unique(out_raw$pair_directional),
                             unique(out_delta$pair_directional))

  list(out_raw = out_raw, out_delta = out_delta, common_pairs = common_pairs,
       df_full = df)
}


# =============================================================================
# PER-MOUSE PLOTS (tile + ridge)
# =============================================================================

plot_outliers_mouse <- function(df_full, out_raw, common_pairs_used,
                                 pair_colors, phase_lines, nm) {
  df_flagged <- df_full %>%
    mutate(is_common = pair_directional %in% common_pairs_used)

  quantile_df <- df_full %>%
    group_by(window) %>%
    summarise(q0 = min(strength), q100 = max(strength), .groups = "drop")

  # Tile plot
  p_tile <- ggplot(df_flagged, aes(x = window, y = strength)) +
    geom_tile(aes(fill = ifelse(is_common, pair_directional, NA_character_)),
              width = 0.9, height = 0.05, colour = "white") +
    scale_fill_manual(
      values   = pair_colors,
      na.value = "darkgrey",
      name     = "Common outlier pair"
    ) +
    geom_line(data = quantile_df, aes(x = window, y = q0),
              colour = "black", linetype = "dashed", linewidth = 0.8) +
    geom_line(data = quantile_df, aes(x = window, y = q100),
              colour = "black", linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = phase_lines, linetype = "dashed",
               colour = "black", linewidth = 0.4) +
    labs(title = paste0(nm, " \u2014 common outlier pairs"),
         x = "Time window", y = "Strength") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "right")

  # Ridge plot
  common_df <- df_flagged %>% filter(is_common)
  palette_used <- pair_colors[names(pair_colors) %in%
                                unique(common_df$pair_directional)]

  p_ridge <- ggplot(df_flagged,
                    aes(x = strength, y = as.factor(window))) +
    geom_density_ridges(scale = 4, rel_min_height = 0.01,
                         colour = "grey30", fill = "grey80", alpha = 0.6) +
    {if (nrow(common_df) > 0)
      geom_point(data = common_df,
                 aes(colour = pair_directional),
                 position = position_jitter(height = 0.1),
                 size = 1.5, alpha = 0.8)} +
    {if (length(palette_used) > 0)
      scale_colour_manual(values = palette_used, name = "Common pair")} +
    labs(title = "Density ridges + common outliers",
         x = "Strength", y = "Time window") +
    theme_ridges(grid = FALSE, center_axis_labels = TRUE)

  p_tile + p_ridge
}


# =============================================================================
# PER-MOUSE RUNNER
# =============================================================================

run_outliers <- function(nm) {
  message("  ", nm, ": loading interactions ...")

  int_path <- here::here("results/tables/DCM", nm,
                           paste0(nm, "_interactions.csv"))
  if (!file.exists(int_path)) {
    message("    interactions.csv not found — skipping.")
    return(list(out_raw = NULL, out_delta = NULL))
  }
  df <- read_csv(int_path, show_col_types = FALSE)

  cpt_path <- here::here("results/tables/DCM", nm,
                           paste0(nm, "_changepoints_geo.csv"))
  phase_lines <- if (file.exists(cpt_path)) {
    cpt <- read_csv(cpt_path, show_col_types = FALSE)
    if (nrow(cpt) > 0) cpt$changepoint[cpt$score >= 0.8] else integer(0)
  } else integer(0)

  message("  ", nm, ": detecting outliers ...")
  res <- analyze_outliers(df, nm)

  out_tbl <- here::here("results/tables/DCM", nm)
  write_csv(res$out_raw,
            file.path(out_tbl, paste0(nm, "_outliers_raw.csv")))
  write_csv(res$out_delta,
            file.path(out_tbl, paste0(nm, "_outliers_delta.csv")))
  message("    Saved: ", nm, "_outliers_raw.csv / _delta.csv")

  list(out_raw = res$out_raw, out_delta = res$out_delta,
       df_full = res$df_full, phase_lines = phase_lines)
}


# =============================================================================
# CROSS-MOUSE AGGREGATION
# =============================================================================

# Deterministic random colors for all directional pairs (archive: set.seed(123))
generate_pair_colors <- function(pairs) {
  set.seed(123L)
  n <- length(pairs)
  cols <- rgb(runif(n), runif(n), runif(n))
  setNames(cols, pairs)
}


# =============================================================================
# MAIN
# =============================================================================

message("--- 05_outliers: cross-mouse outlier analysis ---")

all_samples_bc <- c(cohort_colonized, cohort_colonized_2)

mouse_results <- list()
for (nm in all_samples_bc) {
  mouse_results[[nm]] <- tryCatch(
    run_outliers(nm),
    error = function(e) {
      message("  ERROR in ", nm, ": ", conditionMessage(e))
      list(out_raw = NULL, out_delta = NULL)
    }
  )
}

# ---- Pool all outlier pairs across mice ----
message("  Aggregating cross-mouse pairs ...")

raw_pairs_df <- bind_rows(lapply(names(mouse_results), function(nm) {
  d <- mouse_results[[nm]]$out_raw
  if (is.null(d) || nrow(d) == 0) return(NULL)
  tibble(pair_directional = unique(d$pair_directional), mouse = nm, type = "raw")
}))

delta_pairs_df <- bind_rows(lapply(names(mouse_results), function(nm) {
  d <- mouse_results[[nm]]$out_delta
  if (is.null(d) || nrow(d) == 0) return(NULL)
  tibble(pair_directional = unique(d$pair_directional), mouse = nm, type = "delta")
}))

combined_df <- bind_rows(raw_pairs_df, delta_pairs_df)

if (nrow(combined_df) > 0) {
  pair_summary <- combined_df %>%
    distinct() %>%
    count(pair_directional, type, name = "num_models")

  pair_wide <- pair_summary %>%
    pivot_wider(names_from = type, values_from = num_models,
                values_fill = 0L)

  # Common pairs: in both raw AND delta, in >= dcm_min_mice mice
  common_pairs_tbl <- pair_wide %>%
    filter(raw >= dcm_min_mice & delta >= dcm_min_mice)

  out_tbl <- here::here("results/tables/DCM")
  write_csv(common_pairs_tbl,
            file.path(out_tbl, "common_outlier_pairs.csv"))
  message("  Common pairs (>= ", dcm_min_mice, " mice): ",
          nrow(common_pairs_tbl))

  # ---- Bar plot: cross-mouse pair counts ----
  plot_df <- pair_summary %>%
    filter(pair_directional %in% common_pairs_tbl$pair_directional) %>%
    filter(num_models >= dcm_min_mice)

  if (nrow(plot_df) > 0) {
    p_bar <- ggplot(plot_df,
                    aes(x = reorder(pair_directional, num_models),
                        y = num_models, fill = type)) +
      geom_col(position = position_dodge()) +
      coord_flip() +
      scale_fill_manual(values = c("raw" = "#2b8cbe", "delta" = "#e34a33"),
                        name = "Outlier type") +
      labs(
        title = paste0("Common outlier pairs \u2014 \u2265 ", dcm_min_mice,
                        " mice (raw + delta)"),
        x = "Directional pair",
        y = "Number of mice flagged"
      ) +
      theme_minimal(base_size = 10)

    save_fig_dcm(p_bar, "common_outlier_pairs", nm = NULL,
                 w = 14, h = max(6, 0.3 * nrow(plot_df)))
  }

  # ---- Per-mouse tile + ridge plots ----
  all_pairs     <- unique(combined_df$pair_directional)
  pair_colors   <- generate_pair_colors(all_pairs)
  common_used   <- unique(common_pairs_tbl$pair_directional)

  for (nm in all_samples_bc) {
    res <- mouse_results[[nm]]
    if (is.null(res$df_full)) next

    tryCatch({
      p <- plot_outliers_mouse(
        df_full        = res$df_full,
        out_raw        = res$out_raw,
        common_pairs_used = common_used,
        pair_colors    = pair_colors,
        phase_lines    = res$phase_lines,
        nm             = nm
      )
      save_fig_dcm(p, paste0(nm, "_outliers"), nm = nm, w = 16, h = 7)
    }, error = function(e)
      message("  Plot ERROR in ", nm, ": ", conditionMessage(e)))
  }

} else {
  message("  No outlier pairs found across mice.")
}

message("--- 05_outliers: done ---")
