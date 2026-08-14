# =============================================================================
# Title:   03_dynamics.R — Barcode lineage dynamics plots
# Author:  Melis Gencel
# Date:    2026-03-19
# Input:   samples_bc — named list from 02_process.R
#          metadata/0_config.R  (breaks.hash, labels.hash, limits.hash, palettes)
# Output:  results/figures/barcodes/
#            dynamics_line_<sample>.pdf/.png  — log-scale line, one per sample
#            dynamics_area_<sample>.pdf/.png  — stacked area, one per sample
#            dynamics_line_combined.pdf/.png  — 2x4 panel, all 8 samples
#            dynamics_area_combined.pdf/.png  — 2x4 panel, all 8 samples
# Notes:   Line plot: all barcodes drawn in gray; top barcodes (hex_line !=
#          "#cccccc") overlaid in their paper-specific colors.
#          Area plot: all barcodes colored via hex_area (random palette for
#          rare, specific color for top). No ggnewscale needed.
# =============================================================================

library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

source(here::here("metadata/0_config.R"))
source(here::here("analysis/02_barcode/02_process.R"))

if (!dir.exists(path_results_bc)) dir.create(path_results_bc, recursive = TRUE)


# =============================================================================
# HELPER — COHORT KEY FOR AXIS SETTINGS
# =============================================================================

cohort_key <- function(sample_id) {
  if (sample_id %in% cohort_colonized) "m" else "c2"
}


# =============================================================================
# LOG-SCALE LINE PLOT (one sample)
# Gray background for all barcodes, top barcodes in specific colors.
# =============================================================================

plot_line_sample <- function(df, sample_id) {
  key <- cohort_key(sample_id)
  br  <- breaks.hash[[key]]
  lb  <- labels.hash[[key]]
  lm  <- limits.hash[[key]]

  df_plot <- df %>%
    filter(Freq > 0) %>%
    mutate(ID = factor(ID, levels = levels(ID)))

  # Split into background (gray) and foreground (colored)
  df_gray    <- df_plot %>% filter(hex_line == "#cccccc")
  df_colored <- df_plot %>% filter(hex_line != "#cccccc")

  # Build named color vector for colored barcodes
  color_vals <- df_colored %>%
    distinct(ID, hex_line) %>%
    { setNames(.$hex_line, as.character(.$ID)) }

  p <- ggplot() +
    # Background: all barcodes, thin gray lines
    geom_line(data = df_gray,
              aes(x = as.integer(as.character(Time)),
                  y = Freq, group = ID),
              colour = "#cccccc", linewidth = 0.3, alpha = 0.6) +
    # Foreground: top barcodes, thick colored lines
    geom_line(data = df_colored,
              aes(x = as.integer(as.character(Time)),
                  y = Freq, group = ID, colour = ID),
              linewidth = 1.2) +
    scale_colour_manual(values = color_vals, guide = "none") +
    scale_y_log10(
      limits = c(1e-7, 1),
      breaks = c(1e-7, 1e-5, 1e-3, 1e-1),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_x_continuous(breaks = br, labels = lb, limits = lm) +
    labs(title = sample_id,
         x = "Time (post-gavage)",
         y = "Barcode frequency (log\u2081\u2080)") +
    theme_Publication(base_size = 14, aspect.ratio = 0.75) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p
}


# =============================================================================
# STACKED AREA PLOT (one sample)
# All barcodes colored via hex_area — no ggnewscale needed.
# =============================================================================

plot_area_sample <- function(df, sample_id) {
  key <- cohort_key(sample_id)
  br  <- breaks.hash[[key]]
  lb  <- labels.hash[[key]]
  lm  <- limits.hash[[key]]

  # Build named fill vector: ID (as character) → hex_area
  color_vals <- df %>%
    distinct(ID, hex_area) %>%
    { setNames(.$hex_area, as.character(.$ID)) }

  df_plot <- df %>%
    mutate(Time_int = as.integer(as.character(Time)))

  p <- ggplot(df_plot,
              aes(x = Time_int, y = Freq,
                  group = ID, fill = as.character(ID))) +
    geom_area(position = "stack") +
    scale_fill_manual(values = color_vals, guide = "none") +
    scale_x_continuous(breaks = br, labels = lb, limits = lm) +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format(accuracy = 1)) +
    labs(title = sample_id,
         x = NULL,
         y = "Barcode frequency") +
    theme_Publication(base_size = 14, aspect.ratio = 0.75) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p
}


# =============================================================================
# GENERATE PER-SAMPLE FIGURES
# =============================================================================

message("--- 03_dynamics: generating per-sample figures ---")

for (nm in names(samples_bc)) {
  df <- samples_bc[[nm]]

  p_line <- plot_line_sample(df, nm)
  p_area <- plot_area_sample(df, nm)

  save_fig_bc(p_line, paste0("dynamics_line_", nm), w = 8.25, h = 6)
  save_fig_bc(p_area, paste0("dynamics_area_", nm), w = 8.25, h = 6)
}


# =============================================================================
# COMBINED 2x4 PANELS — line and area
# =============================================================================

# message("--- 03_dynamics: generating combined panels ---")
#
# make_combined <- function(plot_fn, title_label) {
#   # Row 1: colonized (m1-m4)
#   row1 <- lapply(cohort_colonized, function(nm) {
#     plot_fn(samples_bc[[nm]], nm) +
#       theme(plot.title = element_text(size = 10))
#   })
#   # Row 2: colonized_2 (m5-m8)
#   row2 <- lapply(cohort_colonized_2, function(nm) {
#     plot_fn(samples_bc[[nm]], nm) +
#       theme(plot.title = element_text(size = 10))
#   })
#
#   wrap_plots(c(row1, row2), nrow = 2)
# }
#
# save_fig_bc(make_combined(plot_line_sample, "line"),
#             "dynamics_line_combined", w = 33, h = 12)
# save_fig_bc(make_combined(plot_area_sample, "area"),
#             "dynamics_area_combined", w = 33, h = 12)

message("--- 03_dynamics (barcode): done ---")
