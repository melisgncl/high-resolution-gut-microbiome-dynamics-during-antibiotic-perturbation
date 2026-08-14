# =============================================================================
# Title:   04_plot_community.R — Community composition area and line plots
# Author:  Melis Gencel
# Date:    2026-03-19
# Input:   samples_filtered — named list from 03_filter.R
#          metadata/0_config.R (themes, family.colors, breaks.hash, etc.)
# Output:  results/figures/16S/
#            area_colonized.pdf/.png  — stacked area, 4-panel (m1-m4)
#            area_colonized_2.pdf/.png
#            area_controls.pdf/.png
#            log_colonized.pdf/.png   — log-scale line, 4-panel
#            log_colonized_2.pdf/.png
#            log_controls.pdf/.png
# =============================================================================

library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)

source(here::here("metadata/0_config.R"))
source(here::here("analysis/01_16S/03_filter.R"))

# =============================================================================
# OUTPUT DIRECTORY
# =============================================================================

if (!dir.exists(path_results_16S)) dir.create(path_results_16S, recursive = TRUE)


# save_fig() is defined in metadata/0_config.R (Section 12)

# =============================================================================
# STACKED AREA — relative abundance per family
# =============================================================================

plot_area_sample <- function(df, sample_id, breaks, labels, limits) {
  df_plot <- df %>%
    mutate(Family = factor(Family, levels = rev(major.family))) %>%
    arrange(Time, Family)

  ggplot(df_plot, aes(x = Time, y = Abundance.family, fill = Family)) +
    geom_area(position = "stack") +
    scale_fill_manual(values = family.colors, drop = FALSE) +
    scale_x_continuous(limits = limits, breaks = breaks, labels = labels) +
    scale_y_continuous(limits = c(0, 1),
                       labels = scales::percent_format(accuracy = 1)) +
    labs(title = sample_id, x = NULL, y = "Relative abundance") +
    theme_Publication(base_size = 14, aspect.ratio = 0.75) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1))
}


# =============================================================================
# LOG-SCALE LINE — family dynamics over time
# =============================================================================

plot_log_sample <- function(df, sample_id, breaks, labels, limits) {
  df_plot <- df %>%
    filter(Family != "other", Abundance.family > 0) %>%
    mutate(Family = factor(Family, levels = major.family))

  ggplot(df_plot, aes(x = Time, y = Abundance.family,
                      group = Family, color = Family)) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = family.colors, drop = FALSE) +
    scale_x_continuous(limits = limits, breaks = breaks, labels = labels) +
    scale_y_log10(
      limits = c(1e-6, 1),
      breaks = c(1e-6, 1e-4, 1e-3, 1e-2, 1e-1, 1),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    labs(title = sample_id, x = NULL, y = "Relative abundance (log\u2081\u2080)") +
    theme_Publication(base_size = 14, aspect.ratio = 0.75) +
    guides(color = guide_legend(override.aes = list(linewidth = 3))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# =============================================================================
# FOUR-PANEL FIGURE PER COHORT
# =============================================================================

make_cohort_panels <- function(sample_names, cohort_key, plot_fn) {
  br  <- breaks.hash[[cohort_key]]
  lb  <- labels.hash[[cohort_key]]
  lm  <- limits.hash[[cohort_key]]

  panels <- lapply(sample_names, function(nm) {
    plot_fn(df        = samples_filtered[[nm]],
            sample_id = nm,
            breaks    = br,
            labels    = lb,
            limits    = lm)
  })

  wrap_plots(panels, nrow = 1) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}


# =============================================================================
# GENERATE ALL FIGURES
# =============================================================================

message("--- 04_plot_community: generating figures ---")

# Colonized (m1-m4)
save_fig(make_cohort_panels(cohort_colonized, "m",  plot_area_sample), "area_colonized")
save_fig(make_cohort_panels(cohort_colonized, "m",  plot_log_sample),  "log_colonized")

# Colonized 2 (m5-m8)
save_fig(make_cohort_panels(cohort_colonized_2, "c2", plot_area_sample), "area_colonized_2")
save_fig(make_cohort_panels(cohort_colonized_2, "c2", plot_log_sample),  "log_colonized_2")

# Controls (c_m1-c_m4)
save_fig(make_cohort_panels(cohort_controls,  "c",  plot_area_sample), "area_controls")
save_fig(make_cohort_panels(cohort_controls,  "c",  plot_log_sample),  "log_controls")

message("--- 04_plot_community: done ---")
