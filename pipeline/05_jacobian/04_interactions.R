# =============================================================================
# Title:   04_interactions.R — Jacobian interaction classification + density plots
# Author:  Melis Gencel
# Date:    2026-03-20
# Input:   results/tables/DCM/<nm>/<nm>_jacobian_long.csv    (from 01)
#          results/tables/DCM/lag_times.csv                  (from 02)
#          metadata/0_config.R
# Output:  results/tables/DCM/<nm>/<nm>_interactions.csv
#          results/figures/DCM/<nm>/<nm>_density_hgt.png/pdf
#          results/figures/DCM/<nm>/<nm>_interaction_proportions.png/pdf
# Notes:   Exact archive methodology (Gencel 2020):
#          - Interaction types classified from sign pattern of J[i,j] and J[j,i]:
#              (+,+) mutualism | (-,-) competition | (+,-) or (-,+) antagonism
#              (+,0) or (0,+) commensalism | (-,0) or (0,-) amensalism | (0,0) none
#          - Density plots split by before/after HGT using per-mouse lag time.
#          - Exact archive color palette (dcm_density_palette_anchors in config).
#          - Four-panel layout: all / C*-C* / C*-focal / species-species.
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)

source(here::here("metadata/0_config.R"))


# =============================================================================
# INTERACTION CLASSIFICATION
# =============================================================================

# Helper: is label a barcode cluster (C<number>)?
is_cluster <- function(x) grepl("^C[0-9]+$", x)

# Classify a pair by sign pattern of (strength_ij, strength_ji).
classify_interaction <- function(sij, sji) {
  dplyr::case_when(
    sij >  0 & sji >  0 ~ "mutualism",
    sij <  0 & sji <  0 ~ "competition",
    sij >  0 & sji <  0 ~ "antagonism",
    sij <  0 & sji >  0 ~ "antagonism",
    sij >  0 & sji == 0 ~ "commensalism",
    sij == 0 & sji >  0 ~ "commensalism",
    sij <  0 & sji == 0 ~ "amensalism",
    sij == 0 & sji <  0 ~ "amensalism",
    sij == 0 & sji == 0 ~ "none",
    TRUE                  ~ "none"
  )
}

# Build full interaction table with types and pair groups.
# long_df: window, species_i, species_j, strength
classify_all_interactions <- function(long_df, nm) {
  # Create reverse-lookup for J[j,i] strength
  reverse_df <- long_df %>%
    rename(species_i = species_j, species_j = species_i,
           strength_ji = strength) %>%
    select(window, species_i, species_j, strength_ji)

  merged <- long_df %>%
    left_join(reverse_df, by = c("window", "species_i", "species_j")) %>%
    mutate(
      strength_ji      = replace_na(strength_ji, 0),
      pair_directional = paste0(species_i, "\u2192", species_j),
      interaction_type = classify_interaction(strength, strength_ji),
      group = dplyr::case_when(
        is_cluster(species_i) & is_cluster(species_j)         ~ "C-C",
        is_cluster(species_i) & species_j == dcm_hgt_target   ~ "C-focal",
        species_i == dcm_hgt_target & is_cluster(species_j)   ~ "C-focal",
        !is_cluster(species_i) & !is_cluster(species_j)        ~ "species-species",
        TRUE                                                    ~ "other"
      ),
      mouse = nm
    )
  merged
}


# =============================================================================
# DENSITY PLOTS  (before / after HGT)
# =============================================================================

# Exact archive palette: 12 anchors interpolated to n_timepoints colors.
make_time_palette <- function(n) {
  colorRampPalette(dcm_density_palette_anchors)(n)
}

# Four-panel density figure: all / C-C / C-focal / species-species
plot_density_hgt <- function(df, lag_est, nm) {
  n_times      <- length(unique(df$window))
  time_palette <- make_time_palette(n_times)

  time_cutoff  <- if (!is.na(lag_est)) lag_est else max(df$window)

  df <- df %>%
    mutate(
      Time_Group = factor(
        ifelse(window > time_cutoff, "After HGT", "Before HGT"),
        levels = c("Before HGT", "After HGT")
      ),
      window_fac = as.factor(window)
    )

  make_plot <- function(data, title_text) {
    if (nrow(data) == 0) return(ggplot() + theme_void() + ggtitle(title_text))
    ggplot(data, aes(x = strength, colour = window_fac, group = window_fac)) +
      geom_density(linewidth = 0.8) +
      facet_wrap(~Time_Group, nrow = 1, scales = "free_y") +
      scale_colour_manual(values = time_palette, guide = "none") +
      labs(
        title = title_text,
        x     = "Interaction strength",
        y     = "Density"
      ) +
      theme_minimal(base_size = 10) +
      theme(
        strip.text   = element_text(face = "bold"),
        plot.title   = element_text(face = "bold", size = 11)
      )
  }

  p1 <- make_plot(df,                               "All interactions")
  p2 <- make_plot(filter(df, group == "C-C"),        paste0("C\u2013C"))
  p3 <- make_plot(filter(df, group == "C-focal"),    paste0("C\u2013", dcm_hgt_target))
  p4 <- make_plot(filter(df, group == "species-species"), "Species\u2013species")

  (p1 / p2 / p3 / p4) +
    plot_annotation(
      title = paste0(nm, " \u2014 Jacobian interaction strength (before/after HGT)"),
      theme = theme(plot.background = element_rect(fill = "white", colour = NA))
    )
}

# Stacked-bar proportion plot per interaction type over time
plot_interaction_proportions <- function(df, nm) {
  prop_df <- df %>%
    group_by(window, interaction_type) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(window) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup() %>%
    complete(window, interaction_type, fill = list(count = 0, proportion = 0))

  ggplot(prop_df, aes(x = window, y = proportion, fill = interaction_type)) +
    geom_col(position = "stack") +
    scale_fill_brewer(palette = "Set2", name = "Type") +
    labs(
      title = paste0(nm, " \u2014 interaction type proportions over time"),
      x = "Time window", y = "Proportion"
    ) +
    theme_minimal(base_size = 10)
}


# =============================================================================
# PER-MOUSE RUNNER
# =============================================================================

run_interactions <- function(nm, lag_tbl) {
  message("  ", nm, ": loading Jacobian long CSV ...")

  long_path <- here::here("results/tables/DCM", nm,
                            paste0(nm, "_jacobian_long.csv"))
  if (!file.exists(long_path)) {
    message("    Long CSV not found — skipping.")
    return(invisible(NULL))
  }
  long_df <- read_csv(long_path, show_col_types = FALSE)

  lag_est <- lag_tbl$lag[lag_tbl$mouse == nm]
  lag_est <- if (length(lag_est) == 0) NA_integer_ else lag_est[1]

  message("  ", nm, ": classifying interactions (lag = ", lag_est, ") ...")
  interactions <- classify_all_interactions(long_df, nm)

  out_tbl <- here::here("results/tables/DCM", nm)
  write_csv(interactions,
            file.path(out_tbl, paste0(nm, "_interactions.csv")))
  message("    Saved: ", nm, "_interactions.csv")

  # ---- Density plot ----
  p_density <- plot_density_hgt(interactions, lag_est, nm)
  save_fig_dcm(p_density, paste0(nm, "_density_hgt"), nm = nm, w = 14, h = 16)

  # ---- Proportion plot ----
  p_prop <- plot_interaction_proportions(interactions, nm)
  save_fig_dcm(p_prop, paste0(nm, "_interaction_proportions"),
               nm = nm, w = 10, h = 5)

  message("  ", nm, ": done.")
  invisible(interactions)
}


# =============================================================================
# RUN ALL MICE
# =============================================================================

message("--- 04_interactions: interaction classification + density plots ---")

lag_tbl <- tryCatch(
  read_csv(here::here("results/tables/DCM/lag_times.csv"),
           show_col_types = FALSE),
  error = function(e) {
    warning("lag_times.csv not found — lag time will be NA for all mice.")
    data.frame(mouse = character(), lag = integer())
  }
)

all_samples_bc <- c(cohort_colonized, cohort_colonized_2)

for (nm in all_samples_bc) {
  tryCatch(
    run_interactions(nm, lag_tbl),
    error = function(e)
      message("  ERROR in ", nm, ": ", conditionMessage(e))
  )
}

# Controls (16S only — no barcode clusters, but classify anyway)
for (nm in cohort_controls) {
  tryCatch(
    run_interactions(nm, lag_tbl),
    error = function(e)
      message("  ERROR in ", nm, ": ", conditionMessage(e))
  )
}

message("--- 04_interactions: done ---")
