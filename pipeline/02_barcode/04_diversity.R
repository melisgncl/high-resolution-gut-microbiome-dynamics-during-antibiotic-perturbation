# =============================================================================
# Title:   04_diversity.R — Hill diversity indices on barcode frequencies
# Author:  Melis Gencel
# Date:    2026-03-19
# Input:   samples_bc — named list from 02_process.R
#          metadata/0_config.R
# Output:  results/tables/barcodes/diversity_bc_all.csv
#          results/figures/barcodes/diversity_bc_colonized.pdf/.png
#          results/figures/barcodes/diversity_bc_colonized_2.pdf/.png
#
# Hill number definitions (Chao et al. 2014):
#   q=0   : lineage richness (# barcodes with Freq > 0)
#   q=1   : exp(Shannon entropy) — effective # equally frequent lineages
#   q=inf : 1 / max(Freq)       — inverse dominance (Berger-Parker)
#
# Computed on barcode frequencies (not 16S). Complements 01_16S/05_diversity.R.
# =============================================================================

library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)

source(here::here("metadata/0_config.R"))
source(here::here("analysis/02_barcode/02_process.R"))


# =============================================================================
# DIVERSITY CALCULATION — one sample
# =============================================================================

calc_diversity_bc <- function(df) {
  df %>%
    group_by(Time, Sample) %>%
    summarise(
      q_0 = sum(Freq > 0),
      q_1 = {
        p <- Freq[Freq > 0]
        p <- p / sum(p)
        exp(-sum(p * log(p)))
      },
      q_inf = {
        if (all(Freq == 0)) NA_real_
        else 1 / max(Freq)
      },
      .groups = "drop"
    ) %>%
    mutate(Time = as.integer(as.character(Time)))
}


# =============================================================================
# COMPUTE FOR ALL SAMPLES
# =============================================================================

message("--- 04_diversity (barcode): computing Hill indices ---")

diversity_bc <- bind_rows(
  lapply(names(samples_bc), function(nm) calc_diversity_bc(samples_bc[[nm]]))
) %>%
  mutate(
    Cohort = case_when(
      Sample %in% cohort_colonized   ~ "colonized",
      Sample %in% cohort_colonized_2 ~ "colonized_2"
    ),
    Cohort = factor(Cohort, levels = c("colonized", "colonized_2"))
  )

str(diversity_bc)
head(diversity_bc)


# =============================================================================
# SAVE TABLE
# =============================================================================

out_dir_tables <- file.path(path_tables, "barcodes")
if (!dir.exists(out_dir_tables)) dir.create(out_dir_tables, recursive = TRUE)

write_csv(diversity_bc, file.path(out_dir_tables, "diversity_bc_all.csv"))
message("Saved: ", file.path(out_dir_tables, "diversity_bc_all.csv"))


# =============================================================================
# PLOT HELPER — three-panel diversity figure for one cohort
# =============================================================================

make_diversity_panels_bc <- function(cohort_samples, cohort_key, cohort_pal) {
  sub <- diversity_bc %>%
    filter(Sample %in% cohort_samples) %>%
    mutate(Sample = factor(Sample, levels = cohort_samples))

  pal <- setNames(cohort_pal, cohort_samples)
  br  <- breaks.hash[[cohort_key]]
  lb  <- labels.hash[[cohort_key]]
  lm  <- limits.hash[[cohort_key]]

  base_panel <- function(y_var, y_label) {
    ggplot(sub, aes(x = Time, y = .data[[y_var]],
                    group = Sample, colour = Sample)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2) +
      scale_colour_manual(values = pal) +
      scale_x_continuous(limits = lm, breaks = br, labels = lb) +
      labs(y = y_label, x = "Time (post-gavage)") +
      theme_Publication(base_size = 14, aspect.ratio = 0.75) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  p0   <- base_panel("q_0",   "Richness (q = 0)") +
    theme(legend.position = "none")
  p1   <- base_panel("q_1",   "Eff. lineages (q = 1)") +
    theme(legend.position = "none")
  pinf <- base_panel("q_inf", "Inv. dominance (q = \u221e)") +
    theme(legend.position = "right")

  p0 | p1 | pinf
}


# =============================================================================
# GENERATE FIGURES
# =============================================================================

message("--- 04_diversity (barcode): generating figures ---")

save_fig_bc(
  make_diversity_panels_bc(cohort_colonized, "m",  pal1),
  "diversity_bc_colonized"
)
save_fig_bc(
  make_diversity_panels_bc(cohort_colonized_2, "c2", pal2),
  "diversity_bc_colonized_2"
)

message("--- 04_diversity (barcode): done ---")
