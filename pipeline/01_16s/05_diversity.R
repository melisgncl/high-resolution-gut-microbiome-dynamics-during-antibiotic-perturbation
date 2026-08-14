# =============================================================================
# Title:   05_diversity.R — Hill diversity indices (q=0, 1, ∞) for all cohorts
# Author:  Melis Gencel
# Date:    2026-03-19
# Input:   samples_filtered — named list from 03_filter.R
#          metadata/0_config.R
# Output:  results/figures/16S/diversity_<cohort>.pdf/.png
#          results/tables/16S/diversity_all.csv
#
# Hill number definitions (Chao et al. 2014):
#   q=0   : species richness (count of families with rel_abund > 0)
#   q=1   : exp(Shannon entropy) — effective number of equally common families
#   q=∞   : 1 / max(rel_abund)  — inverse Berger-Parker (dominance-weighted)
#
# "other" bin is excluded from all indices: it is not a biological entity.
# Relative abundances are renormalised after exclusion so they sum to 1.
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

source(here::here("metadata/0_config.R"))
source(here::here("analysis/01_16S/03_filter.R"))


# =============================================================================
# DIVERSITY CALCULATION
# =============================================================================

# Compute q_0, q_1, q_inf at every timepoint for one sample.
# Input df: long tibble with Family, Abundance.family, Time, Sample columns.
calc_diversity_sample <- function(df) {
  df %>%
    group_by(Time, Sample) %>%
    summarise(
      q_0 = {
        sum(Family != "other" & Abundance.family > 0)
      },
      q_1 = {
        p <- Abundance.family[Family != "other" & Abundance.family > 0]
        p <- p / sum(p)                          # renormalise to sum = 1
        exp(-sum(p * log(p)))
      },
      q_inf = {
        p_named <- Abundance.family[Family != "other"]
        if (length(p_named) == 0 || all(p_named == 0)) NA_real_
        else 1 / max(p_named)
      },
      .groups = "drop"
    )
}


# =============================================================================
# COMPUTE FOR ALL SAMPLES
# =============================================================================

message("--- 05_diversity: computing Hill indices ---")

diversity_all <- bind_rows(
  lapply(names(samples_filtered), function(nm) {
    calc_diversity_sample(samples_filtered[[nm]])
  })
) %>%
  mutate(
    Cohort = case_when(
      Sample %in% cohort_colonized   ~ "colonized",
      Sample %in% cohort_colonized_2 ~ "colonized_2",
      Sample %in% cohort_controls    ~ "controls"
    ),
    Cohort = factor(Cohort, levels = c("colonized", "colonized_2", "controls"))
  )

str(diversity_all)
head(diversity_all)


# =============================================================================
# SAVE TABLE
# =============================================================================

out_dir_tables <- file.path(path_tables, "16S")
if (!dir.exists(out_dir_tables)) dir.create(out_dir_tables, recursive = TRUE)

write_csv(diversity_all, file.path(out_dir_tables, "diversity_all.csv"))
message("Saved: ", file.path(out_dir_tables, "diversity_all.csv"))


# =============================================================================
# PLOT HELPERS
# =============================================================================

# save_fig() is defined in metadata/0_config.R (Section 12)

if (!dir.exists(path_results_16S)) dir.create(path_results_16S, recursive = TRUE)

# Build a three-panel diversity figure (q0 | q1 | q_inf) for one cohort.
# Each panel shows one line per sample.
make_diversity_panels <- function(cohort_samples, cohort_key, cohort_pal) {
  sub <- diversity_all %>%
    filter(Sample %in% cohort_samples) %>%
    mutate(Sample = factor(Sample, levels = cohort_samples))

  pal <- setNames(cohort_pal, cohort_samples)
  br  <- breaks.hash[[cohort_key]]
  lb  <- labels.hash[[cohort_key]]
  lm  <- limits.hash[[cohort_key]]

  base_aes <- function(y_var, y_label) {
    ggplot(sub, aes(x = Time, y = .data[[y_var]],
                    group = Sample, color = Sample)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 2) +
      scale_color_manual(values = pal) +
      scale_x_continuous(limits = lm, breaks = br, labels = lb) +
      labs(y = y_label, x = "Time (post-gavage)") +
      theme_Publication(base_size = 14, aspect.ratio = 0.75) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }

  p0    <- base_aes("q_0",   "Richness (q = 0)") +
    theme(legend.position = "none")
  p1    <- base_aes("q_1",   "Eff. families (q = 1)") +
    theme(legend.position = "none")
  p_inf <- base_aes("q_inf", "Inv. dominance (q = \u221e)") +
    theme(legend.position = "right")

  p0 | p1 | p_inf
}


# =============================================================================
# GENERATE DIVERSITY FIGURES
# =============================================================================

message("--- 05_diversity: generating figures ---")

save_fig(make_diversity_panels(cohort_colonized,   "m",  pal1), "diversity_colonized")
save_fig(make_diversity_panels(cohort_colonized_2, "c2", pal2), "diversity_colonized_2")
save_fig(make_diversity_panels(cohort_controls,    "c",  pal3), "diversity_controls")

message("--- 05_diversity: done ---")
