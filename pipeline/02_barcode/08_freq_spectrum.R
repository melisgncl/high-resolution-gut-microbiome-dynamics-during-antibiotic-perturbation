# =============================================================================
# Title:   08_freq_spectrum.R — Barcode frequency spectrum & rank-abundance
# Author:  Melis Gencel
# Date:    2026-03-19
# Input:   samples_bc — named list from 02_process.R
#          metadata/0_config.R
# Output:  results/figures/barcodes/freq_spectrum_<cohort>.pdf/.png
#          results/figures/barcodes/rank_abundance.pdf/.png
#          results/tables/barcodes/freq_spectrum_summary.csv
#
# Analyses:
#   1. Frequency spectrum — histogram of log10(Freq) at each timepoint.
#      Under neutral drift, expect a 1/f-like distribution (many rare, few
#      abundant). Deviations (bumps at high freq) indicate selection.
#
#   2. Rank-abundance (Whittaker) plot — barcodes ranked by frequency at
#      each timepoint. Tracks how the abundance hierarchy shifts over time.
#      Steep early rank curves → high dominance; flat → even community.
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readr)
library(scales)

source(here::here("metadata/0_config.R"))
source(here::here("analysis/02_barcode/02_process.R"))


# =============================================================================
# SECTION 1 — FREQUENCY SPECTRUM PER COHORT
# =============================================================================

message("--- 08_freq_spectrum: frequency spectrum analysis ---")

# Focus on non-zero frequencies to avoid log(0)
freq_long <- bind_rows(
  lapply(names(samples_bc), function(nm) {
    samples_bc[[nm]] %>%
      filter(Freq > 0) %>%
      mutate(
        Time_int = as.integer(as.character(Time)),
        Cohort   = if (nm %in% cohort_colonized) "colonized" else "colonized_2"
      ) %>%
      select(Sample, Cohort, Time_int, ID, Freq)
  })
) %>%
  mutate(Cohort = factor(Cohort, levels = c("colonized", "colonized_2")))

# Representative timepoints to show (avoid overplotting)
key_times_m  <- c(1, 3, 5, 9, 13, 17, 19)
key_times_nm <- c(1, 3, 5, 9, 13, 17, 18)

plot_spectrum <- function(cohort_label, key_times, cohort_pal_vec, nm_vec) {
  sub <- freq_long %>%
    filter(Cohort == cohort_label,
           Time_int %in% key_times,
           Sample %in% nm_vec) %>%
    mutate(
      time_label = paste0("t=", Time_int),
      time_label = factor(time_label,
                           levels = paste0("t=", sort(unique(Time_int))))
    )

  pal_times <- scales::hue_pal()(length(key_times))
  names(pal_times) <- paste0("t=", sort(key_times))

  ggplot(sub, aes(x = log10(Freq), colour = time_label)) +
    geom_density(linewidth = 0.8, fill = NA) +
    geom_vline(xintercept = log10(1e-3), linetype = "dotted",
               colour = "gray60", linewidth = 0.5) +
    scale_colour_manual(values = pal_times, name = "Timepoint") +
    scale_x_continuous(
      breaks = c(-7, -5, -3, -1, 0),
      labels = c("10\u207b\u2077", "10\u207b\u2075", "10\u207b\u00b3",
                 "10\u207b\u00b9", "1")
    ) +
    facet_wrap(~ Sample, nrow = 1) +
    labs(
      title    = paste0("Frequency spectrum — ", cohort_label),
      subtitle = "Density of log\u2081\u2080(barcode frequency) at selected timepoints",
      x        = "Barcode frequency (log\u2081\u2080)",
      y        = "Density"
    ) +
    theme_Publication(base_size = 12, aspect.ratio = 0.9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right")
}

p_spec_col <- plot_spectrum("colonized", key_times_m,  pal1, cohort_colonized)
p_spec_c2  <- plot_spectrum("colonized_2", key_times_nm, pal2, cohort_colonized_2)

save_fig_bc(p_spec_col, "freq_spectrum_colonized", w = 25, h = 8)
save_fig_bc(p_spec_c2,  "freq_spectrum_colonized_2", w = 25, h = 8)


# =============================================================================
# SECTION 2 — RANK-ABUNDANCE (WHITTAKER) PLOT
# =============================================================================

message("--- 08_freq_spectrum: rank-abundance curves ---")

# For a set of representative timepoints, rank barcodes by frequency
plot_rank_abundance <- function(cohort_label, key_times, nm_vec) {
  sub <- freq_long %>%
    filter(Cohort == cohort_label,
           Time_int %in% key_times,
           Sample %in% nm_vec,
           Freq > 0) %>%
    group_by(Sample, Time_int) %>%
    mutate(rank = rank(-Freq, ties.method = "first")) %>%
    ungroup() %>%
    mutate(
      time_label = paste0("t=", Time_int),
      time_label = factor(time_label, levels = paste0("t=", sort(key_times)))
    )

  pal_times <- scales::hue_pal()(length(key_times))
  names(pal_times) <- paste0("t=", sort(key_times))

  ggplot(sub, aes(x = rank, y = log10(Freq),
                  colour = time_label, group = time_label)) +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    scale_colour_manual(values = pal_times, name = "Timepoint") +
    facet_wrap(~ Sample, nrow = 1) +
    labs(
      title = paste0("Rank-abundance curves — ", cohort_label),
      x     = "Rank (most → least abundant)",
      y     = "Frequency (log\u2081\u2080)"
    ) +
    theme_Publication(base_size = 12, aspect.ratio = 0.9) +
    theme(legend.position = "right")
}

p_ra_col <- plot_rank_abundance("colonized", key_times_m,  cohort_colonized)
p_ra_c2  <- plot_rank_abundance("colonized_2", key_times_nm, cohort_colonized_2)

save_fig_bc(p_ra_col | p_ra_c2, "rank_abundance", w = 30, h = 8)


# =============================================================================
# SECTION 3 — SUMMARY TABLE: mean, median, max frequency per timepoint
# =============================================================================

freq_summary <- freq_long %>%
  group_by(Sample, Cohort, Time_int) %>%
  summarise(
    n_detected = n(),
    mean_freq  = mean(Freq),
    median_freq= median(Freq),
    max_freq   = max(Freq),
    frac_above_1pct = mean(Freq > 0.01),
    .groups    = "drop"
  )

out_dir_tables <- file.path(path_tables, "barcodes")
if (!dir.exists(out_dir_tables)) dir.create(out_dir_tables, recursive = TRUE)

write_csv(freq_summary,
          file.path(out_dir_tables, "freq_spectrum_summary.csv"))
message("Saved: ", file.path(out_dir_tables, "freq_spectrum_summary.csv"))

message("--- 08_freq_spectrum (barcode): done ---")
