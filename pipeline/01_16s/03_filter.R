# =============================================================================
# Title:   03_filter.R — Filter families per cohort, bin "other", write TSVs
# Author:  Melis Gencel
# Date:    2026-03-19
# Input:   samples_family   — named list from 02_process.R
#          metadata/0_config.R
# Output:  samples_filtered — named list of 12 tibbles (in memory)
#          results/tables/16S/coclustering/<sample>_16S_taxa.tsv
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(zoo)

source(here::here("metadata/0_config.R"))
source(here::here("analysis/01_16S/02_process.R"))

# =============================================================================
# COHORT-LEVEL FAMILY FILTER
# =============================================================================

# Families with mean Abundance.family > 1e-5 across all timepoints and all
# samples in the cohort are retained as "named" families.
# Threshold matches archive (1_family_analysis_new.R lines 121, 127, 134).
get_main_families <- function(cohort_samples, threshold = 1e-5) {
  bind_rows(cohort_samples) %>%
    group_by(Family) %>%
    summarise(mean_abund = mean(Abundance.family, na.rm = TRUE), .groups = "drop") %>%
    filter(mean_abund > threshold) %>%
    pull(Family)
}

m_main   <- get_main_families(samples_family[cohort_colonized])
c2_main  <- get_main_families(samples_family[cohort_colonized_2])
c_main   <- get_main_families(samples_family[cohort_controls])

message("Main families — colonized_1: ", length(m_main),
        " | colonized_2: ", length(c2_main),
        " | controls: ", length(c_main))


# =============================================================================
# BIN RARE FAMILIES AS "other"
# =============================================================================

bin_to_other <- function(df, main_families) {
  df %>%
    mutate(Family = if_else(Family %in% main_families, Family, "other")) %>%
    group_by(Time, Family, Sample) %>%
    summarise(
      Abundance        = sum(Abundance),
      Abundance.family = sum(Abundance.family),
      .groups = "drop"
    )
}


# =============================================================================
# APPLY FILTER PER COHORT
# =============================================================================

filter_cohort <- function(sample_names, main_families) {
  lapply(sample_names, function(nm) {
    bin_to_other(samples_family[[nm]], main_families)
  }) %>% setNames(sample_names)
}

filtered_colonized   <- filter_cohort(cohort_colonized,   m_main)
filtered_colonized_2 <- filter_cohort(cohort_colonized_2, c2_main)
filtered_controls    <- filter_cohort(cohort_controls,    c_main)

# Drop bT0 (time index 0) from colonized_2 — pre-colonisation baseline
# (colonized_1 cohort time 0 was already dropped in 01_import.R)
filtered_colonized_2 <- lapply(filtered_colonized_2, function(df) {
  filter(df, Time != 0L)
})

samples_filtered <- c(filtered_colonized, filtered_colonized_2, filtered_controls)

# Sanity check
for (nm in names(samples_filtered)) {
  df <- samples_filtered[[nm]]
  message(nm, ": ", length(unique(df$Family)), " families | ",
          length(unique(df$Time)), " timepoints")
}


# =============================================================================
# CO-CLUSTERING TSV OUTPUT
# =============================================================================

# Pivot to Time × Family wide format. Missing timepoints within the observed
# range are filled by linear interpolation (zoo::na.approx) — replicates the
# per-mouse row-index interpolation in the archive's 1_family_analysis_new.R
# but generalised to handle any missing-timepoint pattern.

out_dir_coclust <- file.path(path_tables, "16S", "coclustering")
if (!dir.exists(out_dir_coclust)) dir.create(out_dir_coclust, recursive = TRUE)

write_coclustering_tsv <- function(df, sample_id, main_families, full_time_range) {
  # Use main families only (exclude the "other" bin)
  df_main <- df %>%
    filter(Family %in% main_families) %>%
    mutate(Time = as.integer(Time))

  wide <- df_main %>%
    dplyr::select(Time, Family, Abundance.family) %>%
    pivot_wider(names_from = Family, values_from = Abundance.family)

  # Ensure complete time range; interpolate gaps linearly
  full_times    <- tibble(Time = as.integer(full_time_range))
  wide_complete <- full_times %>%
    left_join(wide, by = "Time") %>%
    arrange(Time) %>%
    mutate(across(-Time, ~ zoo::na.approx(.x, x = Time, na.rm = FALSE))) %>%
    mutate(across(-Time, ~ if_else(is.na(.x), 0, .x)))

  out_path <- file.path(out_dir_coclust, paste0(sample_id, "_16S_taxa.tsv"))
  write_tsv(wide_complete, out_path)
  message("  Wrote: ", out_path)
  invisible(wide_complete)
}

message("--- 03_filter: writing co-clustering TSVs ---")

# Colonized (m1-m4): times 1-19
for (nm in cohort_colonized) {
  write_coclustering_tsv(samples_family[[nm]], nm, m_main, 1:19)
}

# Colonized 2 (m5-m8): times 1-18
for (nm in cohort_colonized_2) {
  write_coclustering_tsv(samples_family[[nm]], nm, c2_main, 1:18)
}

# Controls (c_m1-c_m4): use observed time range (varies per mouse)
for (nm in cohort_controls) {
  t_range <- seq(
    min(samples_family[[nm]]$Time),
    max(samples_family[[nm]]$Time)
  )
  write_coclustering_tsv(samples_family[[nm]], nm, c_main, t_range)
}

message("Co-clustering TSVs written to ", out_dir_coclust)
