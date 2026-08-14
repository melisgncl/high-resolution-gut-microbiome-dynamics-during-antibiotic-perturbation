# =============================================================================
# Title:   02_process.R — Rarefy, compute relative abundance, interpolate
# Author:  Melis Gencel
# Date:    2026-03-19
# Input:   samples_raw      — named list from 01_import.R
#          metadata/0_config.R  (rarefaction_thresholds, interpolation_targets)
# Output:  samples_family   — named list of 12 tibbles (in memory)
#          results/tables/16S/family/<sample>_family.csv  (intermediate saves)
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(vegan)
library(zoo)

source(here::here("metadata/0_config.R"))
source(here::here("analysis/01_16S/01_import.R"))

# =============================================================================
# SANITY CHECK — inspect raw input structure before processing
# =============================================================================

message("--- 02_process: input inspection ---")
for (nm in names(samples_raw)[1:3]) {
  message(nm, ": ", nrow(samples_raw[[nm]]), " rows")
  message("  columns: ", paste(colnames(samples_raw[[nm]]), collapse = ", "))
  print(head(samples_raw[[nm]], 3))
}

# =============================================================================
# RAREFY
# =============================================================================

# Rarefy each timepoint independently to the cohort threshold.
# OTUs with total reads <= min_abundance across ALL timepoints are dropped first
# (consistent with archive: removes noise, not biologically meaningful OTUs).
# Returns tibble: Time, Family, Abundance, Abundance.family
rarefy_sample <- function(df, threshold, min_abundance = 10) {
  # Filter low-abundance OTUs
  df_clean <- df %>%
    dplyr::select(OTU, Abundance, Time, Family) %>%
    group_by(OTU) %>%
    mutate(total = sum(Abundance)) %>%
    filter(total > min_abundance) %>%
    ungroup() %>%
    dplyr::select(-total)

  # Rarefy each timepoint independently
  rarefied <- lapply(unique(df_clean$Time), function(tp) {
    sub        <- df_clean %>% filter(Time == tp)
    total_seqs <- sum(sub$Abundance)

    if (total_seqs >= threshold) {
      wide <- sub %>%
        dplyr::select(OTU, Abundance) %>%
        pivot_wider(names_from = "OTU", values_from = "Abundance",
                    values_fill = 0L) %>%
        as.data.frame()
      rare <- vegan::rrarefy(wide, sample = threshold)
      as_tibble(rare) %>%
        pivot_longer(cols = everything(),
                     names_to = "OTU", values_to = "Abundance") %>%
        mutate(Time = tp)
    } else {
      message("  WARNING: Time ", tp, " has ", total_seqs,
              " reads (< threshold ", threshold, ") — kept un-rarefied")
      sub %>% dplyr::select(OTU, Abundance, Time)
    }
  })

  rarefied_data <- bind_rows(rarefied)

  # Re-attach Family annotations then summarise at family level
  family_lookup <- df_clean %>%
    dplyr::select(OTU, Family) %>%
    distinct()

  rarefied_data %>%
    left_join(family_lookup, by = "OTU") %>%
    group_by(Time, Family) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    group_by(Time) %>%
    mutate(Abundance.family = Abundance / sum(Abundance)) %>%
    ungroup()
}


# =============================================================================
# INTERPOLATE
# =============================================================================

# Linearly interpolate missing timepoints for target families via zoo::na.approx.
# Non-target families: their zeros/NAs are left as-is (zeroed after the step).
# Families where ALL Abundance.family values are zero/NA are removed (absent taxa).
#
# NOTE: This function does NOT apply an abundance threshold filter.
#       Threshold filtering is deferred to 03_filter.R so that all families
#       are available for cohort-level decisions.
interpolate_sample <- function(df, targets) {
  no_targets <- length(targets) == 0 || all(targets == "")

  df_out <- df %>%
    # Zeros → NA so interpolation can bridge true missing timepoints
    mutate(Abundance.family = na_if(Abundance.family, 0)) %>%
    # Remove families that are entirely absent (all NA across all times)
    group_by(Family) %>%
    filter(!all(is.na(Abundance.family))) %>%
    ungroup()

  if (!no_targets) {
    df_out <- df_out %>%
      arrange(Time) %>%
      group_by(Family) %>%
      group_modify(~ {
        if (.y$Family %in% targets &&
            any(is.na(.x$Abundance.family)) &&
            sum(!is.na(.x$Abundance.family)) >= 2) {
          .x$Abundance.family <- zoo::na.approx(
            .x$Abundance.family, x = .x$Time, na.rm = FALSE
          )
        }
        .x
      }) %>%
      ungroup()
  }

  df_out %>%
    mutate(Abundance.family = if_else(is.na(Abundance.family), 0, Abundance.family))
}


# =============================================================================
# PROCESS ALL 12 SAMPLES
# =============================================================================

process_all_samples <- function(samples_raw) {
  message("--- 02_process: rarefaction + interpolation ---")

  result <- lapply(names(samples_raw), function(nm) {
    message("  Processing ", nm, "...")
    df  <- samples_raw[[nm]]
    thr <- rarefaction_thresholds[nm]
    tgt <- interpolation_targets[[nm]]

    family_df         <- rarefy_sample(df, threshold = thr)
    family_df         <- interpolate_sample(family_df, targets = tgt)
    family_df$Sample  <- nm
    family_df
  })
  names(result) <- names(samples_raw)
  result
}


# =============================================================================
# RUN
# =============================================================================

samples_family <- process_all_samples(samples_raw)

# Save intermediate CSVs — downstream scripts can load these directly
# without re-running the expensive rarefaction step.
out_dir <- file.path(path_tables, "16S", "family")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

for (nm in names(samples_family)) {
  write_csv(samples_family[[nm]],
            file.path(out_dir, paste0(nm, "_family.csv")))
}
message("Saved family CSVs to ", out_dir)

# Sanity check
for (nm in names(samples_family)) {
  df <- samples_family[[nm]]
  message(nm, ": ", nrow(df), " rows | times: ",
          paste(sort(unique(df$Time)), collapse = " "),
          " | families: ", length(unique(df$Family)))
}
