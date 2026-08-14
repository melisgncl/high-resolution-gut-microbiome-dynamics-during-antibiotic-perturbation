# =============================================================================
# Title:   00_filter.R — Barcode pre-filtering for clustering
# Author:  Melis Gencel
# Date:    2026-03-20
# Input:   samples_bc — named list from 02_process.R
#          metadata/0_config.R
# Output:  results/tables/clustering/<nm>/<nm>_filtered.csv
#            Wide matrix: rows = barcodes, cols = timepoints + mean + points
#            Only barcodes meeting frequency + persistence thresholds.
# Notes:   freq_threshold: minimum mean frequency across all timepoints
#          time_threshold: minimum number of timepoints with Freq > 0
#          Noise IDs and noise timepoints already removed in samples_bc.
#          Creates per-mouse subdirectories for tables and figures.
#
# TO ADD A NEW SAMPLE: edit metadata/0_config.R only — see the instructions
#          at the top of that file. Do not edit this script.
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(readr)

source(here::here("metadata/0_config.R"))
source(here::here("analysis/02_barcode/02_process.R"))


# =============================================================================
# PARAMETERS
# =============================================================================

freq_threshold <- 1e-5   # minimum mean frequency across all timepoints

# Per-sample minimum number of timepoints with Freq > 0.
# Default 11; override specific samples.
time_threshold_default <- 11L
time_threshold_override <- c(
  m1 = NA_integer_, m2 = NA_integer_, m3 = NA_integer_, m4 = NA_integer_,
  m5 = NA_integer_, m6 = 10L,         m7 = 10L,         m8 = NA_integer_
)

resolve_time_threshold <- function(nm) {
  ov <- time_threshold_override[[nm]]
  if (!is.na(ov)) ov else time_threshold_default
}


# =============================================================================
# FILTER & RESHAPE
# =============================================================================

message("--- 00_filter (clustering): filtering barcodes for clustering ---")
message("  freq_threshold    = ", freq_threshold)
message("  time_threshold    = ", time_threshold_default, " (default)")

filter_for_clustering <- function(nm) {
  tt <- resolve_time_threshold(nm)
  message("  ", nm, " (time_threshold = ", tt, "):")

  df <- samples_bc[[nm]] %>%
    mutate(Time_int = as.integer(as.character(Time)))

  # Wide: rows = barcodes, cols = timepoints (integer names), values = Freq
  wide <- df %>%
    select(ID, Time_int, Freq) %>%
    pivot_wider(names_from  = Time_int,
                values_from = Freq,
                values_fill = 0)

  time_cols <- setdiff(names(wide), "ID")

  wide <- wide %>%
    mutate(
      mean   = rowMeans(select(., all_of(time_cols)), na.rm = TRUE),
      points = rowSums(select(., all_of(time_cols)) > 0, na.rm = TRUE)
    )

  filtered <- wide %>%
    filter(points >= tt, mean >= freq_threshold)

  message("    ", nrow(wide), " barcodes → ", nrow(filtered),
          " pass freq + persistence filter")

  if (nrow(filtered) == 0)
    stop(nm, ": 0 barcodes pass thresholds — adjust freq_threshold or time_threshold.")

  filtered
}


# =============================================================================
# CREATE DIRECTORIES AND SAVE
# =============================================================================

all_samples_bc <- c(cohort_colonized, cohort_colonized_2)

filtered_bc <- setNames(
  lapply(all_samples_bc, function(nm) {

    # Per-mouse output directories
    tbl_dir <- here::here("results/tables/clustering", nm)
    fig_dir <- here::here("results/figures/clustering", nm)
    if (!dir.exists(tbl_dir)) dir.create(tbl_dir, recursive = TRUE)
    if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

    f <- filter_for_clustering(nm)

    out_path <- file.path(tbl_dir, paste0(nm, "_filtered.csv"))
    write_csv(f, out_path)
    message("    Saved: ", out_path)

    f
  }),
  all_samples_bc
)

message("--- 00_filter (clustering): done ---")
