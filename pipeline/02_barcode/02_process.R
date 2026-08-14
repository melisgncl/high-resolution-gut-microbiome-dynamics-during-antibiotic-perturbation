# =============================================================================
# Title:   02_process.R — Normalize barcode reads; assign plot colors
# Author:  Melis Gencel
# Date:    2026-03-19
# Input:   samples_bc_raw — named list from 01_import.R
#          metadata/all_top_max2.csv  (top barcode → hex color mapping)
#          metadata/top_colors3.csv   (ordered palette for rare barcodes)
#          metadata/0_config.R
# Output:  samples_bc — named list (in memory), one tibble per sample:
#            Sample, Time, ID, Reads, Center, Freq,
#            max_freq, mean_freq, final_freq, start_freq,
#            hex_line (specific color or "#cccccc"),
#            hex_area (specific color or random palette color)
#          results/tables/barcodes/processed/<nm>_processed.csv per sample
# Notes:   Color strategy:
#            Top barcodes  (Center in all_top_max2): use pre-assigned hex
#              → same colors as original paper figures
#            Rare barcodes (line plot): "#cccccc" (gray background)
#            Rare barcodes (area plot): random color from top_colors2 palette
#              → set.seed(42) in metadata/0_config.R ensures reproducibility
# =============================================================================

library(here)
library(dplyr)
library(readr)
library(tidyr)

source(here::here("metadata/0_config.R"))
source(here::here("analysis/02_barcode/01_import.R"))


# =============================================================================
# NORMALIZATION — frequency per timepoint (column normalization)
# =============================================================================

normalize_sample <- function(df) {
  df %>%
    group_by(Time) %>%
    mutate(
      total_reads = sum(Reads, na.rm = TRUE),
      Freq        = if_else(total_reads > 0, Reads / total_reads, 0)
    ) %>%
    ungroup() %>%
    select(-total_reads)
}


# =============================================================================
# SUMMARY STATS — per barcode across all timepoints
# =============================================================================

summarise_barcodes <- function(df) {
  df %>%
    group_by(ID, Center, Sample) %>%
    summarise(
      max_freq   = max(Freq,  na.rm = TRUE),
      mean_freq  = mean(Freq, na.rm = TRUE),
      final_freq = Freq[which.max(Time)],
      start_freq = Freq[which.min(Time)],
      .groups    = "drop"
    )
}


# =============================================================================
# COLOR ASSIGNMENT
# =============================================================================
# hex_line: top barcodes → specific paper color; rare → gray
# hex_area: top barcodes → specific paper color; rare → random palette color
#
# all.top.max and top_colors2 are loaded in metadata/0_config.R Section 10.
# set.seed(42) is called there before building long.color.list.random.

assign_colors <- function(stats_df) {
  # Build lookup: Center → hex from all_top_max2
  top_hex <- all.top.max %>%
    select(Center, hex) %>%
    distinct()

  # Join
  stats_df <- stats_df %>%
    left_join(top_hex, by = "Center")

  # hex_line: top barcodes keep their color; rare get gray
  stats_df <- stats_df %>%
    mutate(hex_line = if_else(!is.na(hex), hex, "#cccccc"))

  # hex_area: rare barcodes get a random color from the palette
  # (deterministic: set.seed(42) already called in 0_config.R)
  rare_idx  <- which(is.na(stats_df$hex))
  n_rare    <- length(rare_idx)
  if (n_rare > 0) {
    rare_cols <- long.color.list[seq_len(n_rare)]  # already shuffled via set.seed
    stats_df$hex_area <- stats_df$hex_line         # start with hex_line values
    stats_df$hex_area[rare_idx] <- rare_cols
  } else {
    stats_df$hex_area <- stats_df$hex_line
  }

  stats_df %>% select(-hex)
}


# =============================================================================
# PROCESS ONE SAMPLE
# =============================================================================

process_sample_bc <- function(df) {
  df_norm  <- normalize_sample(df)
  stats    <- summarise_barcodes(df_norm)
  stats_c  <- assign_colors(stats)

  # Join colors and stats back onto full time-series
  df_norm %>%
    left_join(stats_c %>% select(ID, max_freq, mean_freq, final_freq,
                                  start_freq, hex_line, hex_area),
              by = "ID") %>%
    # Factor-order ID by max_freq (important for stacked area ordering)
    mutate(ID = factor(ID, levels = stats_c$ID[order(stats_c$max_freq)]))
}


# =============================================================================
# PROCESS ALL SAMPLES
# =============================================================================

message("--- 02_process (barcode): normalizing and assigning colors ---")

samples_bc <- lapply(names(samples_bc_raw), function(nm) {
  message("  Processing ", nm, " ...")
  process_sample_bc(samples_bc_raw[[nm]])
})
names(samples_bc) <- names(samples_bc_raw)


# =============================================================================
# SAVE PROCESSED TABLES
# =============================================================================

out_dir_proc <- file.path(path_tables, "barcodes", "processed")
if (!dir.exists(out_dir_proc)) dir.create(out_dir_proc, recursive = TRUE)

for (nm in names(samples_bc)) {
  out_path <- file.path(out_dir_proc, paste0(nm, "_processed.csv"))
  write_csv(samples_bc[[nm]] %>% mutate(ID = as.integer(as.character(ID))),
            out_path)
}
message("Saved processed CSVs to: ", out_dir_proc)


# =============================================================================
# SANITY CHECKS
# =============================================================================

message("\n--- Sanity checks ---")
for (nm in names(samples_bc)) {
  df <- samples_bc[[nm]]
  n_top <- sum(df$hex_line != "#cccccc") / n_distinct(df$Time)
  message(nm, ": ",
          n_distinct(df$ID), " barcodes | ",
          round(n_top), " top-colored | ",
          "max freq = ", round(max(df$max_freq, na.rm = TRUE), 4))
}

message("--- 02_process (barcode): done ---")
