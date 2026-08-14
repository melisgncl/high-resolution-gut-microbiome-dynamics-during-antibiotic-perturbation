# =============================================================================
# Title:   01_import.R — Import raw barcode clustering data (8 samples)
# Author:  Melis Gencel
# Date:    2026-03-19
# Input:   data/raw/barcodes/processed_sample/M/M{5-8}/Sample_M*_clustering.txt
#          data/raw/barcodes/processed_sample/NM/M{5-8}/Sample_M*_clustering.txt
#          data/raw/barcodes/processed_sample/{M,NM}/M*/M*_cluster.csv
#          metadata/0_config.R  (bc_file_map, bc_noise_ids, bc_noise_times)
# Output:  samples_bc_raw — named list of 8 tibbles (m1-m8), in memory
#            columns: Time (int), ID (int), Reads (int), Center (chr), Sample (chr)
# Notes:   Controls (c_m1-c_m4) have no barcode data — this module covers
#          colonized_1 (m1-m4) and colonized_2 (m5-m8) cohorts only.
#          Noise filters documented in metadata/0_config.R Section 10b.
# =============================================================================

library(here)
library(readr)
library(dplyr)

source(here::here("metadata/0_config.R"))


# =============================================================================
# HELPER — RESOLVE FILE PATHS FOR ONE SAMPLE
# =============================================================================

get_bc_paths <- function(sample_id) {
  info      <- bc_file_map[[sample_id]]
  base_dir  <- if (info$cohort == "M") path_bc_m else path_bc_nm
  sample_dir <- file.path(base_dir, info$file)
  list(
    clustering = file.path(sample_dir, paste0("Sample_", info$file, "_clustering.txt")),
    cluster    = file.path(sample_dir, paste0(info$file, "_cluster.csv"))
  )
}


# =============================================================================
# READER — IMPORT AND FILTER ONE SAMPLE
# =============================================================================

import_sample_bc <- function(sample_id) {
  paths <- get_bc_paths(sample_id)

  # --- validate files exist ---
  for (p in unlist(paths)) {
    if (!file.exists(p)) stop("File not found: ", p)
  }

  # --- read clustering data (Time, Reads, ID) ---
  df <- read_tsv(paths$clustering, show_col_types = FALSE,
                 col_types = cols(Time = col_integer(),
                                  Reads = col_integer(),
                                  ID    = col_integer()))

  # --- read cluster lookup (Cluster.ID → barcode sequence Center) ---
  cluster_map <- read_csv(paths$cluster, show_col_types = FALSE,
                          col_types = cols(Cluster.ID = col_integer(),
                                           Center     = col_character()))

  # --- apply noise filters (IDs) ---
  bad_ids <- bc_noise_ids[[sample_id]]
  if (!is.null(bad_ids) && length(bad_ids) > 0) {
    n_before <- nrow(df)
    df <- df %>% filter(!ID %in% bad_ids)
    message("  ", sample_id, ": removed ", n_before - nrow(df),
            " rows (noise barcode IDs)")
  }

  # --- apply noise filters (timepoints) ---
  bad_times <- bc_noise_times[[sample_id]]
  if (!is.null(bad_times) && length(bad_times) > 0) {
    n_before <- nrow(df)
    df <- df %>% filter(!Time %in% bad_times)
    message("  ", sample_id, ": removed ", n_before - nrow(df),
            " rows (noise timepoints: ", paste(bad_times, collapse = ", "), ")")
  }

  # --- merge to add barcode sequence ---
  df <- df %>%
    left_join(cluster_map, by = c("ID" = "Cluster.ID")) %>%
    mutate(Sample = sample_id)

  # --- sanity checks ---
  n_na_center <- sum(is.na(df$Center))
  if (n_na_center > 0) {
    warning(sample_id, ": ", n_na_center,
            " rows have no matching Center sequence in cluster CSV.")
  }

  df %>% select(Sample, Time, ID, Reads, Center)
}


# =============================================================================
# LOAD ALL 8 SAMPLES
# =============================================================================

load_all_bc <- function() {
  message("--- 01_import (barcode): loading all samples ---")

  bc_samples <- c(cohort_colonized, cohort_colonized_2)   # m1-m8; no barcodes for controls

  samples_bc_raw <- lapply(bc_samples, function(nm) {
    message("  Loading ", nm, " ...")
    import_sample_bc(nm)
  })
  names(samples_bc_raw) <- bc_samples

  message("Loaded ", length(samples_bc_raw), " samples: ",
          paste(names(samples_bc_raw), collapse = ", "))
  samples_bc_raw
}

samples_bc_raw <- load_all_bc()


# =============================================================================
# SANITY CHECKS
# =============================================================================

message("\n--- Sanity checks ---")
for (nm in names(samples_bc_raw)) {
  df <- samples_bc_raw[[nm]]
  message(nm, ": ", nrow(df), " rows | ",
          length(unique(df$ID)), " unique barcodes | ",
          "timepoints: ", paste(sort(unique(df$Time)), collapse = " "))
}

message("--- 01_import (barcode): done ---")
