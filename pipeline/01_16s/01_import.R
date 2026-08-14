# =============================================================================
# Title:   01_import.R — Import and time-reformat raw 16S data
# Author:  Melis Gencel
# Date:    2026-03-19
# Input:   data/raw/16S/formatted/*_16S_abundances.csv
#          data/raw/16S/control_SPN.tsv
# Output:  samples_raw — named list of 12 tibbles, one per sample (in memory)
# =============================================================================

library(here)
library(readr)
library(dplyr)
library(tidyr)

source(here::here("metadata/0_config.R"))

# =============================================================================
# TIME REFORMATTING
# =============================================================================

# Colonized cohort (m1-m4): CSV stores time in days, multiply × 24 for hours,
# then map hours → integer index used throughout the project.
reformat_time_numeric <- function(df) {
  time_map <- c(
    "0"   =  0, "3"   =  1, "6"   =  2, "12"  =  3, "24"  =  4,
    "48"  =  5, "72"  =  6, "96"  =  7, "120" =  8, "144" =  9,
    "168" = 10, "192" = 11, "216" = 12, "240" = 13, "264" = 14,
    "288" = 15, "312" = 16, "336" = 17, "360" = 18, "384" = 19
  )
  # Time is already in hours (read_colonized_sample did the ×24). Just round and map.
  df %>%
    mutate(
      Time = as.character(round(Time)),
      Time = as.integer(time_map[Time])
    )
}

# Colonized cohort 2 (m5-m8): CSV stores time as string labels.
reformat_time_string <- function(df) {
  time_map <- c(
    "bT0" =  0, "bT1" =  1, "bT2" =  2,
    "t1"  =  3, "t2"  =  4, "t3"  =  5, "t4"  =  6,
    "d2"  =  7, "d3"  =  8, "d4"  =  9, "d5"  = 10,
    "d6"  = 11, "d7"  = 12, "d8"  = 13, "d9"  = 14,
    "d10" = 15, "d11" = 16, "d12" = 17, "d13" = 18,
    "d14" = 19, "d15" = 20
  )
  df %>%
    mutate(Time = as.integer(time_map[as.character(Time)]))
}

# =============================================================================
# SAMPLE READERS
# =============================================================================

# Colonized cohort (m1-m4 = files M5-M8)
# Time is numeric (days in CSV), needs × 24 then reindex. Drop time 0 (pre-exp).
read_colonized_sample <- function(file_name) {
  file_path <- file.path(path_16S_raw, file_name)
  read_csv(file_path,
           col_types    = cols(data = col_skip()),
           show_col_types = FALSE) %>%
    mutate(Time = Time * 24) %>%
    reformat_time_numeric() %>%
    filter(!is.na(Time), Time != 0L)
}

# Colonized cohort 2 (m5-m8 = files NM5-NM8)
# Time is string labels. Drop pre-colonisation baselines (indices 1, 2 = bT1,
# bT2), then shift remaining indices down by 2 so post-colonisation starts at 1.
# Time 0 (bT0) is retained here and dropped later in 03_filter.R.
read_colonized_2_sample <- function(file_name) {
  file_path <- file.path(path_16S_raw, file_name)
  read_csv(file_path,
           col_types    = cols(Sample = col_skip(), Species = col_skip()),
           show_col_types = FALSE) %>%
    reformat_time_string() %>%
    filter(!is.na(Time), !(Time %in% c(1L, 2L))) %>%
    mutate(Time = if_else(Time == 0L, 0L, Time - 2L))
}

# Controls (c_m1-c_m4 = subjects M5-M8 in control_SPN.tsv).
# Each control mouse has different missing timepoints — handled individually.
read_controls <- function() {
  c_M <- read_tsv(path_16S_control, show_col_types = FALSE)

  # c_m1 (M5): no time 5; time is already numeric
  c_m1 <- c_M %>%
    filter(Subject == "M5") %>%
    mutate(Time = as.numeric(as.character(Time)))

  # c_m2 (M6): time is already numeric
  c_m2 <- c_M %>%
    filter(Subject == "M6") %>%
    mutate(Time = as.numeric(as.character(Time)))

  # c_m3 (M7): has "B" (before) timepoint; drop time 0, 1, 6 (low coverage)
  c_m3 <- c_M %>%
    filter(Subject == "M7") %>%
    mutate(
      Time = if_else(Time == "B", "0", as.character(Time)),
      Time = as.numeric(Time)
    ) %>%
    filter(!Time %in% c(0, 1, 6))

  # c_m4 (M8): has "B" timepoint; drop time 0 only
  c_m4 <- c_M %>%
    filter(Subject == "M8") %>%
    mutate(
      Time = if_else(Time == "B", "0", as.character(Time)),
      Time = as.numeric(Time)
    ) %>%
    filter(!Time %in% c(0))

  list(c_m1 = c_m1, c_m2 = c_m2, c_m3 = c_m3, c_m4 = c_m4)
}

# =============================================================================
# LOAD ALL 12 SAMPLES
# =============================================================================

load_all_samples <- function() {
  message("--- 01_import: loading all samples ---")

  message("Colonized cohort (m1-m4)...")
  colonized <- list(
    m1 = read_colonized_sample("M5_16S_abundances.csv"),
    m2 = read_colonized_sample("M6_16S_abundances.csv"),
    m3 = read_colonized_sample("M7_16S_abundances.csv"),
    m4 = read_colonized_sample("M8_16S_abundances.csv")
  )

  message("Colonized cohort 2 (m5-m8)...")
  colonized_2 <- list(
    m5 = read_colonized_2_sample("NM5_16S_abundances.csv"),
    m6 = read_colonized_2_sample("NM6_16S_abundances.csv"),
    m7 = read_colonized_2_sample("NM7_16S_abundances.csv"),
    m8 = read_colonized_2_sample("NM8_16S_abundances.csv")
  )

  message("Controls (c_m1-c_m4)...")
  controls <- read_controls()

  all_samples <- c(colonized, colonized_2, controls)
  message("Loaded ", length(all_samples), " samples: ",
          paste(names(all_samples), collapse = ", "))
  all_samples
}

# =============================================================================
# RUN
# =============================================================================

samples_raw <- load_all_samples()

# Quick sanity check
for (nm in names(samples_raw)) {
  df <- samples_raw[[nm]]
  message(nm, ": ", nrow(df), " rows | timepoints: ",
          paste(sort(unique(df$Time)), collapse = " "))
}
