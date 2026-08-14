# =============================================================================
# Title:   00_qc.R — Barcode sequencing depth QC (run before any processing)
# Author:  Melis Gencel
# Date:    2026-03-19
# Input:   samples_bc_raw  — named list from 01_import.R
#          metadata/0_config.R
# Output:  results/tables/barcodes/qc_depth_bc.csv  — depth per sample per timepoint
#          results/figures/barcodes/qc_depth_bc.pdf/.png
# =============================================================================

library(here)
library(dplyr)
library(ggplot2)
library(scales)

source(here::here("metadata/0_config.R"))
source(here::here("analysis/02_barcode/01_import.R"))


# =============================================================================
# COMPUTE DEPTH PER SAMPLE PER TIMEPOINT
# =============================================================================

depth_bc <- bind_rows(
  lapply(names(samples_bc_raw), function(nm) {
    samples_bc_raw[[nm]] %>%
      group_by(Time) %>%
      summarise(Depth = sum(Reads, na.rm = TRUE),
                N_barcodes = n_distinct(ID),
                .groups = "drop") %>%
      mutate(
        Sample = nm,
        Cohort = case_when(
          nm %in% cohort_colonized   ~ "colonized",
          nm %in% cohort_colonized_2 ~ "colonized_2"
        ),
        Cohort = factor(Cohort, levels = c("colonized", "colonized_2"))
      )
  })
)

message("--- Depth summary ---")
str(depth_bc)
head(depth_bc)


# =============================================================================
# FLAG LOW-COVERAGE TIMEPOINTS (< 1000 total reads)
# =============================================================================

low_cov_threshold <- 1000L

low_cov <- depth_bc %>% filter(Depth < low_cov_threshold)
if (nrow(low_cov) > 0) {
  message("WARNING: ", nrow(low_cov),
          " timepoint(s) below ", low_cov_threshold, " reads:")
  print(low_cov %>% select(Sample, Time, Depth, N_barcodes))
} else {
  message("All timepoints above ", low_cov_threshold, " reads threshold.")
}


# =============================================================================
# SAVE QC TABLE
# =============================================================================

out_dir_tables <- file.path(path_tables, "barcodes")
if (!dir.exists(out_dir_tables)) dir.create(out_dir_tables, recursive = TRUE)

write.csv(depth_bc,
          file.path(out_dir_tables, "qc_depth_bc.csv"),
          row.names = FALSE)
message("Saved QC table: ", file.path(out_dir_tables, "qc_depth_bc.csv"))


# =============================================================================
# PLOT — depth per timepoint, faceted by cohort
# =============================================================================

pal_bc <- c(
  m1 = "#045a8d", m2 = "#2b8cbe", m3 = "#74a9cf", m4 = "#bdc9e1",
  m5 = "#7a0177", m6 = "#c51b8a", m7 = "#9138A7", m8 = "#552586"
)

p_qc <- ggplot(depth_bc,
               aes(x = Time, y = Depth, group = Sample, color = Sample)) +
  geom_line(linewidth = 1) +
  geom_point(aes(shape = Depth < low_cov_threshold), size = 3) +
  scale_color_manual(values = pal_bc) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 4),
    labels = c("FALSE" = "Above threshold", "TRUE" = "Below threshold"),
    name   = ""
  ) +
  scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  facet_wrap(~ Cohort, scales = "free_x") +
  labs(
    title = "Barcode sequencing depth per timepoint (before processing)",
    y     = "Total reads (log\u2081\u2080)",
    x     = "Time index",
    color = "Sample"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        strip.text      = element_text(face = "bold"))

save_fig_bc(p_qc, "qc_depth_bc", w = 14, h = 6)

message("--- 00_qc (barcode): done ---")
