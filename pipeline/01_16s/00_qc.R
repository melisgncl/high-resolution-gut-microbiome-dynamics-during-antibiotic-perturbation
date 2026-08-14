# =============================================================================
# Title:   00_qc.R — Sequencing depth QC report (run before rarefaction)
# Author:  Melis Gencel
# Date:    2026-03-19
# Input:   samples_raw      — named list from 01_import.R
#          metadata/0_config.R (rarefaction_thresholds, cohort vectors)
# Output:  results/figures/16S/qc_depth.pdf/.png  — per-sample depth plot
#          results/tables/16S/qc_depth.csv         — depth table with flags
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

source(here::here("metadata/0_config.R"))
source(here::here("analysis/01_16S/01_import.R"))

# =============================================================================
# COMPUTE SEQUENCING DEPTH PER SAMPLE PER TIMEPOINT
# =============================================================================

# Handles both "Abundance" and "Reads" column names defensively.
get_depth_col <- function(df) {
  col <- intersect(c("Abundance", "Reads"), names(df))[1]
  if (is.na(col)) stop("No Abundance or Reads column found. Columns: ",
                        paste(names(df), collapse = ", "))
  col
}

depth_all <- bind_rows(
  lapply(names(samples_raw), function(nm) {
    df  <- samples_raw[[nm]]
    col <- get_depth_col(df)
    df %>%
      group_by(Time) %>%
      summarise(Depth = sum(.data[[col]]), .groups = "drop") %>%
      mutate(
        Sample    = nm,
        Cohort    = case_when(
          nm %in% cohort_colonized   ~ "colonized",
          nm %in% cohort_colonized_2 ~ "colonized_2",
          nm %in% cohort_controls    ~ "controls"
        ),
        Threshold = rarefaction_thresholds[nm]
      )
  })
) %>%
  mutate(
    Below_threshold = Depth < Threshold,
    Cohort          = factor(Cohort, levels = c("colonized", "colonized_2", "controls"))
  )

str(depth_all)
head(depth_all)


# =============================================================================
# FLAG TIMEPOINTS BELOW THRESHOLD
# =============================================================================

below <- depth_all %>% filter(Below_threshold)
if (nrow(below) > 0) {
  message("WARNING: ", nrow(below),
          " timepoint(s) below rarefaction threshold:")
  print(below %>% dplyr::select(Sample, Time, Depth, Threshold))
} else {
  message("All timepoints meet rarefaction thresholds.")
}


# =============================================================================
# SAVE QC TABLE
# =============================================================================

out_dir_tables <- file.path(path_tables, "16S")
if (!dir.exists(out_dir_tables)) dir.create(out_dir_tables, recursive = TRUE)

write.csv(depth_all,
          file.path(out_dir_tables, "qc_depth.csv"),
          row.names = FALSE)
message("Saved QC table: ", file.path(out_dir_tables, "qc_depth.csv"))


# =============================================================================
# PLOT — sequencing depth before rarefaction, faceted by cohort
# =============================================================================

if (!dir.exists(path_results_16S)) dir.create(path_results_16S, recursive = TRUE)

p_qc <- ggplot(depth_all,
               aes(x = Time, y = Depth, group = Sample, color = Sample)) +
  geom_line(linewidth = 1) +
  geom_point(aes(shape = Below_threshold), size = 3) +
  geom_hline(aes(yintercept = Threshold, color = Sample),
             linetype = "dashed", alpha = 0.35) +
  scale_y_log10(
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 4),
    labels = c("FALSE" = "Above threshold", "TRUE" = "Below threshold"),
    name   = ""
  ) +
  facet_wrap(~ Cohort, scales = "free_x") +
  labs(
    title  = "Sequencing depth per timepoint (before rarefaction)",
    y      = "Total reads (log\u2081\u2080)",
    x      = "Time index",
    color  = "Sample"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "right",
        strip.text      = element_text(face = "bold"))

for (fmt in fig_formats) {
  path <- file.path(path_results_16S, paste0("qc_depth.", fmt))
  ggsave(path, plot = p_qc, width = 14, height = 6,
         dpi    = fig_dpi,
         device = if (fmt == "pdf") cairo_pdf else NULL)
  message("Saved: ", path)
}

message("--- 00_qc: done ---")
