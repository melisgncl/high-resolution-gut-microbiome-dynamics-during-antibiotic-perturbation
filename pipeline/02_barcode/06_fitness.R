# STATUS: TODO — fitting 500k barcodes is too slow.
# Needs optimization before running.
# Options to explore later:
#   - Subsample to top N barcodes only (e.g. top 1000)
#   - Vectorize with matrix operations instead of per-barcode lm()
#   - Use data.table for speed
#   - Fit only barcodes above a frequency threshold
#
# =============================================================================
# Title:   06_fitness.R — Per-lineage fitness estimation (selection coefficients)
# Author:  Melis Gencel
# Date:    2026-03-19
# Input:   samples_bc — named list from 02_process.R
#          metadata/0_config.R
# Output:  results/tables/barcodes/fitness_bc.csv
#          results/figures/barcodes/fitness_distribution.pdf/.png
#          results/figures/barcodes/fitness_top_lineages.pdf/.png
#
# Method:  For each barcode in the early colonisation window (time 1–TMAX_FIT),
#          fit a log-linear model: log10(Freq) ~ Time.
#          Slope = per-timepoint log10 growth rate ≈ selection coefficient s.
#          s > 0: expanding lineage; s < 0: contracting.
#          Only barcodes detected at >= 3 timepoints in window are fitted.
#          TMAX_FIT = 9 (up to ~6d post-gavage) — the initial colonisation phase
#          before community competition dominates.
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readr)

source(here::here("metadata/0_config.R"))
source(here::here("analysis/02_barcode/02_process.R"))

TMAX_FIT   <- 9L    # last timepoint included in fitness window (~6d post-gavage)
MIN_DETECT <- 3L    # minimum timepoints with Freq > 0 required for model fit


# =============================================================================
# FIT LOG-LINEAR MODEL FOR ONE BARCODE
# Returns slope (s), R², p-value; NA if insufficient data.
# =============================================================================

fit_fitness <- function(time_vec, freq_vec) {
  mask  <- freq_vec > 0
  if (sum(mask) < MIN_DETECT) return(tibble(s = NA_real_, r2 = NA_real_, pval = NA_real_))

  t_fit <- time_vec[mask]
  f_fit <- log10(freq_vec[mask])

  tryCatch({
    fit  <- lm(f_fit ~ t_fit)
    smry <- summary(fit)
    tibble(
      s    = coef(fit)[["t_fit"]],
      r2   = smry$r.squared,
      pval = coef(smry)["t_fit", "Pr(>|t|)"]
    )
  }, error = function(e) tibble(s = NA_real_, r2 = NA_real_, pval = NA_real_))
}


# =============================================================================
# COMPUTE FITNESS FOR ALL SAMPLES
# =============================================================================

message("--- 06_fitness: estimating per-lineage selection coefficients ---")
message("  Fitness window: time 1 to ", TMAX_FIT,
        " | min detections: ", MIN_DETECT)

fitness_all <- bind_rows(
  lapply(names(samples_bc), function(nm) {
    df <- samples_bc[[nm]] %>%
      filter(as.integer(as.character(Time)) <= TMAX_FIT) %>%
      mutate(Time_int = as.integer(as.character(Time)))

    # One row per barcode
    bc_list <- df %>%
      distinct(ID, Center, max_freq, hex_line, Sample)

    results <- bind_rows(
      lapply(seq_len(nrow(bc_list)), function(i) {
        bc_id <- bc_list$ID[i]
        sub   <- df %>% filter(ID == bc_id)
        fit_fitness(sub$Time_int, sub$Freq) %>%
          mutate(
            ID       = bc_id,
            Center   = bc_list$Center[i],
            max_freq = bc_list$max_freq[i],
            hex_line = bc_list$hex_line[i],
            Sample   = nm
          )
      })
    )
    results
  })
) %>%
  mutate(
    Cohort = case_when(
      Sample %in% cohort_colonized   ~ "colonized",
      Sample %in% cohort_colonized_2 ~ "colonized_2"
    ),
    Cohort = factor(Cohort, levels = c("colonized", "colonized_2")),
    sig    = !is.na(pval) & pval < 0.05   # flag significant fitness estimates
  )

str(fitness_all)
message("  Fitted: ", sum(!is.na(fitness_all$s)), " barcodes | ",
        "Significant (p < 0.05): ", sum(fitness_all$sig, na.rm = TRUE))


# =============================================================================
# SAVE TABLE
# =============================================================================

out_dir_tables <- file.path(path_tables, "barcodes")
if (!dir.exists(out_dir_tables)) dir.create(out_dir_tables, recursive = TRUE)

write_csv(fitness_all %>% mutate(ID = as.integer(as.character(ID))),
          file.path(out_dir_tables, "fitness_bc.csv"))
message("Saved: ", file.path(out_dir_tables, "fitness_bc.csv"))


# =============================================================================
# PLOT 1 — DISTRIBUTION OF SELECTION COEFFICIENTS PER COHORT
# =============================================================================

p_dist <- fitness_all %>%
  filter(!is.na(s)) %>%
  ggplot(aes(x = s, fill = Cohort, colour = Cohort)) +
  geom_density(alpha = 0.4, linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  scale_fill_manual(
    values = c(colonized = c1, colonized_2 = c2),
    labels = c(colonized = "Colonized 1 (m1-m4)", colonized_2 = "Colonized 2 (m5-m8)")
  ) +
  scale_colour_manual(
    values = c(colonized = c1, colonized_2 = c2),
    labels = c(colonized = "Colonized 1 (m1-m4)", colonized_2 = "Colonized 2 (m5-m8)")
  ) +
  facet_wrap(~ Sample, nrow = 2) +
  labs(
    title = paste0("Selection coefficient distribution (early window, t=1\u2013", TMAX_FIT, ")"),
    x     = "Selection coefficient s (log\u2081\u2080 frequency change per timepoint)",
    y     = "Density",
    fill  = "Cohort", colour = "Cohort"
  ) +
  theme_Publication(base_size = 12, aspect.ratio = 0.75) +
  theme(legend.position = "bottom")

save_fig_bc(p_dist, "fitness_distribution", w = 25, h = 12)


# =============================================================================
# PLOT 2 — TOP EXPANDING LINEAGES PER SAMPLE (s > 0, significant, top 10)
# =============================================================================

top_expanding <- fitness_all %>%
  filter(sig, s > 0) %>%
  group_by(Sample) %>%
  slice_max(s, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    label = paste0("ID:", as.integer(as.character(ID))),
    label = factor(label, levels = rev(unique(label[order(s)])))
  )

p_top <- ggplot(top_expanding,
                aes(x = s, y = label, fill = Cohort)) +
  geom_col() +
  scale_fill_manual(values = c(colonized = c1, colonized_2 = c2)) +
  facet_wrap(~ Sample, scales = "free_y", nrow = 2) +
  labs(
    title = "Top 10 fastest-expanding lineages per mouse (significant fits)",
    x     = "Selection coefficient s",
    y     = NULL,
    fill  = "Cohort"
  ) +
  theme_Publication(base_size = 11, aspect.ratio = 1.2) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 7))

save_fig_bc(p_top, "fitness_top_lineages", w = 25, h = 14)

message("--- 06_fitness (barcode): done ---")
