# =============================================================================
# Title:   09_crosscohort.R — Cross-cohort lineage fate comparison
# Author:  Melis Gencel
# Date:    2026-03-19
# Input:   samples_bc — named list from 02_process.R
#          metadata/0_config.R
# Output:  results/tables/barcodes/crosscohort_shared.csv
#          results/figures/barcodes/crosscohort_scatter.pdf/.png
#          results/figures/barcodes/crosscohort_trajectories.pdf/.png
#          results/figures/barcodes/crosscohort_cohort_maxfreq.pdf/.png
#
# Rationale:
#   The same inoculum was used for both cohorts (colonized and colonized_2).
#   Barcodes shared across cohorts are the same founding lineages that ended
#   up in different gut environments (intact microbiome vs disrupted).
#   If a barcode thrives in one cohort but not the other, that is direct
#   evidence of microbiome-mediated competition (or release from competition).
#
# Analyses:
#   1. Per-barcode max_freq in each cohort — scatter (colonized vs colonized_2)
#   2. Spearman correlation of max_freq between cohorts
#   3. Trajectory comparison for top shared barcodes (time-series overlay)
#   4. Cohort effect: are shared lineages systematically favored in one cohort?
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

N_TRAJ_TOP <- 20L   # top shared barcodes to show trajectory plots for


# =============================================================================
# BUILD PER-BARCODE SUMMARY ACROSS ALL MICE
# =============================================================================

message("--- 09_crosscohort: building cross-cohort barcode table ---")

bc_summary_all <- bind_rows(
  lapply(names(samples_bc), function(nm) {
    samples_bc[[nm]] %>%
      distinct(ID, Center, max_freq, mean_freq, final_freq, Sample, hex_line) %>%
      mutate(
        Cohort = if (nm %in% cohort_colonized) "colonized" else "colonized_2"
      )
  })
)

# Aggregate by Center + Cohort: mean max_freq across mice in cohort
bc_cohort <- bc_summary_all %>%
  group_by(Center, Cohort) %>%
  summarise(
    mean_max_freq  = mean(max_freq,   na.rm = TRUE),
    max_max_freq   = max(max_freq,    na.rm = TRUE),
    n_mice         = n_distinct(Sample),
    .groups        = "drop"
  )

# Pivot to wide: one row per Center, columns for colonized / colonized_2
bc_wide <- bc_cohort %>%
  pivot_wider(
    names_from  = Cohort,
    values_from = c(mean_max_freq, max_max_freq, n_mice),
    values_fill = 0
  )

# Only barcodes present in BOTH cohorts (shared across the inoculum)
shared <- bc_wide %>%
  filter(n_mice_colonized >= 1, n_mice_colonized_2 >= 1)

message("  Total unique barcodes: ",  nrow(bc_wide))
message("  Shared between cohorts: ", nrow(shared))

write_csv(shared,
          file.path(path_tables, "barcodes", "crosscohort_shared.csv"))
message("Saved: ", file.path(path_tables, "barcodes", "crosscohort_shared.csv"))


# =============================================================================
# ANALYSIS 1 — SCATTER: colonized max_freq vs colonized_2 max_freq
# =============================================================================

message("--- 09_crosscohort: scatter plot ---")

# Spearman correlation
if (nrow(shared) > 5) {
  cor_test <- cor.test(shared$mean_max_freq_colonized,
                        shared$mean_max_freq_colonized_2,
                        method = "spearman", exact = FALSE)
  rho <- round(cor_test$estimate, 3)
  pv  <- signif(cor_test$p.value, 2)
  cor_label <- paste0("Spearman \u03c1 = ", rho, ", p = ", pv)
} else {
  cor_label <- "Insufficient shared barcodes for correlation"
}
message("  ", cor_label)

# Add color: use hex from all_top_max2 if available
shared_plot <- shared %>%
  left_join(all.top.max %>% select(Center, hex), by = "Center") %>%
  mutate(
    hex      = if_else(is.na(hex), "#cccccc", hex),
    is_top   = hex != "#cccccc",
    n_mice_t = n_mice_colonized + n_mice_colonized_2
  )

p_scatter <- ggplot(shared_plot,
                    aes(x = log10(mean_max_freq_colonized + 1e-8),
                        y = log10(mean_max_freq_colonized_2  + 1e-8))) +
  geom_point(data = filter(shared_plot, !is_top),
             colour = "#cccccc", size = 1, alpha = 0.4) +
  geom_point(data = filter(shared_plot,  is_top),
             aes(colour = hex), size = 2.5, alpha = 0.9) +
  scale_colour_identity() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "black", linewidth = 0.7) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = cor_label, size = 4) +
  scale_x_continuous(
    breaks = c(-8, -5, -3, -1),
    labels = c("0", "10\u207b\u2075", "10\u207b\u00b3", "10\u207b\u00b9")
  ) +
  scale_y_continuous(
    breaks = c(-8, -5, -3, -1),
    labels = c("0", "10\u207b\u2075", "10\u207b\u00b3", "10\u207b\u00b9")
  ) +
  labs(
    title    = "Cross-cohort lineage fate — shared barcodes",
    subtitle = "Same founding lineages, different gut environments",
    x        = "Mean peak frequency in colonized mice (log\u2081\u2080)",
    y        = "Mean peak frequency in colonized_2 mice (log\u2081\u2080)",
    caption  = "Dashed line = equal fate in both cohorts"
  ) +
  theme_Publication(base_size = 14, aspect.ratio = 1)

save_fig_bc(p_scatter, "crosscohort_scatter", w = 10, h = 10)


# =============================================================================
# ANALYSIS 2 — TRAJECTORY COMPARISON FOR TOP SHARED BARCODES
# =============================================================================

message("--- 09_crosscohort: trajectory overlay for top shared barcodes ---")

# Top shared barcodes by combined max_freq
top_shared_centers <- shared %>%
  mutate(combined_max = max_max_freq_colonized + max_max_freq_colonized_2) %>%
  slice_max(combined_max, n = N_TRAJ_TOP, with_ties = FALSE) %>%
  pull(Center)

traj_df <- bind_rows(
  lapply(names(samples_bc), function(nm) {
    samples_bc[[nm]] %>%
      filter(Center %in% top_shared_centers) %>%
      mutate(
        Time_int = as.integer(as.character(Time)),
        Cohort   = if (nm %in% cohort_colonized) "colonized" else "colonized_2"
      ) %>%
      select(Sample, Cohort, Center, Time_int, Freq, hex_line)
  })
) %>%
  mutate(Cohort = factor(Cohort, levels = c("colonized", "colonized_2")))

# Mean trajectory per Center per Cohort
traj_mean <- traj_df %>%
  group_by(Center, Cohort, Time_int) %>%
  summarise(mean_Freq = mean(Freq, na.rm = TRUE), .groups = "drop") %>%
  left_join(all.top.max %>% select(Center, hex), by = "Center") %>%
  mutate(hex = if_else(is.na(hex), "#888888", hex))

p_traj <- ggplot(traj_mean,
                 aes(x = Time_int, y = log10(mean_Freq + 1e-8),
                     group = Center, colour = hex)) +
  geom_line(linewidth = 0.9) +
  scale_colour_identity() +
  facet_wrap(~ Cohort, ncol = 1) +
  labs(
    title = paste0("Top-", N_TRAJ_TOP,
                   " shared lineages — mean trajectory by cohort"),
    x     = "Time index",
    y     = "Mean frequency (log\u2081\u2080)"
  ) +
  theme_Publication(base_size = 12, aspect.ratio = 0.5) +
  theme(legend.position = "none")

save_fig_bc(p_traj, "crosscohort_trajectories", w = 14, h = 12)


# =============================================================================
# ANALYSIS 3 — COHORT EFFECT: favored vs disfavored in colonized environment
# =============================================================================

message("--- 09_crosscohort: cohort effect ---")

# Log2 fold-change in max_freq: colonized vs colonized_2
shared_lfc <- shared %>%
  filter(mean_max_freq_colonized > 0, mean_max_freq_colonized_2 > 0) %>%
  mutate(
    lfc = log2(mean_max_freq_colonized / mean_max_freq_colonized_2)
  ) %>%
  left_join(all.top.max %>% select(Center, hex), by = "Center") %>%
  mutate(hex = if_else(is.na(hex), "#cccccc", hex))

p_lfc <- shared_lfc %>%
  slice_max(abs(lfc), n = 40, with_ties = FALSE) %>%
  arrange(lfc) %>%
  mutate(label = substr(Center, 1, 8),
         label = factor(label, levels = unique(label))) %>%
  ggplot(aes(x = lfc, y = label, fill = lfc > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = c1, "FALSE" = c2),
                    labels = c("TRUE" = "Favored in colonized",
                               "FALSE"= "Favored in colonized_2"),
                    name   = "") +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  labs(
    title = "Cohort effect: top 40 differentially abundant shared lineages",
    x     = "log\u2082(colonized max_freq / colonized_2 max_freq)",
    y     = "Barcode (first 8 nt shown)"
  ) +
  theme_Publication(base_size = 11, aspect.ratio = 1.8) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 7, family = "mono"))

save_fig_bc(p_lfc, "crosscohort_cohort_maxfreq", w = 10, h = 14)

message("--- 09_crosscohort (barcode): done ---")
