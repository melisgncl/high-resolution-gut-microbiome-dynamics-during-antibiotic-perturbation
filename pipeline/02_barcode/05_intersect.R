# =============================================================================
# Title:   05_intersect.R — Barcode sharing within and between cohorts
# Author:  Melis Gencel
# Date:    2026-03-19
# Input:   samples_bc — named list from 02_process.R
#          metadata/0_config.R
# Output:  results/tables/barcodes/top_barcodes_per_sample.csv
#          results/tables/barcodes/barcode_sharing_matrix.csv
#          results/figures/barcodes/intersect_heatmap.pdf/.png
#          results/figures/barcodes/intersect_sharing_bar.pdf/.png
# Notes:   Top N barcodes per sample defined by max_freq rank.
#          Intersection analysis uses presence/absence heatmap (ggplot2 only)
#          and a bar chart of sharing counts — no UpSetR dependency.
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

source(here::here("metadata/0_config.R"))
source(here::here("analysis/02_barcode/02_process.R"))

N_TOP <- 2000L   # top barcodes per sample by max_freq (matches original paper)


# =============================================================================
# GET TOP N BARCODES PER SAMPLE
# =============================================================================

message("--- 05_intersect: selecting top ", N_TOP, " barcodes per sample ---")

top_per_sample <- bind_rows(
  lapply(names(samples_bc), function(nm) {
    samples_bc[[nm]] %>%
      distinct(ID, Center, max_freq, Sample) %>%
      slice_max(max_freq, n = N_TOP, with_ties = FALSE) %>%
      mutate(Sample = nm)
  })
)

head(top_per_sample)


# =============================================================================
# SAVE TOP BARCODES TABLE
# =============================================================================

out_dir_tables <- file.path(path_tables, "barcodes")
if (!dir.exists(out_dir_tables)) dir.create(out_dir_tables, recursive = TRUE)

write_csv(top_per_sample,
          file.path(out_dir_tables, "top_barcodes_per_sample.csv"))
message("Saved: ", file.path(out_dir_tables, "top_barcodes_per_sample.csv"))


# =============================================================================
# PRESENCE / ABSENCE MATRIX
# =============================================================================

# Wide matrix: rows = unique Center sequences, cols = samples, values = 0/1
sharing_mat <- top_per_sample %>%
  select(Center, Sample) %>%
  distinct() %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = Sample, values_from = present,
              values_fill = 0L)

write_csv(sharing_mat,
          file.path(out_dir_tables, "barcode_sharing_matrix.csv"))
message("Saved: ", file.path(out_dir_tables, "barcode_sharing_matrix.csv"))

# Long format for plotting
sharing_long <- sharing_mat %>%
  pivot_longer(-Center, names_to = "Sample", values_to = "Present") %>%
  mutate(
    Sample = factor(Sample, levels = c(cohort_colonized, cohort_colonized_2)),
    Cohort = case_when(
      Sample %in% cohort_colonized   ~ "colonized",
      Sample %in% cohort_colonized_2 ~ "colonized_2"
    )
  )

# Summarise: for each barcode, how many mice per cohort carry it
sharing_summary <- sharing_mat %>%
  mutate(
    n_colonized   = rowSums(select(., all_of(cohort_colonized))),
    n_colonized_2 = rowSums(select(., all_of(cohort_colonized_2))),
    n_total       = n_colonized + n_colonized_2
  )

message("\n--- Sharing summary (top ", N_TOP, " barcodes per sample) ---")
message("Total unique barcodes across all mice: ", nrow(sharing_mat))
message("Shared in >= 2 colonized_1 mice: ",
        sum(sharing_summary$n_colonized >= 2))
message("Shared in >= 2 colonized_2 mice: ",
        sum(sharing_summary$n_colonized_2 >= 2))
message("Shared between both cohorts (>= 1 each): ",
        sum(sharing_summary$n_colonized >= 1 & sharing_summary$n_colonized_2 >= 1))


# =============================================================================
# PLOT 1 — HEATMAP: presence / absence for barcodes shared in >= 2 mice
# =============================================================================

# Subset to barcodes present in at least 2 mice (interesting)
shared_centers <- sharing_summary %>%
  filter(n_total >= 2) %>%
  pull(Center)


message("Barcodes shared in >= 2 mice: ", length(shared_centers),
        " (shown in heatmap)")

heatmap_df <- sharing_long %>%
  filter(Center %in% shared_centers)

# Order barcodes by total sharing count
bc_order <- sharing_summary %>%
  filter(n_total >= 2) %>%
  arrange(desc(n_total), desc(n_colonized)) %>%
  pull(Center)

sample_order <- c(cohort_colonized, cohort_colonized_2)

heatmap_df <- heatmap_df %>%
  mutate(
    Center = factor(Center, levels = bc_order),
    Sample = factor(Sample, levels = sample_order)
  )

# Cap display at top 200 for visual clarity
if (length(bc_order) > 200) {
  heatmap_df <- heatmap_df %>%
    filter(Center %in% bc_order[seq_len(200)])
  message("(Heatmap shows top 200 shared barcodes)")
}

p_heatmap <- ggplot(heatmap_df,
                    aes(x = Sample, y = Center, fill = factor(Present))) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = c("0" = "#f0f0f0", "1" = "#252525"),
                    labels = c("0" = "Absent", "1" = "Present"),
                    name   = "") +
  facet_grid(. ~ Cohort, scales = "free_x", space = "free_x") +
  labs(title = paste0("Top-", N_TOP, " barcode sharing across mice"),
       x = NULL, y = "Barcode (Center sequence)") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.y     = element_blank(),
    axis.ticks.y    = element_blank(),
    strip.text      = element_text(face = "bold"),
    legend.position = "bottom"
  )

save_fig_bc(p_heatmap, "intersect_heatmap", w = 10, h = 12)


# =============================================================================
# PLOT 2 — BAR CHART: sharing counts by number of mice
# =============================================================================

# Count how many barcodes are shared in exactly k colonized_1 / k colonized_2 mice
count_sharing <- function(vec, label) {
  tbl <- table(vec)
  tibble(n_mice = as.integer(names(tbl)),
         count  = as.integer(tbl),
         Cohort = label) %>%
    filter(n_mice >= 1)
}

bar_df <- bind_rows(
  count_sharing(sharing_summary$n_colonized,   "colonized"),
  count_sharing(sharing_summary$n_colonized_2, "colonized_2")
) %>%
  mutate(Cohort = factor(Cohort, levels = c("colonized", "colonized_2")))

pal_cohort <- c(colonized = c1, colonized_2 = c2)

p_bar <- ggplot(bar_df, aes(x = factor(n_mice), y = count, fill = Cohort)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = pal_cohort,
                    labels = c(colonized   = "Colonized 1 (m1-m4)",
                               colonized_2 = "Colonized 2 (m5-m8)")) +
  labs(
    title = paste0("Barcode sharing within cohorts (top-", N_TOP, " per mouse)"),
    x     = "Number of mice sharing barcode",
    y     = "Number of unique barcodes",
    fill  = "Cohort"
  ) +
  theme_Publication(base_size = 14, aspect.ratio = 0.75) +
  theme(legend.position = "right")

save_fig_bc(p_bar, "intersect_sharing_bar", w = 10, h = 6)

message("--- 05_intersect (barcode): done ---")
