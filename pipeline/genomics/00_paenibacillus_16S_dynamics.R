# =============================================================================
# Title:   Paenibacillaceae 16S Dynamics — Colonized vs Control Mice
# Author:  Melis Gencel
# Date:    2026-03-22
# Input:   results/tables/16S/family/<sample>_family.csv (all 12 mice)
# Output:  results/genomics/figures/paenibacillus_RA_all_mice.pdf
#          results/genomics/figures/paenibacillus_RA_all_mice.png
#          results/genomics/figures/paenibacillus_RA_controls_only.pdf
# Purpose: Phase 0.2 of execution_plan.md — test H10 (ecological release) and
#          partially H4 (intrinsic resistance). Does Paenibacillaceae expand in
#          control mice (antibiotic-only, no E. coli) as well as in colonized mice?
#          If yes: spectinomycin alone drives the expansion (H10/H4 viable).
#          If no: E. coli colonization or its interaction is required (H10 weakened).
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)

source(here::here("metadata/0_config.R"))

# ---------------------------------------------------------------------------
# 1. Load all 12 family tables and filter for Paenibacillaceae
# ---------------------------------------------------------------------------

load_paeni <- function(sample_id) {
  path <- here::here("results/tables/16S/family",
                     paste0(sample_id, "_family.csv"))
  df <- read_csv(path, show_col_types = FALSE)
  str(df)  # inspect before filtering
  df |>
    filter(Family == "Paenibacillaceae") |>
    select(Time, RA = Abundance.family, Sample)
}

paeni_list <- lapply(all_samples, load_paeni)
paeni_all  <- bind_rows(paeni_list)

dim(paeni_all)
head(paeni_all)

# ---------------------------------------------------------------------------
# 2. Annotate group membership and apply consistent factor ordering
# ---------------------------------------------------------------------------

paeni_all <- paeni_all |>
  mutate(
    Group = case_when(
      Sample %in% cohort_colonized   ~ "Colonized cohort 1 (m1–m4)",
      Sample %in% cohort_colonized_2 ~ "Colonized cohort 2 (m5–m8)",
      Sample %in% cohort_controls    ~ "Antibiotic-only controls (c_m1–c_m4)"
    ),
    Group = factor(Group, levels = c(
      "Colonized cohort 1 (m1–m4)",
      "Colonized cohort 2 (m5–m8)",
      "Antibiotic-only controls (c_m1–c_m4)"
    )),
    # Per-cohort time labels (for x-axis)
    TimeLabel = case_when(
      Sample %in% cohort_colonized   ~ Time,
      Sample %in% cohort_colonized_2 ~ Time,
      Sample %in% cohort_controls    ~ Time
    )
  )

# ---------------------------------------------------------------------------
# 3. Summary stats per group × time — for ribbon plot
# ---------------------------------------------------------------------------

paeni_summary <- paeni_all |>
  group_by(Group, Time) |>
  summarise(
    mean_RA = mean(RA, na.rm = TRUE),
    sd_RA   = sd(RA, na.rm = TRUE),
    n       = n(),
    se_RA   = sd_RA / sqrt(n),
    .groups = "drop"
  )

# ---------------------------------------------------------------------------
# 4. Panel A — individual mouse trajectories, faceted by group
# ---------------------------------------------------------------------------

group_colors <- c(
  "Colonized cohort 1 (m1–m4)"          = c1,   # blue
  "Colonized cohort 2 (m5–m8)"          = c2,   # purple
  "Antibiotic-only controls (c_m1–c_m4)" = c3   # orange
)

pA <- ggplot(paeni_all, aes(x = Time, y = RA, colour = Group, group = Sample)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 1.5, alpha = 0.9) +
  facet_wrap(~Group, ncol = 3, scales = "free_x") +
  scale_colour_manual(values = group_colors, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, NA)) +
  labs(
    title = "Paenibacillaceae relative abundance — all mice",
    x     = "Observation timepoint (index)",
    y     = "Relative abundance"
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  )

# ---------------------------------------------------------------------------
# 5. Panel B — mean ± SE ribbon for colonized vs controls (same time index)
#    Controls have fewer timepoints (10 obs) vs colonized (18-19 obs).
#    Overlay both cohorts + controls on one plot for direct comparison.
# ---------------------------------------------------------------------------

pB <- ggplot(paeni_summary,
             aes(x = Time, y = mean_RA,
                 colour = Group, fill = Group, group = Group)) +
  geom_ribbon(aes(ymin = pmax(0, mean_RA - se_RA),
                  ymax = mean_RA + se_RA),
              alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                     limits = c(0, NA)) +
  labs(
    title  = "Paenibacillaceae — mean ± SE per group",
    x      = "Observation timepoint (index)",
    y      = "Mean relative abundance",
    colour = NULL, fill = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

# ---------------------------------------------------------------------------
# 6. Save combined figure
# ---------------------------------------------------------------------------

out_dir <- here::here("results/genomics/figures")

combined <- pA / pB + plot_annotation(
  caption = paste0(
    "Phase 0.2 — Execution plan decision: does Paenibacillaceae expand in ",
    "antibiotic-only control mice?\n",
    "Controls: c_m1–c_m4 (spectinomycin only, no E. coli inoculation). ",
    "Colonized: m1–m8 (E. coli + spectinomycin)."
  )
)

ggsave(file.path(out_dir, "paenibacillus_RA_all_mice.pdf"),
       plot   = combined,
       width  = 18, height = 14,
       device = cairo_pdf)

ggsave(file.path(out_dir, "paenibacillus_RA_all_mice.png"),
       plot   = combined,
       width  = 18, height = 14,
       dpi    = 300)

message("Saved: paenibacillus_RA_all_mice.pdf / .png")

# ---------------------------------------------------------------------------
# 7. Controls-only close-up (smaller panel for easy visual inspection)
# ---------------------------------------------------------------------------

paeni_controls <- paeni_all |> filter(Group == "Antibiotic-only controls (c_m1–c_m4)")

pC <- ggplot(paeni_controls, aes(x = Time, y = RA, colour = Sample, group = Sample)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_colour_manual(values = setNames(pal3, cohort_controls)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1),
                     limits = c(0, NA)) +
  labs(
    title  = "Paenibacillaceae — antibiotic-only controls",
    x      = "Observation timepoint (index; 1=day 1, 10=day 10)",
    y      = "Relative abundance",
    colour = "Mouse"
  ) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "paenibacillus_RA_controls_only.pdf"),
       plot   = pC,
       width  = 10, height = 6,
       device = cairo_pdf)

ggsave(file.path(out_dir, "paenibacillus_RA_controls_only.png"),
       plot   = pC,
       width  = 10, height = 6,
       dpi    = 300)

message("Saved: paenibacillus_RA_controls_only.pdf / .png")

# ---------------------------------------------------------------------------
# 8. Numerical summary for inline interpretation (print to console)
# ---------------------------------------------------------------------------

cat("\n=== Paenibacillaceae RA summary by group ===\n")

paeni_all |>
  group_by(Group) |>
  summarise(
    max_RA   = max(RA, na.rm = TRUE),
    mean_RA  = mean(RA, na.rm = TRUE),
    n_mice   = n_distinct(Sample),
    n_zeros  = sum(RA == 0, na.rm = TRUE),
    n_obs    = n()
  ) |>
  print()

cat("\n=== Controls — per-mouse max Paenibacillaceae RA ===\n")

paeni_controls |>
  group_by(Sample) |>
  summarise(
    max_RA  = max(RA, na.rm = TRUE),
    last_RA = RA[which.max(Time)]
  ) |>
  print()

cat("\n=== Colonized — per-mouse max Paenibacillaceae RA ===\n")

paeni_all |>
  filter(Group != "Antibiotic-only controls (c_m1–c_m4)") |>
  group_by(Sample) |>
  summarise(
    max_RA  = max(RA, na.rm = TRUE),
    last_RA = RA[which.max(Time)]
  ) |>
  print()
