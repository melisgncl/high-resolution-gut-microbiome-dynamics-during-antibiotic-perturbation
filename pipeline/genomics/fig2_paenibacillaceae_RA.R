# ============================================================
# Figure 2: Paenibacillaceae relative abundance dynamics
# Author: Melis Gencel
# Date: 2026-03-23
# Input: results/tables/16S/family/{m1..m8,c_m1..c_m4}_family.csv
# Output: results/figures/genomics/fig2_paenibacillaceae_RA.png + .pdf
# ============================================================

library(tidyverse)
library(patchwork)

# Source project config for theme
source(here::here("metadata", "0_config.R"))

# Read all 12 family abundance CSVs
samples <- c(paste0("m", 1:8), paste0("c_m", 1:4))
all_data <- map_dfr(samples, function(s) {
  path <- here::here("results", "tables", "16S", "family", paste0(s, "_family.csv"))
  read_csv(path, show_col_types = FALSE)
})

# Filter to Paenibacillaceae
paeni <- all_data %>%
  filter(Family == "Paenibacillaceae") %>%
  mutate(
    Group = ifelse(str_starts(Sample, "c_"), "Control\n(antibiotic only)", "Colonised\n(E. coli + antibiotic)"),
    Mouse = Sample
  )

# Create plot
p <- ggplot(paeni, aes(x = Time, y = Abundance.family, colour = Mouse, group = Mouse)) +
  geom_line(linewidth = 0.8, alpha = 0.8) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~Group, scales = "free_y") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_colour_manual(
    values = c(
      m1 = "#1B9E77", m2 = "#D95F02", m3 = "#7570B3", m4 = "#E7298A",
      m5 = "#66A61E", m6 = "#E6AB02", m7 = "#A6761D", m8 = "#666666",
      c_m1 = "#FF7F00", c_m2 = "#FF4444", c_m3 = "#AA44AA", c_m4 = "#44AAAA"
    )
  ) +
  labs(
    title = "Paenibacillaceae relative abundance over time",
    subtitle = "Colonised mice reach 50-67% RA; controls max 1.4%",
    x = "Time (days)",
    y = "Relative abundance",
    colour = "Mouse"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10, colour = "grey40"),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# Save
out_dir <- here::here("results", "figures", "genomics")
ggsave(file.path(out_dir, "fig2_paenibacillaceae_RA.png"), p,
       width = 10, height = 5, dpi = 300)
ggsave(file.path(out_dir, "fig2_paenibacillaceae_RA.pdf"), p,
       width = 10, height = 5)

cat("Saved fig2_paenibacillaceae_RA.png + .pdf\n")
