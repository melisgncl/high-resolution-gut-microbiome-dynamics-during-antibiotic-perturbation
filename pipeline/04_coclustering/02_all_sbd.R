# =============================================================================
# Title:   02_all_sbd.R — Cross-mouse global SBD heatmap + violin analysis
# Author:  Melis Gencel
# Date:    2026-03-20
# Input:   results/tables/clustering/<nm>/<nm>_clustered_loess.csv
#          results/tables/16S/family/<nm>_family.csv
#          results/tables/clustering/similarity_all_overlap.csv  (step 04)
#          metadata/0_config.R
# Output:  results/figures/coclustering/all_sbd_heatmap.pdf/.png
#          results/figures/coclustering/all_sbd_violin.pdf/.png
#          results/figures/coclustering/all_overlap_heatmap.pdf/.png
#          results/figures/coclustering/all_overlap_violin.pdf/.png
#          results/tables/coclustering/all_sbd_dist.csv
#          results/tables/coclustering/all_sbd_pairs.csv
# Notes:   Combines all 8 mice into one all-vs-all SBD distance matrix.
#          Series labels: m1.C1, m1.Akkermansiaceae, m2.C1, ...
#          Heatmap: ggplot tiles + ggdendro dendrogram + clone/taxon strip.
#          Violin: SBD distance (and overlap if available) by pair type
#          (clone-clone, clone-species, species-species) within vs between
#          the k_global HC co-clusters.
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(patchwork)
library(ggdendro)
library(dtwclust)
library(proxy)
library(ggsignif)
library(ggbeeswarm)
library(dendextend)

source(here::here("metadata/0_config.R"))

# Number of global co-clusters to cut from HC (user-configurable)
k_global <- 8L


# =============================================================================
# HELPERS  (identical logic to 01_coclustering.R)
# =============================================================================

load_loess_wide <- function(nm) {
  df <- read_csv(
    here::here("results/tables/clustering", nm,
               paste0(nm, "_clustered_loess.csv")),
    show_col_types = FALSE
  ) %>%
    filter(!is.na(loess_value)) %>%
    pivot_wider(id_cols = time, names_from = cluster,
                values_from = loess_value) %>%
    arrange(time)
  mat           <- as.matrix(df[, -1])
  colnames(mat) <- paste0("C", colnames(mat))
  mat
}

load_16S_wide <- function(nm) {
  df <- read_csv(
    here::here("results/tables/16S/family", paste0(nm, "_family.csv")),
    show_col_types = FALSE
  ) %>%
    filter(Sample == nm) %>%
    select(Family, Time, Abundance.family) %>%
    pivot_wider(id_cols = Time, names_from = Family,
                values_from = Abundance.family, values_fill = 0) %>%
    arrange(Time)
  as.matrix(df[, -1])
}

transform_taxa <- function(taxa_mat, target_len) {
  keep     <- apply(taxa_mat, 2, function(col) sum(col == 0) < 8)
  taxa_mat <- taxa_mat[, keep, drop = FALSE]
  ts_list  <- dtwclust::tslist(t(log10(taxa_mat + 1e-6)))
  dtwclust::reinterpolate(ts_list, new.length = target_len)
}

# Colour dendrogram: barcode cluster labels → cluster.colors[k], families → black
color_dend <- function(dend) {
  leaf_lab <- labels(dend)
  col_vec  <- setNames(rep("black", length(leaf_lab)), leaf_lab)
  bc_idx   <- grep("(\\.C|^C)[0-9]+$", leaf_lab)
  for (i in bc_idx) {
    k <- as.integer(regmatches(leaf_lab[i],
                               regexpr("[0-9]+$", leaf_lab[i])))
    if (!is.na(k) && k >= 1L && k <= length(cluster.colors))
      col_vec[i] <- cluster.colors[k]
  }
  labels_colors(dend) <- col_vec
  dend
}

# Classify a pair of series labels as clone-clone / clone-species / species-species
pair_type <- function(l1, l2) {
  is_clone <- function(x) grepl("(\\.C|^C)[0-9]+$", x)
  if (is_clone(l1) && is_clone(l2)) return("clone-clone")
  if (!is_clone(l1) && !is_clone(l2)) return("species-species")
  "clone-species"
}


# =============================================================================
# COMBINE ALL SERIES  (all 8 mice, with mouse-prefix labels)
# =============================================================================

message("--- 02_all_sbd: loading and combining series ---")

all_samples_bc  <- c(cohort_colonized, cohort_colonized_2)
combined_series <- list()

for (nm in all_samples_bc) {
  message("  loading ", nm, " ...")
  clust_mat   <- load_loess_wide(nm)
  taxa_mat    <- load_16S_wide(nm)
  target_len  <- nrow(clust_mat)
  taxa_series  <- transform_taxa(taxa_mat, target_len)
  clust_series <- dtwclust::tslist(t(clust_mat))
  series_nm    <- append(clust_series, taxa_series)
  names(series_nm) <- paste0(nm, ".", names(series_nm))
  combined_series  <- c(combined_series, series_nm)
}

message("  Total series: ", length(combined_series))


# =============================================================================
# GLOBAL SBD DISTANCE + AVERAGE-LINKAGE HC
# =============================================================================

message("--- 02_all_sbd: computing global SBD distance matrix ---")

sbd_all <- proxy::dist(combined_series, method = "SBD", znorm = TRUE)
attr(sbd_all, "Labels") <- names(combined_series)

sbd_mat <- as.matrix(sbd_all)
rownames(sbd_mat) <- colnames(sbd_mat) <- names(combined_series)

hc_all        <- hclust(sbd_all, method = "average")
dend_all      <- color_dend(as.dendrogram(hc_all))
ordered_names <- names(combined_series)[hc_all$order]

out_tbl <- here::here("results/tables/coclustering")
out_fig <- here::here("results/figures/coclustering")
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)

write_csv(
  tibble(series = rownames(sbd_mat), as.data.frame(sbd_mat)),
  file.path(out_tbl, "all_sbd_dist.csv")
)
message("  Saved: all_sbd_dist.csv")


# =============================================================================
# SBD HEATMAP  (ggplot + ggdendro, style matches 03_cross_sample.R)
# =============================================================================

message("--- 02_all_sbd: building SBD heatmap ---")

n <- length(combined_series)

# Panel 1: dendrogram
dend_data <- ggdendro::dendro_data(hc_all, type = "rectangle")

p_dend <- ggplot() +
  geom_segment(data = dend_data$segments,
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.3, colour = "grey30") +
  scale_x_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  ggdendro::theme_dendro() +
  theme(plot.margin = margin(4, 5, 0, 5))

# Panel 2: annotation strip — clone (mouse colour) vs taxon (grey)
mouse_col_map <- c(setNames(pal1, cohort_colonized),
                   setNames(pal2, cohort_colonized_2))

annot_df <- data.frame(
  name  = factor(ordered_names, levels = ordered_names),
  mouse = sub("\\..*", "", ordered_names),
  type  = ifelse(grepl("(\\.C)[0-9]+$", ordered_names), "clone", "taxon")
) %>%
  mutate(fill = ifelse(type == "clone", mouse_col_map[mouse], "#cccccc"))

p_annot <- ggplot(annot_df, aes(x = name, y = 1, fill = fill)) +
  geom_tile() +
  scale_fill_identity() +
  scale_x_discrete(expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(0, 5, 0, 5))

# Panel 3: heatmap tiles
sbd_long <- as.data.frame(sbd_mat) %>%
  mutate(row = rownames(sbd_mat)) %>%
  pivot_longer(-row, names_to = "col", values_to = "sbd") %>%
  mutate(
    row = factor(row, levels = rev(ordered_names)),
    col = factor(col, levels = ordered_names)
  )

p_heat <- ggplot(sbd_long, aes(x = col, y = row, fill = sbd)) +
  geom_tile(colour = "white", linewidth = 0.03) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick",
                       midpoint = 1, limits = c(0, 2), name = "SBD") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 5) +
  theme(
    axis.text.x    = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4),
    axis.text.y    = element_text(size = 4),
    legend.position = "right",
    plot.margin    = margin(0, 5, 5, 5)
  )

p_sbd_full <- (p_dend / p_annot / p_heat) +
  plot_layout(heights = c(1.5, 0.12, 5)) +
  plot_annotation(
    title    = "Cross-mouse SBD co-clustering \u2014 all mice (m1\u2013m8)",
    subtitle = "Shape-Based Distance on z-normalized LOESS + 16S series, average linkage HC",
    theme    = theme(plot.title    = element_text(size = 9, face = "bold"),
                     plot.subtitle = element_text(size = 7))
  )

for (fmt in fig_formats) {
  out_path <- file.path(out_fig, paste0("all_sbd_heatmap.", fmt))
  ggplot2::ggsave(out_path, plot = p_sbd_full,
                  width = 18, height = 20, dpi = fig_dpi,
                  limitsize = FALSE,
                  device = if (fmt == "pdf") cairo_pdf else NULL)
  message("  Saved: ", out_path)
}


# =============================================================================
# OVERLAP HEATMAP  (same HC ordering, from step 04)
# =============================================================================

message("--- 02_all_sbd: building overlap heatmap ---")

overlap_path <- here::here("results/tables/clustering/similarity_all_overlap.csv")

if (file.exists(overlap_path)) {
  ovlp_df  <- read_csv(overlap_path, show_col_types = FALSE)
  ovlp_mat <- as.matrix(ovlp_df[, -1])
  rownames(ovlp_mat) <- colnames(ovlp_mat) <- ovlp_df$cluster

  # Pad with zeros for series not in overlap matrix (16S families have no overlap)
  missing_labs <- setdiff(ordered_names, rownames(ovlp_mat))
  if (length(missing_labs) > 0) {
    n_old   <- nrow(ovlp_mat)
    n_new   <- n_old + length(missing_labs)
    pad_mat <- matrix(0, n_new, n_new,
                      dimnames = list(c(rownames(ovlp_mat), missing_labs),
                                      c(colnames(ovlp_mat), missing_labs)))
    pad_mat[seq_len(n_old), seq_len(n_old)] <- ovlp_mat
    ovlp_mat <- pad_mat
  }
  ovlp_ord <- ovlp_mat[ordered_names, ordered_names]

  ovlp_long <- as.data.frame(ovlp_ord) %>%
    mutate(row = rownames(ovlp_ord)) %>%
    pivot_longer(-row, names_to = "col", values_to = "overlap") %>%
    mutate(
      row = factor(row, levels = rev(ordered_names)),
      col = factor(col, levels = ordered_names)
    )

  p_ovlp <- ggplot(ovlp_long, aes(x = col, y = row, fill = overlap)) +
    geom_tile(colour = "white", linewidth = 0.03) +
    scale_fill_gradient(low = "white", high = "blue",
                        limits = c(0, 1), name = "Overlap\ncoef.") +
    labs(
      title    = "Barcode overlap \u2014 all mice (same HC order as SBD)",
      subtitle = "Szymkiewicz\u2013Simpson overlap coefficient",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 5) +
    theme(
      axis.text.x    = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4),
      axis.text.y    = element_text(size = 4),
      legend.position = "right"
    )

  for (fmt in fig_formats) {
    out_path <- file.path(out_fig, paste0("all_overlap_heatmap.", fmt))
    ggplot2::ggsave(out_path, plot = p_ovlp,
                    width = 18, height = 17, dpi = fig_dpi,
                    limitsize = FALSE,
                    device = if (fmt == "pdf") cairo_pdf else NULL)
    message("  Saved: ", out_path)
  }
} else {
  message("  Overlap matrix not found — skipping (run 04_cluster_similarity.R first).")
  ovlp_mat <- NULL
}


# =============================================================================
# PAIR ANALYSIS
# =============================================================================

message("--- 02_all_sbd: pair analysis ---")

co_group <- setNames(cutree(hc_all, k = k_global), names(combined_series))

n_s     <- nrow(sbd_mat)
tri_idx <- which(lower.tri(matrix(0, n_s, n_s)), arr.ind = TRUE)

pairs_df <- data.frame(
  series1 = rownames(sbd_mat)[tri_idx[, 1]],
  series2 = colnames(sbd_mat)[tri_idx[, 2]],
  sbd     = sbd_mat[tri_idx]
) %>%
  mutate(
    cogroup1     = co_group[series1],
    cogroup2     = co_group[series2],
    same_cogroup = cogroup1 == cogroup2,
    pair_type    = mapply(pair_type, series1, series2),
    relationship = ifelse(same_cogroup,
                          "within co-cluster", "between co-clusters")
  )

# Add overlap values (0 if missing — 16S pairs have no barcode overlap)
if (!is.null(ovlp_mat)) {
  pairs_df$overlap <- mapply(function(l1, l2) {
    if (l1 %in% rownames(ovlp_mat) && l2 %in% colnames(ovlp_mat))
      ovlp_mat[l1, l2]
    else 0
  }, pairs_df$series1, pairs_df$series2)
}

write_csv(pairs_df, file.path(out_tbl, "all_sbd_pairs.csv"))
message("  Saved: all_sbd_pairs.csv  (",
        nrow(filter(pairs_df, same_cogroup)), " within-cogroup pairs)")


# =============================================================================
# VIOLIN PLOTS
# =============================================================================

message("--- 02_all_sbd: violin plots ---")

do_violin <- function(df_within, y_var, y_lab, title_str, stem) {
  types_present  <- unique(df_within$pair_type)
  all_comps      <- list(c("clone-clone",   "clone-species"),
                         c("clone-species", "species-species"),
                         c("clone-clone",   "species-species"))
  valid_comps    <- Filter(function(c) all(c %in% types_present), all_comps)

  p <- ggplot(df_within, aes(x = pair_type, y = .data[[y_var]])) +
    geom_violin(trim = FALSE, fill = "#b0c4de", colour = "grey40") +
    geom_quasirandom(size = 0.4, alpha = 0.5, colour = "grey30") +
    geom_signif(
      comparisons      = valid_comps,
      test             = "wilcox.test",
      map_signif_level = TRUE,
      step_increase    = 0.12,
      textsize         = 4,
      tip_length       = 0
    ) +
    labs(
      title    = title_str,
      subtitle = paste0("Within-cogroup pairs  (k = ", k_global, " HC groups)"),
      x        = "Pair type",
      y        = y_lab
    ) +
    theme_Publication(base_size = 12, aspect.ratio = 0.75)

  for (fmt in fig_formats) {
    out_path <- file.path(out_fig, paste0(stem, ".", fmt))
    ggplot2::ggsave(out_path, plot = p, width = 8, height = 8, dpi = fig_dpi,
                    device = if (fmt == "pdf") cairo_pdf else NULL)
    message("  Saved: ", out_path)
  }
}

within_df <- filter(pairs_df, same_cogroup)

do_violin(within_df, "sbd", "SBD distance",
          "SBD distance within co-clusters \u2014 all mice",
          "all_sbd_violin")

if ("overlap" %in% names(within_df)) {
  # Only clone-involved pairs have meaningful overlap values
  do_violin(
    filter(within_df, pair_type != "species-species"),
    "overlap", "Overlap coefficient",
    "Barcode overlap within co-clusters \u2014 all mice",
    "all_overlap_violin"
  )
}

message("--- 02_all_sbd: done ---")
