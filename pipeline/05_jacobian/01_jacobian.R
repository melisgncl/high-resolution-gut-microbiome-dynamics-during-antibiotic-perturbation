# =============================================================================
# Title:   01_jacobian.R — Expanding-window Jacobian + eigenvalues + KPCA +
#          changepoints (geo and cpt methods)
# Author:  Melis Gencel
# Date:    2026-03-20
# Input:   results/tables/clustering/<nm>/<nm>_clustered_loess.csv
#          results/tables/16S/family/<nm>_family.csv
#          metadata/0_config.R
# Output:  results/tables/DCM/<nm>/<nm>_jacobians.rds
#          results/tables/DCM/<nm>/<nm>_jacobian_long.csv
#          results/tables/DCM/<nm>/<nm>_eigenvalues.csv
#          results/tables/DCM/<nm>/<nm>_kpca.csv
#          results/tables/DCM/<nm>/<nm>_changepoints_geo.csv
#          results/tables/DCM/<nm>/<nm>_changepoints_cpt.csv
#          results/figures/DCM/<nm>/<nm>_eigen_evolution.png
#          results/figures/DCM/<nm>/<nm>_kpca_trajectory.png
# Notes:   Exact archive methodology (Gencel 2020):
#          - Series: barcode LOESS (log10) + 16S families (LOESS-smoothed,
#            log10) on same dense grid.
#          - Controls: 16S only (no barcode), dense grid = 10x obs timepoints.
#          - Jacobian: J[j,i] = cov(dz_i/dt [1:t], z_j [1:t]), expanding window.
#          - Derivatives: central differences, diff(y,lag=2)/diff(xx,lag=2).
#          - Eigenvalues: sorted descending by conjugate magnitude |lambda|.
#          - KPCA: kernlab::kpca with rbf kernel, sigma from sigest() median.
#          - Changepoints: two separate methods saved to separate CSVs.
#            find_changepoints_geo() uses changepoint.geo::geomcp (2D joint).
#            find_changepoints_cpt() uses changepoint::cpt.np (1D per component).
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(kernlab)
library(changepoint)
library(changepoint.geo)

source(here::here("metadata/0_config.R"))

set.seed(170L)


# =============================================================================
# DATA LOADING
# =============================================================================

# Barcode LOESS → wide matrix (rows = dense grid, cols = C1, C2, ...)
# Returns list(mat, time_dense)
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
  list(mat = mat, time_dense = df$time)
}

# 16S family abundances → wide matrix (rows = timepoints, cols = families)
# Returns list(mat, time_obs)
load_16S_family <- function(nm) {
  df <- read_csv(
    here::here("results/tables/16S/family", paste0(nm, "_family.csv")),
    show_col_types = FALSE
  ) %>%
    filter(Sample == nm) %>%
    select(Family, Time, Abundance.family) %>%
    pivot_wider(id_cols = Time, names_from = Family,
                values_from = Abundance.family, values_fill = 0) %>%
    arrange(Time)

  list(mat = as.matrix(df[, -1]), time_obs = df$Time)
}

# Curated per-mouse 16S taxa (archive 7_DCM / 2_plot_DCM input).
# These wide TSVs (Time + one column per family) are the same hand-curated
# family set the original DCM used, so the Jacobian reproduces the archive.
# Colonised mice only (m1-m8); controls fall back to load_16S_family.
load_16S_curated <- function(nm) {
  df <- read_tsv(
    here::here("data/16S/taxa_long_format", paste0(nm, "_16S_taxa.tsv")),
    show_col_types = FALSE
  ) %>%
    arrange(Time)

  list(mat = as.matrix(df[, setdiff(names(df), "Time")]), time_obs = df$Time)
}


# =============================================================================
# 16S PREPROCESSING
# =============================================================================

# Drop families with >= n_zero_threshold zero timepoints (archive: countZeroes)
filter_16S_zeros <- function(mat, n_zero_threshold = dcm_zero_filter) {
  keep <- apply(mat, 2, function(col) sum(col == 0) < n_zero_threshold)
  mat[, keep, drop = FALSE]
}

# Drop families with mean abundance < dcm_mean_filter (archive: filter_columns_by_mean)
filter_16S_mean <- function(mat, threshold = dcm_mean_filter) {
  keep <- colMeans(mat) >= threshold
  mat[, keep, drop = FALSE]
}

# LOESS-smooth each 16S family (log10-transformed) to a target dense time grid.
# obs_time: actual measurement timepoints (same length as nrow(mat))
# dense_time: target grid (output length)
# Returns smoothed matrix (rows = dense_time, cols = families), log10 scale.
loess_smooth_16S <- function(mat, obs_time, dense_time,
                              span = dcm_loess_span) {
  log_mat <- log10(mat + 1e-6)
  smoothed <- vapply(seq_len(ncol(log_mat)), function(j) {
    predict(loess(log_mat[, j] ~ obs_time, span = span),
            newdata = data.frame(obs_time = dense_time), se = FALSE)
  }, numeric(length(dense_time)))
  colnames(smoothed) <- colnames(mat)
  smoothed
}


# =============================================================================
# JACOBIAN COMPUTATION
# =============================================================================

# Central differences: diff(y, lag=2) / diff(xx, lag=2)
# Returns matrix with nrow(mat) - 2 rows.
compute_derivatives <- function(mat, xx) {
  apply(mat, 2, function(y) diff(y, lag = 2) / diff(xx, lag = 2))
}

# Expanding-window Jacobian.
# J[j, i] = cov(dz_i/dt[1:t], z_j[1:t]) for each window endpoint t.
# Returns list: jacobians (list of matrices), spots (window endpoints)
compute_expanding_jacobians <- function(abundance_mat, xx,
                                        window = dcm_window) {
  deriv_mat <- compute_derivatives(abundance_mat, xx)
  n_deriv   <- nrow(deriv_mat)
  n_sp      <- ncol(abundance_mat)
  labels    <- colnames(abundance_mat)

  spots             <- seq(from = window, to = n_deriv, by = window)
  spots[length(spots)] <- n_deriv
  spots             <- unique(spots)

  jacobians <- vector("list", length(spots))
  for (k in seq_along(spots)) {
    t_end <- spots[k]
    J     <- matrix(NA_real_, n_sp, n_sp, dimnames = list(labels, labels))
    for (i in seq_len(n_sp)) {
      for (j in seq_len(n_sp)) {
        J[j, i] <- cov(deriv_mat[seq_len(t_end), i],
                       abundance_mat[seq_len(t_end), j])
      }
    }
    jacobians[[k]] <- J
  }
  list(jacobians = jacobians, spots = spots)
}

# Flatten jacobian list → long data frame (window, species_i, species_j, strength)
flatten_jacobians <- function(jacobians) {
  rows <- vector("list", length(jacobians))
  for (k in seq_along(jacobians)) {
    J    <- jacobians[[k]]
    long <- reshape2::melt(t(J))
    colnames(long) <- c("species_i", "species_j", "strength")
    long$window    <- k
    rows[[k]]      <- long
  }
  bind_rows(rows) %>%
    mutate(species_i = as.character(species_i),
           species_j = as.character(species_j)) %>%
    select(window, species_i, species_j, strength)
}


# =============================================================================
# EIGENVALUE EXTRACTION
# =============================================================================

# Returns long data frame: window, rank, R (real), I (imaginary), C (magnitude)
extract_eigenvalues <- function(jacobians) {
  rows <- vector("list", length(jacobians))
  for (k in seq_along(jacobians)) {
    J    <- jacobians[[k]]
    ev   <- eigen(J, symmetric = FALSE)
    R    <- Re(ev$values)
    Iv   <- Im(ev$values)
    Cv   <- sqrt(R^2 + Iv^2)
    ord  <- order(Cv, decreasing = TRUE)
    rows[[k]] <- data.frame(
      window = k,
      rank   = seq_along(ord),
      R      = R[ord],
      I      = Iv[ord],
      C      = Cv[ord]
    )
  }
  bind_rows(rows)
}


# =============================================================================
# KPCA
# =============================================================================

# Wide eigenspace: one row per window, cols = R_1, R_2, ..., I_1, I_2, ...
# Runs kernlab KPCA with rbf kernel, sigma from sigest() median.
# Returns data frame: window, kpca1, kpca2
compute_kpca <- function(eigen_df) {
  eigen_wide <- eigen_df %>%
    select(window, rank, R, I) %>%
    pivot_longer(c(R, I), names_to = "part", values_to = "value") %>%
    mutate(col_name = paste0(part, "_", rank)) %>%
    select(window, col_name, value) %>%
    pivot_wider(names_from = col_name, values_from = value) %>%
    arrange(window)

  df_num <- as.data.frame(eigen_wide %>% select(-window))
  df_num[is.na(df_num)] <- 0
  df_num <- as.data.frame(scale(df_num, center = TRUE, scale = FALSE))

  sigma_est <- tryCatch(
    as.numeric(kernlab::sigest(as.matrix(df_num), frac = 1)[2]),
    error = function(e) 0.1
  )

  kpca_res <- kernlab::kpca(~., data = df_num, kernel = "rbfdot",
                             kpar = list(sigma = sigma_est), features = 2)

  data.frame(
    window = eigen_wide$window,
    kpca1  = kpca_res@rotated[, 1],
    kpca2  = if (ncol(kpca_res@rotated) >= 2) kpca_res@rotated[, 2] else 0
  )
}


# =============================================================================
# CHANGEPOINTS — TWO METHODS
# =============================================================================

# Method 1: changepoint.geo::geomcp — joint 2D changepoint detection
# Pools ang.cpts + dist.cpts across all penalty types and nquantiles.
find_changepoints_geo <- function(kpca_df) {
  penalties       <- c("MBIC", "SIC", "BIC", "Hannan-Quinn")
  nq_range        <- seq(2, ceiling(nrow(kpca_df)))
  xy              <- as.matrix(kpca_df[, c("kpca1", "kpca2")])
  all_cpts        <- integer(0)

  for (penalty in penalties) {
    for (nq in nq_range) {
      res <- tryCatch(
        changepoint.geo::geomcp(xy, penalty = penalty, pen.value = 10,
                                 test.stat = "Empirical", nquantiles = nq),
        error = function(e) NULL
      )
      if (!is.null(res)) {
        all_cpts <- c(all_cpts,
                      unique(sort(c(res@ang.cpts, res@dist.cpts))))
      }
    }
  }

  if (length(all_cpts) == 0)
    return(data.frame(changepoint = integer(), frequency = integer(),
                      score = numeric()))

  freq_tbl <- table(all_cpts)
  out <- data.frame(
    changepoint = as.integer(names(freq_tbl)),
    frequency   = as.integer(freq_tbl)
  )
  out$score <- out$frequency / max(out$frequency)
  out[order(-out$frequency), ]
}

# Method 2: changepoint::cpt.np — univariate non-parametric changepoints
# Runs on kpca1 and kpca2 separately, pools results.
find_changepoints_cpt <- function(kpca_df) {
  penalties <- c("MBIC", "SIC", "BIC", "Hannan-Quinn")
  nq_range  <- seq(2, ceiling(nrow(kpca_df)))
  x1        <- as.numeric(kpca_df$kpca1)
  x2        <- as.numeric(kpca_df$kpca2)
  all_cpts  <- integer(0)

  for (penalty in penalties) {
    for (nq in nq_range) {
      m1 <- tryCatch(
        changepoint::cpt.np(x1, penalty = penalty, pen.value = 10,
                             test.stat = "empirical_distribution",
                             nquantiles = nq, minseglen = 1),
        error = function(e) NULL
      )
      m2 <- tryCatch(
        changepoint::cpt.np(x2, penalty = penalty, pen.value = 10,
                             test.stat = "empirical_distribution",
                             nquantiles = nq, minseglen = 1),
        error = function(e) NULL
      )
      if (!is.null(m1)) all_cpts <- c(all_cpts, m1@cpts)
      if (!is.null(m2)) all_cpts <- c(all_cpts, m2@cpts)
    }
  }

  # CROPS penalty
  for (nq in nq_range) {
    m1 <- tryCatch(
      changepoint::cpt.np(x1, penalty = "CROPS", pen.value = c(1, 100),
                           nquantiles = nq, minseglen = 1),
      error = function(e) NULL
    )
    m2 <- tryCatch(
      changepoint::cpt.np(x2, penalty = "CROPS", pen.value = c(1, 100),
                           nquantiles = nq, minseglen = 1),
      error = function(e) NULL
    )
    if (!is.null(m1)) {
      v <- as.vector(m1@cpts.full)
      all_cpts <- c(all_cpts, v[!is.na(v)])
    }
    if (!is.null(m2)) {
      v <- as.vector(m2@cpts.full)
      all_cpts <- c(all_cpts, v[!is.na(v)])
    }
  }

  if (length(all_cpts) == 0)
    return(data.frame(changepoint = integer(), frequency = integer(),
                      score = numeric()))

  freq_tbl <- table(all_cpts)
  out <- data.frame(
    changepoint = as.integer(names(freq_tbl)),
    frequency   = as.integer(freq_tbl)
  )
  out$score <- out$frequency / max(out$frequency)
  out[order(-out$frequency), ]
}


# =============================================================================
# FIGURES
# =============================================================================

# Eigenvalue evolution in complex plane (Re vs Im, coloured by rank, per window)
plot_eigen_evolution <- function(eigen_df, nm) {
  n_windows <- max(eigen_df$window)
  p <- ggplot(eigen_df, aes(x = R, y = I, colour = factor(rank))) +
    geom_point(alpha = 0.4, size = 2) +
    geom_path(aes(group = rank), alpha = 0.5, linewidth = 0.8,
              arrow = arrow(angle = 15, length = unit(0.15, "cm"),
                            type = "closed")) +
    geom_vline(xintercept = 0, colour = "grey30") +
    geom_hline(yintercept = 0, colour = "grey30") +
    scale_colour_manual(values = eigen.colors, guide = "none") +
    labs(title = paste0(nm, " \u2014 eigenvalue evolution"),
         x = expression("Re(" * lambda * ")"),
         y = expression("Im(" * lambda * ")")) +
    theme_minimal(base_size = 10)
  p
}

# KPCA trajectory over time (kpca1 and kpca2 vs window index)
plot_kpca_trajectory <- function(kpca_df, cpt_geo, cpt_cpt, nm) {
  sig_geo <- if (nrow(cpt_geo) > 0)
    cpt_geo$changepoint[cpt_geo$score >= 0.8] else integer(0)
  sig_cpt <- if (nrow(cpt_cpt) > 0)
    cpt_cpt$changepoint[cpt_cpt$score >= 0.8] else integer(0)

  p <- ggplot(kpca_df, aes(x = window)) +
    geom_line(aes(y = kpca1), colour = "#a00000", linewidth = 1.2) +
    geom_point(aes(y = kpca1), colour = "#a00000", size = 3) +
    geom_line(aes(y = kpca2), colour = "#2066a8", linewidth = 1.2) +
    geom_point(aes(y = kpca2), colour = "#2066a8", size = 3) +
    geom_vline(xintercept = sig_geo, linetype = "dashed",
               colour = "darkred", linewidth = 0.8) +
    geom_vline(xintercept = sig_cpt, linetype = "dotted",
               colour = "darkblue", linewidth = 0.8) +
    labs(title = paste0(nm, " \u2014 KPCA trajectory"),
         subtitle = "Red dashed = geo changepoints | Blue dotted = cpt changepoints",
         x = "Time window", y = "KPCA component") +
    theme_minimal(base_size = 10)
  p
}


# =============================================================================
# PER-MOUSE RUNNER
# =============================================================================

run_jacobian <- function(nm, is_control = FALSE) {
  message("  ", nm, ": loading data ...")

  # ---- Output directories ----
  out_tbl <- here::here("results/tables/DCM", nm)
  dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
  dir.create(here::here("results/figures/DCM", nm),
             recursive = TRUE, showWarnings = FALSE)

  # ---- Build combined series matrix ----
  if (is_control) {
    # Controls: 16S only, no barcode
    s16 <- load_16S_family(nm)
    taxa_mat <- filter_16S_mean(s16$mat)
    n_obs    <- nrow(taxa_mat)
    xx       <- seq(1, n_obs, length.out = n_obs * 10 - 9)
    # Use sequential indices so LOESS training range matches the dense grid
    # (avoids extrapolation for mice with missing early timepoints).
    # Enforce minimum span so ceiling(span * n) >= 3 (required for degree-2
    # LOESS); controls have few timepoints (7-10) so span=0.2 is too small.
    span_eff <- max(dcm_loess_span, 3 / n_obs)
    smoothed <- loess_smooth_16S(taxa_mat, seq_len(n_obs), xx, span = span_eff)
    combined <- smoothed
  } else {
    # Colonised mice: barcode LOESS + LOESS-smoothed 16S (curated taxa, matches archive)
    bc  <- load_loess_wide(nm)
    s16 <- load_16S_curated(nm)

    # Dense grid in the same units as barcode time column
    xx <- bc$time_dense

    # Subset 16S to time range covered by barcode dense grid
    t_min <- min(xx)
    t_max <- max(xx)
    row_keep <- s16$time_obs >= t_min & s16$time_obs <= t_max
    taxa_mat  <- s16$mat[row_keep, , drop = FALSE]
    taxa_time <- s16$time_obs[row_keep]

    taxa_mat <- filter_16S_zeros(taxa_mat)
    smoothed <- loess_smooth_16S(taxa_mat, taxa_time, xx)
    combined <- cbind(bc$mat, smoothed)
  }

  message("    Series: ", ncol(combined),
          "  |  Dense grid: ", nrow(combined), " points")

  # ---- Jacobians ----
  message("  ", nm, ": computing expanding-window Jacobians ...")
  jac_out   <- compute_expanding_jacobians(combined, xx)
  jacobians <- jac_out$jacobians
  spots     <- jac_out$spots
  message("    ", length(jacobians), " Jacobian matrices")

  # Save as RDS (used by downstream scripts)
  saveRDS(jacobians,
          file.path(out_tbl, paste0(nm, "_jacobians.rds")))

  # Save flattened long CSV
  write_csv(
    flatten_jacobians(jacobians),
    file.path(out_tbl, paste0(nm, "_jacobian_long.csv"))
  )
  message("    Saved: Jacobians (RDS + long CSV)")

  # ---- Eigenvalues ----
  message("  ", nm, ": extracting eigenvalues ...")
  eigen_df <- extract_eigenvalues(jacobians)
  write_csv(eigen_df, file.path(out_tbl, paste0(nm, "_eigenvalues.csv")))
  message("    Saved: ", nm, "_eigenvalues.csv")

  p_eigen <- plot_eigen_evolution(eigen_df, nm)
  save_fig_dcm(p_eigen, paste0(nm, "_eigen_evolution"), nm = nm, w = 10, h = 7)

  # ---- KPCA ----
  message("  ", nm, ": running KPCA ...")
  kpca_df <- tryCatch(
    compute_kpca(eigen_df),
    error = function(e) {
      message("    KPCA failed: ", conditionMessage(e))
      data.frame(window = seq_along(jacobians), kpca1 = NA_real_, kpca2 = NA_real_)
    }
  )
  write_csv(kpca_df, file.path(out_tbl, paste0(nm, "_kpca.csv")))
  message("    Saved: ", nm, "_kpca.csv")

  # ---- Changepoints ----
  message("  ", nm, ": changepoints (geo) ...")
  cpt_geo <- tryCatch(
    find_changepoints_geo(kpca_df),
    error = function(e) {
      message("    geo failed: ", conditionMessage(e))
      data.frame(changepoint = integer(), frequency = integer(), score = numeric())
    }
  )
  write_csv(cpt_geo,
            file.path(out_tbl, paste0(nm, "_changepoints_geo.csv")))
  message("    ", nrow(cpt_geo), " geo changepoints found")

  message("  ", nm, ": changepoints (cpt) ...")
  cpt_cpt <- tryCatch(
    find_changepoints_cpt(kpca_df),
    error = function(e) {
      message("    cpt failed: ", conditionMessage(e))
      data.frame(changepoint = integer(), frequency = integer(), score = numeric())
    }
  )
  write_csv(cpt_cpt,
            file.path(out_tbl, paste0(nm, "_changepoints_cpt.csv")))
  message("    ", nrow(cpt_cpt), " cpt changepoints found")

  # ---- KPCA trajectory figure ----
  p_kpca <- plot_kpca_trajectory(kpca_df, cpt_geo, cpt_cpt, nm)
  save_fig_dcm(p_kpca, paste0(nm, "_kpca_trajectory"), nm = nm, w = 10, h = 5)

  message("  ", nm, ": done.")
  invisible(list(jacobians = jacobians, eigen_df = eigen_df,
                 kpca_df = kpca_df, cpt_geo = cpt_geo, cpt_cpt = cpt_cpt))
}


# =============================================================================
# RUN ALL MICE
# =============================================================================

message("--- 01_jacobian: per-mouse Jacobian analysis ---")

all_samples_bc <- c(cohort_colonized, cohort_colonized_2)

for (nm in all_samples_bc) {
  tryCatch(
    run_jacobian(nm, is_control = FALSE),
    error = function(e) message("  ERROR in ", nm, ": ", conditionMessage(e))
  )
}

message("--- 01_jacobian: control mice (16S only) ---")

for (nm in cohort_controls) {
  tryCatch(
    run_jacobian(nm, is_control = TRUE),
    error = function(e) message("  ERROR in ", nm, ": ", conditionMessage(e))
  )
}

message("--- 01_jacobian: done ---")
