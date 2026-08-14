# =============================================================================
# Title:   03_stability.R — Paenibacillaceae stability contribution (Δλ_max)
# Author:  Melis Gencel
# Date:    2026-03-20
# Input:   results/tables/DCM/<nm>/<nm>_jacobians.rds          (from 01)
#          results/tables/DCM/<nm>/<nm>_changepoints_geo.csv   (from 01)
#          results/tables/DCM/lag_times.csv                    (from 02)
#          metadata/0_config.R
# Output:  results/tables/DCM/<nm>/<nm>_delta_lambda.csv
#          results/figures/DCM/<nm>/<nm>_stability.png/pdf
#          results/figures/DCM/stability_all_mice.png/pdf
# Notes:   Exact archive methodology (Gencel 2020):
#          - lambda_full:    Re(lambda_max) of full Jacobian at each window.
#          - lambda_without: Re(lambda_max) after dropping dcm_hgt_target row/col.
#          - delta_lambda = lambda_full - lambda_without.
#          - Role: Stabilizing if delta_lambda < 0, Destabilizing otherwise.
#          - Phase lines from find_changepoints_geo (score >= 0.8).
#          - Establishment time from 02_lag_time as dark-red dashed line.
#          - Additional figure: Re(lambda_max) all mice overlaid.
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)

source(here::here("metadata/0_config.R"))


# =============================================================================
# EIGENVALUE HELPERS
# =============================================================================

# Dominant eigenvalue: Re(lambda) with largest real part (> 1e-12),
# fallback to first eigenvalue. Matches archive's get_first_eigen_real().
get_lambda_max <- function(J, drop_species = NULL) {
  J <- as.matrix(J)
  if (!is.null(drop_species)) {
    keep <- !rownames(J) %in% drop_species
    J    <- J[keep, keep, drop = FALSE]
  }
  ev         <- eigen(J, symmetric = FALSE)
  real_parts <- Re(ev$values)
  max_real   <- real_parts[which.max(real_parts)]
  if (max_real > 1e-12) max_real else real_parts[1]
}


# =============================================================================
# PER-MOUSE STABILITY ANALYSIS
# =============================================================================

run_stability <- function(nm) {
  message("  ", nm, ": loading ...")

  jac_path <- here::here("results/tables/DCM", nm,
                          paste0(nm, "_jacobians.rds"))
  if (!file.exists(jac_path)) {
    message("    Jacobian RDS not found — skipping.")
    return(invisible(NULL))
  }
  jacobians <- readRDS(jac_path)

  cpt_path <- here::here("results/tables/DCM", nm,
                          paste0(nm, "_changepoints_geo.csv"))
  cpt_geo  <- if (file.exists(cpt_path))
    read_csv(cpt_path, show_col_types = FALSE) else
    data.frame(changepoint = integer(), score = numeric())

  lag_tbl <- read_csv(here::here("results/tables/DCM/lag_times.csv"),
                      show_col_types = FALSE)
  lag_est  <- lag_tbl$lag[lag_tbl$mouse == nm]
  if (length(lag_est) == 0) lag_est <- NA_integer_

  # ---- Compute lambda_full and lambda_without per window ----
  message("  ", nm, ": computing delta lambda ...")
  lambda_full    <- sapply(jacobians, get_lambda_max)
  lambda_without <- sapply(jacobians, get_lambda_max,
                            drop_species = dcm_hgt_target)
  delta_lambda   <- lambda_full - lambda_without

  df <- data.frame(
    window         = seq_along(jacobians),
    lambda_full    = lambda_full,
    lambda_without = lambda_without,
    delta_lambda   = delta_lambda,
    role           = ifelse(delta_lambda < 0, "Stabilizing", "Destabilizing")
  )

  # Save table
  out_tbl <- here::here("results/tables/DCM", nm)
  write_csv(df, file.path(out_tbl, paste0(nm, "_delta_lambda.csv")))
  message("    Saved: ", nm, "_delta_lambda.csv")

  # Phase lines: changepoints with score >= 0.8
  phase_lines <- if (nrow(cpt_geo) > 0)
    cpt_geo$changepoint[cpt_geo$score >= 0.8] else integer(0)

  # ---- Plot: delta lambda over time ----
  tt         <- dcm_time_type[nm]
  eff_breaks <- breaks.hash[[tt]]
  eff_labels <- labels.hash[[tt]]
  n_windows  <- nrow(df)

  # Build segments so color applies to the segment after each point
  df_seg <- df %>%
    mutate(
      win_end   = lead(window),
      delta_end = lead(delta_lambda),
      role_next = lead(role)
    ) %>%
    filter(!is.na(win_end))

  p <- ggplot(df_seg) +
    geom_segment(aes(x = window, xend = win_end,
                     y = delta_lambda, yend = delta_end,
                     colour = role_next),
                 linewidth = 1) +
    geom_point(aes(x = window, y = delta_lambda), size = 2, colour = "black") +
    scale_colour_manual(values = c("Stabilizing"   = "blue",
                                   "Destabilizing" = "red"),
                        name = "Role") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = phase_lines, linetype = "dotted",
               colour = "black", linewidth = 0.6) +
    {if (!is.na(lag_est))
      geom_vline(xintercept = lag_est, linetype = "dashed",
                 colour = "darkred", linewidth = 1)} +
    scale_x_continuous(
      breaks = seq_len(n_windows),
      labels = if (length(eff_breaks) >= n_windows)
        eff_labels[seq_len(n_windows)] else seq_len(n_windows),
      limits = c(1, n_windows)
    ) +
    labs(
      title    = paste0(nm, " \u2014 stability contribution of ",
                        dcm_hgt_target),
      subtitle = "Dark-red dashed = establishment time | Black dotted = geo changepoints",
      x        = "Time window",
      y        = expression(Delta * "Re(" * lambda[max] * ")")
    ) +
    theme_minimal(base_size = 11)

  save_fig_dcm(p, paste0(nm, "_stability"), nm = nm, w = 10, h = 5)

  invisible(df)
}


# =============================================================================
# RUN ALL MICE + COMBINED FIGURE
# =============================================================================

message("--- 03_stability: delta-lambda analysis ---")

all_samples_bc <- c(cohort_colonized, cohort_colonized_2)

delta_all <- list()
for (nm in all_samples_bc) {
  delta_all[[nm]] <- tryCatch(
    run_stability(nm),
    error = function(e) {
      message("  ERROR in ", nm, ": ", conditionMessage(e))
      NULL
    }
  )
}

# ---- Combined: Re(lambda_max) for all mice overlaid ----
message("  Building combined Re(lambda_max) figure ...")

mouse_col_map <- c(
  setNames(pal1, cohort_colonized),
  setNames(pal2, cohort_colonized_2)
)

combined_df <- bind_rows(lapply(names(delta_all), function(nm) {
  d <- delta_all[[nm]]
  if (is.null(d)) return(NULL)
  data.frame(window = d$window, lambda_max = d$lambda_full, mouse = nm)
}))

lag_tbl <- tryCatch(
  read_csv(here::here("results/tables/DCM/lag_times.csv"),
           show_col_types = FALSE),
  error = function(e) data.frame(mouse = character(), lag = integer())
)

if (nrow(combined_df) > 0) {
  p_all <- ggplot(combined_df, aes(x = window, y = lambda_max,
                                    colour = mouse, group = mouse)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = mouse_col_map, name = "Mouse") +
    labs(
      title    = "Re(\u03bb\u209a\u2093) over time \u2014 all colonised mice",
      subtitle = "Positive = unstable, Negative = stable",
      x        = "Time window",
      y        = expression("Re(" * lambda[max] * ")")
    ) +
    theme_minimal(base_size = 11)

  save_fig_dcm(p_all, "stability_all_mice", nm = NULL, w = 12, h = 6)
}

message("--- 03_stability: done ---")
