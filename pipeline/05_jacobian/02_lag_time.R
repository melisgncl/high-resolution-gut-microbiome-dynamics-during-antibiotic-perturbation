# =============================================================================
# Title:   02_lag_time.R — Paenibacillaceae establishment time (Baranyi model)
# Author:  Melis Gencel
# Date:    2026-03-20
# Input:   results/tables/16S/family/<nm>_family.csv
#          metadata/0_config.R
# Output:  results/tables/DCM/lag_times.csv
#          results/figures/DCM/lag_time_cohort1.png/pdf
#          results/figures/DCM/lag_time_cohort2.png/pdf
# Notes:   Exact archive methodology (Gencel 2020):
#          - Model: Baranyi/Buchanan 3-phase growth model on log10 scale.
#          - Fitting: minpack.lm::nlsLM with grid search over starting values
#            lag in {1,...,8}, mumax in seq(0.01, 1, 0.05). First convergence wins.
#          - The Baranyi formula is defined inline (no nlsMicrobio dependency).
#          - Output: lag (rounded to integer), mumax, LOG10N0, LOG10Nmax per mouse.
#          - Time axis labels from config breaks/labels hashes.
# =============================================================================

library(here)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(minpack.lm)

source(here::here("metadata/0_config.R"))


# =============================================================================
# BARANYI MODEL (inline — no nlsMicrobio dependency)
# =============================================================================

# Baranyi/Buchanan 3-phase kinetics on log10 scale.
# Parameters: lag, mumax, LOG10N0, LOG10Nmax
# Reference: Baranyi & Roberts (1994), Int J Food Microbiol
baranyi_formula <- LOG10N ~ LOG10Nmax +
  log10(
    (-1 + exp(mumax * lag) + exp(mumax * t)) /
    (exp(mumax * t) - 1 + exp(mumax * lag) * 10^(LOG10Nmax - LOG10N0))
  )


# =============================================================================
# FITTING FUNCTION
# =============================================================================

# Fits the Baranyi model to Paenibacillaceae abundance for one mouse.
# Returns list: plot, lag_estimate, df (raw data), fitted_df (dense fitted curve).
# remove_zero_first: if TRUE, drops the first timepoint if abundance == 0
#   (applied to m1 in archive where first sample had 0 Paeni).
find_lag_time <- function(nm, remove_zero_first = FALSE) {
  df <- read_csv(
    here::here("results/tables/16S/family", paste0(nm, "_family.csv")),
    show_col_types = FALSE
  ) %>%
    filter(Sample == nm, Family == dcm_hgt_target) %>%
    select(t = Time, Abundance.family)

  if (nrow(df) == 0) {
    warning("  ", nm, ": no ", dcm_hgt_target, " data found.")
    return(NULL)
  }

  if (remove_zero_first) {
    min_t <- min(df$t)
    if (df$Abundance.family[df$t == min_t] == 0)
      df <- df[df$t != min_t, ]
  }

  df$LOG10N <- log10(df$Abundance.family + 1e-6)
  df        <- na.omit(df)

  LOG10N0_start  <- min(df$LOG10N)
  LOG10Nmax_start <- max(df$LOG10N) + 0.1

  # Robust grid search: scan all starting values, reject non-physical fits
  # (lag < 0, LOG10Nmax above ~1 i.e. >100% abundance, mumax <= 0), keep the
  # lowest residual-SS fit. Replaces the old first-convergence-wins approach,
  # which landed on degenerate optima (e.g. m8 lag = -7, LOG10Nmax = 262).
  fit      <- NULL
  best_rss <- Inf
  lag_grid <- seq(0, 8, by = 1)
  mu_grid  <- seq(0.01, 1, by = 0.05)
  max_t       <- max(df$t)
  log_ceiling <- 0.5
  log_floor   <- max(df$LOG10N) - 1.5

  for (lag_start in lag_grid) {
    for (mu_start in mu_grid) {
      fit_try <- tryCatch(
        nlsLM(
          formula = baranyi_formula,
          data    = df,
          start   = list(lag     = lag_start,
                         mumax   = mu_start,
                         LOG10N0 = LOG10N0_start,
                         LOG10Nmax = LOG10Nmax_start),
          control = nls.lm.control(maxiter = 1000)
        ),
        error = function(e) NULL
      )
      if (is.null(fit_try)) next
      cf <- coef(fit_try)
      if (cf["lag"] < -0.5 || cf["lag"] > max_t)        next
      if (cf["LOG10Nmax"] > log_ceiling ||
          cf["LOG10Nmax"] < log_floor)                  next
      if (cf["mumax"] <= 0)                             next
      rss <- sum(residuals(fit_try)^2)
      if (rss < best_rss) { best_rss <- rss; fit <- fit_try }
    }
  }

  if (is.null(fit)) {
    warning("  ", nm, ": Baranyi model did not converge to a physical fit.")
    return(NULL)
  }

  cf              <- coef(fit)
  lag_est         <- round(cf["lag"])
  mumax_est       <- cf["mumax"]
  message("    ", nm, ": lag = ", lag_est, "  mumax = ", round(mumax_est, 3))

  # Dense fitted curve (archive uses manual formula since predict() was unreliable)
  time_seq <- seq(min(df$t), max(df$t), length.out = 100)
  pred_log10N <- cf["LOG10Nmax"] +
    log10(
      (-1 + exp(cf["mumax"] * cf["lag"]) + exp(cf["mumax"] * time_seq)) /
      (exp(cf["mumax"] * time_seq) - 1 +
         exp(cf["mumax"] * cf["lag"]) * 10^(cf["LOG10Nmax"] - cf["LOG10N0"]))
    )
  fitted_df <- data.frame(Time = time_seq, Log10N = pred_log10N)

  # Time axis for this mouse's cohort type
  tt <- dcm_time_type[nm]
  eff_breaks <- breaks.hash[[tt]]
  eff_labels <- labels.hash[[tt]]

  p <- ggplot() +
    geom_point(data = df, aes(x = t, y = LOG10N),
               colour = "black", size = 2.5) +
    geom_line(data = df, aes(x = t, y = LOG10N),
              colour = "purple", linewidth = 1) +
    geom_line(data = fitted_df, aes(x = Time, y = Log10N),
              colour = "blue", linewidth = 1.2) +
    geom_vline(xintercept = lag_est, colour = "darkgreen",
               linetype = "dashed", linewidth = 1) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("mumax = ", round(mumax_est, 3)),
             hjust = 1.1, vjust = 2, size = 4, colour = "black") +
    annotate("text", x = Inf, y = Inf,
             label = paste0("lag = ", lag_est),
             hjust = 1.1, vjust = 4, size = 4, colour = "darkgreen") +
    labs(title = nm,
         x     = "Time",
         y     = expression(log[10](N))) +
    ylim(c(-7, 0)) +
    scale_x_continuous(breaks = eff_breaks, labels = eff_labels) +
    theme_Publication(base_size = 12, aspect.ratio = 0.75)

  list(
    plot          = p,
    lag_estimate  = as.integer(lag_est),
    mumax         = as.numeric(mumax_est),
    LOG10N0       = as.numeric(cf["LOG10N0"]),
    LOG10Nmax     = as.numeric(cf["LOG10Nmax"]),
    df            = df,
    fitted_df     = fitted_df
  )
}


# =============================================================================
# RUN ALL COLONISED MICE
# =============================================================================

message("--- 02_lag_time: fitting Baranyi model for ", dcm_hgt_target, " ---")

all_samples_bc <- c(cohort_colonized, cohort_colonized_2)

# m1 starts with 0 Paeni at first timepoint → remove_zero_first = TRUE
# (archive: p5=find_lag_time(..., TRUE) for m1)
zero_first_mice <- c("m1")

results <- list()
for (nm in all_samples_bc) {
  message("  ", nm, ": fitting ...")
  results[[nm]] <- tryCatch(
    find_lag_time(nm, remove_zero_first = nm %in% zero_first_mice),
    error = function(e) {
      message("  ERROR in ", nm, ": ", conditionMessage(e))
      NULL
    }
  )
}

# Remove failed mice
results <- Filter(Negate(is.null), results)


# =============================================================================
# SAVE LAG TIME TABLE
# =============================================================================

lag_table <- bind_rows(lapply(names(results), function(nm) {
  r <- results[[nm]]
  data.frame(
    mouse     = nm,
    lag       = r$lag_estimate,
    mumax     = r$mumax,
    LOG10N0   = r$LOG10N0,
    LOG10Nmax = r$LOG10Nmax
  )
}))

out_tbl <- here::here("results/tables/DCM")
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
write_csv(lag_table, file.path(out_tbl, "lag_times.csv"))
message("  Saved: lag_times.csv")
print(lag_table)


# =============================================================================
# FIGURES — 4-PANEL PER COHORT
# =============================================================================

message("  Building cohort panels ...")

make_panel <- function(mice, stem) {
  plots <- lapply(mice, function(nm) {
    if (!is.null(results[[nm]])) results[[nm]]$plot else NULL
  })
  plots <- Filter(Negate(is.null), plots)
  if (length(plots) == 0) return(invisible(NULL))
  p_panel <- wrap_plots(plots, nrow = 1) +
    plot_annotation(
      title = paste0(dcm_hgt_target, " \u2014 Baranyi lag time"),
      theme = theme(plot.background = element_rect(fill = "white", colour = NA))
    )
  save_fig_dcm(p_panel, stem, nm = NULL, w = 26, h = 7.5)
}

make_panel(cohort_colonized,   "lag_time_cohort1")
make_panel(cohort_colonized_2, "lag_time_cohort2")

message("--- 02_lag_time: done ---")
