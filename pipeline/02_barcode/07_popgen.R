# =============================================================================
# Title:   07_popgen.R — Effective population size (Ne) + neutral drift test
# Author:  Melis Gencel
# Date:    2026-03-19
# Input:   diversity_bc — from 04_diversity.R (sourced)
#          samples_bc   — from 02_process.R  (sourced via 04_diversity.R)
#          metadata/0_config.R
# Output:  results/tables/barcodes/Ne_estimates.csv
#          results/tables/barcodes/neutral_drift_test.csv
#          results/figures/barcodes/Ne_q0_decay.pdf/.png
#          results/figures/barcodes/neutral_drift_variance.pdf/.png
#
# — Ne ESTIMATION —
#   Under Wright-Fisher drift, barcode richness decays as:
#     q_0(t) = q_0(0) * (1 - 1/(2*Ne))^t
#   Log-linearise: log(q_0) = log(q_0(0)) + t * log(1 - 1/(2*Ne))
#   Fit slope b → Ne = -0.5 / (exp(b) - 1)
#   Uses full time series per sample.
#
# — NEUTRAL DRIFT TEST —
#   Under neutrality, expected variance of Δfreq per step:
#     E[Var(Δf)] ≈ f*(1-f) / Ne
#   We compute the observed variance of per-barcode frequency changes
#   across consecutive timepoints and compare to WF expectation.
#   Overdispersion factor β = Var_obs / Var_WF:
#     β ≈ 1   → consistent with neutral drift
#     β >> 1  → excess variance → selection present
# =============================================================================

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readr)

source(here::here("metadata/0_config.R"))
source(here::here("analysis/02_barcode/04_diversity.R"))   # provides diversity_bc & samples_bc


# =============================================================================
# SECTION A — EFFECTIVE POPULATION SIZE (Ne)
# =============================================================================

message("--- 07_popgen: estimating Ne from barcode richness decay ---")

estimate_Ne <- function(div_df) {
  # Log-linear fit of q_0 ~ Time
  # Exclude Time = 0 if present; require at least 5 timepoints
  d <- div_df %>%
    filter(q_0 > 0) %>%
    arrange(Time)

  if (nrow(d) < 5) return(tibble(Ne = NA_real_, slope = NA_real_,
                                  r2 = NA_real_, n_tp = nrow(d)))

  fit  <- lm(log(q_0) ~ Time, data = d)
  smry <- summary(fit)
  b    <- coef(fit)[["Time"]]

  Ne_est <- if (b < 0 && exp(b) < 1) -0.5 / (exp(b) - 1) else NA_real_

  tibble(
    Ne    = Ne_est,
    slope = b,
    r2    = smry$r.squared,
    n_tp  = nrow(d)
  )
}

Ne_results <- diversity_bc %>%
  group_by(Sample, Cohort) %>%
  group_modify(~ estimate_Ne(.x)) %>%
  ungroup()

message("\n--- Ne estimates ---")
print(Ne_results %>% select(Sample, Cohort, Ne, slope, r2))

out_dir_tables <- file.path(path_tables, "barcodes")
if (!dir.exists(out_dir_tables)) dir.create(out_dir_tables, recursive = TRUE)

write_csv(Ne_results, file.path(out_dir_tables, "Ne_estimates.csv"))
message("Saved: ", file.path(out_dir_tables, "Ne_estimates.csv"))


# =============================================================================
# PLOT A — q_0 decay with fitted line per sample
# =============================================================================

# Predicted values for plotting
div_pred <- diversity_bc %>%
  left_join(Ne_results %>% select(Sample, slope), by = "Sample") %>%
  group_by(Sample) %>%
  mutate(
    q0_start  = q_0[which.min(Time)],
    q0_fitted = q0_start * exp(slope * (Time - min(Time)))
  ) %>%
  ungroup()

pal_all <- c(setNames(pal1, cohort_colonized),
             setNames(pal2, cohort_colonized_2))

p_decay <- ggplot(div_pred, aes(x = Time, colour = Sample)) +
  geom_line(aes(y = q_0), linewidth = 1.2) +
  geom_point(aes(y = q_0), size = 2) +
  geom_line(aes(y = q0_fitted), linetype = "dashed", linewidth = 0.8) +
  scale_colour_manual(values = pal_all) +
  facet_wrap(~ Cohort, scales = "free_x") +
  labs(
    title    = "Barcode richness decay — observed (solid) vs WF fit (dashed)",
    subtitle = "Slope of dashed line gives Ne estimate",
    x        = "Time index",
    y        = "q\u2080 (barcode richness)",
    colour   = "Sample"
  ) +
  theme_Publication(base_size = 12, aspect.ratio = 0.75) +
  theme(legend.position = "right")

save_fig_bc(p_decay, "Ne_q0_decay", w = 18, h = 7)


# =============================================================================
# SECTION B — NEUTRAL DRIFT TEST
# =============================================================================

message("\n--- 07_popgen: neutral drift overdispersion test ---")

compute_drift_test <- function(df, ne_val) {
  # For each pair of consecutive timepoints, compute Δfreq per barcode
  df_wide <- df %>%
    select(ID, Time, Freq) %>%
    mutate(Time = as.integer(as.character(Time))) %>%
    pivot_wider(names_from = Time, values_from = Freq, values_fill = 0)

  time_cols <- sort(as.integer(setdiff(names(df_wide), "ID")))

  if (length(time_cols) < 2) return(tibble(
    transition = character(0), var_obs = numeric(0),
    var_wf = numeric(0), beta = numeric(0)
  ))

  results <- lapply(seq_len(length(time_cols) - 1), function(i) {
    t1  <- time_cols[i]
    t2  <- time_cols[i + 1]
    f1  <- df_wide[[as.character(t1)]]
    f2  <- df_wide[[as.character(t2)]]

    delta     <- f2 - f1
    f_bar     <- (f1 + f2) / 2
    # WF expected variance: f*(1-f)/Ne; use mean across barcodes
    var_wf    <- mean(f_bar * (1 - f_bar), na.rm = TRUE) /
                   max(ne_val, 1, na.rm = TRUE)
    var_obs   <- var(delta, na.rm = TRUE)
    beta      <- if (!is.na(var_wf) && var_wf > 0) var_obs / var_wf else NA_real_

    tibble(transition = paste0("t", t1, "\u2192t", t2),
           var_obs    = var_obs,
           var_wf     = var_wf,
           beta       = beta)
  })
  bind_rows(results)
}

drift_results <- bind_rows(
  lapply(names(samples_bc), function(nm) {
    ne_val <- Ne_results$Ne[Ne_results$Sample == nm]
    ne_val <- if (length(ne_val) == 0 || is.na(ne_val)) 1000 else ne_val

    compute_drift_test(samples_bc[[nm]], ne_val) %>%
      mutate(
        Sample = nm,
        Cohort = if (nm %in% cohort_colonized) "colonized" else "colonized_2"
      )
  })
) %>%
  mutate(Cohort = factor(Cohort, levels = c("colonized", "colonized_2")))

message("\n--- Drift test: mean overdispersion factor (beta) per sample ---")
print(drift_results %>%
        group_by(Sample, Cohort) %>%
        summarise(mean_beta = round(mean(beta, na.rm = TRUE), 2),
                  .groups = "drop"))

write_csv(drift_results,
          file.path(out_dir_tables, "neutral_drift_test.csv"))
message("Saved: ", file.path(out_dir_tables, "neutral_drift_test.csv"))


# =============================================================================
# PLOT B — Overdispersion factor β over time
# =============================================================================

p_drift <- drift_results %>%
  filter(!is.na(beta)) %>%
  mutate(t_idx = seq_along(transition),
         Sample = factor(Sample,
                         levels = c(cohort_colonized, cohort_colonized_2))) %>%
  ggplot(aes(x = transition, y = beta,
             group = Sample, colour = Sample)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black") +
  scale_colour_manual(values = c(setNames(pal1, cohort_colonized),
                                  setNames(pal2, cohort_colonized_2))) +
  facet_wrap(~ Cohort, scales = "free_x", nrow = 1) +
  labs(
    title    = "Neutral drift test — overdispersion factor \u03b2 over time",
    subtitle = "\u03b2 \u2248 1: neutral | \u03b2 >> 1: excess variance (selection)",
    x        = "Consecutive timepoint transition",
    y        = "Overdispersion factor \u03b2 = Var\u2080\u1d47\u209b / Var\u1d42\u1da0",
    colour   = "Sample"
  ) +
  theme_Publication(base_size = 12, aspect.ratio = 0.75) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
        legend.position = "right")

save_fig_bc(p_drift, "neutral_drift_variance", w = 22, h = 8)

message("--- 07_popgen (barcode): done ---")
