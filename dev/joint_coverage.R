#!/usr/bin/env Rscript
# Issue #32: date-held-out LOO calibration check — does cross-analyte coupling
# improve empirical coverage of the combined AmsPAF?
#
# Holds out every 4th grab date, re-runs amspaf_daily() with and without
# couple_residuals, and records whether the hold-out day's "truth" (the
# deterministic centre from the full data) falls inside the 90% CI.
#
# Also runs a quick λ sensitivity check: how does ridge regularisation
# strength (λ ∈ {0.05, 0.10, 0.20}) affect coupling strength?
#
#   Rscript dev/joint_coverage.R  (first run: slow; ~2-3 min for the LOO loop)

suppressMessages({ library(dplyr); library(ggplot2); devtools::load_all(".", quiet = TRUE) })
options(hydroSense.guideline_dir = "guideline data")

CACHE   <- "test data/bs01_v3_cache.qs2"
CENTRE  <- "dev/bs01_kalman_centre.qs2"
LOO_OUT <- "dev/loo_coverage_coupled.qs2"

N_DRAWS  <- 15L    # kept low: LOO calls are expensive on B.S01
SEED     <- 42L
INTERVAL <- 0.9

cc <- qs2::qs_read(CACHE)
da <- cc$daily_args

## ── Deterministic centre (full data, point mode) ──────────────────────────────
if (file.exists(CENTRE)) {
  ctr_cache <- qs2::qs_read(CENTRE)
  ## CENTRE stores a list with ara/total; use ara (ARA path).
  centre_full <- ctr_cache$ara
} else {
  cat("Computing deterministic centre (point mode) ...\n")
  centre_full <- suppressMessages(do.call(amspaf_daily, da)) |>
    dplyr::select("date", "amspaf")
}

## ── Hold-out grab dates ───────────────────────────────────────────────────────
grab_dates <- sort(unique(as.Date(da$df$datetime)))
holdout    <- grab_dates[seq(4L, length(grab_dates), by = 7L)]
cat(sprintf("Hold-out dates: %d of %d grab dates\n", length(holdout), length(grab_dates)))

## ── LOO loop (cached) ─────────────────────────────────────────────────────────
if (file.exists(LOO_OUT)) {
  cat("Loading cached LOO results from", LOO_OUT, "\n")
  loo_rows <- qs2::qs_read(LOO_OUT)
} else {
  cat("Running LOO loop (", length(holdout), "dates) ...\n")

  loo_rows <- purrr::map(holdout, function(d) {
    truth_row <- dplyr::filter(centre_full, .data$date == .env$d)
    if (nrow(truth_row) == 0L) {
      cat(sprintf("  skip %s (not in centre)\n", d))
      return(NULL)
    }
    truth <- truth_row$amspaf[[1L]]

    da_loo         <- da
    da_loo$df      <- dplyr::filter(da$df, as.Date(.data$datetime) != .env$d)
    da_loo$ndraws  <- N_DRAWS
    da_loo$seed    <- SEED
    da_loo$return  <- "summary"
    da_loo$interval <- INTERVAL

    run_loo <- function(couple) {
      args <- da_loo
      args$couple_residuals <- couple
      suppressMessages(do.call(amspaf_daily, args))
    }

    out_indep  <- run_loo(FALSE)
    out_coup   <- run_loo(TRUE)

    row_indep <- dplyr::filter(out_indep, .data$date == .env$d)
    row_coup  <- dplyr::filter(out_coup,  .data$date == .env$d)

    if (nrow(row_indep) == 0L || nrow(row_coup) == 0L) {
      cat(sprintf("  skip %s (date not in predictive output)\n", d))
      return(NULL)
    }

    tibble::tibble(
      date           = d,
      truth          = truth,
      indep_lower    = row_indep$amspaf_lower[[1L]],
      indep_upper    = row_indep$amspaf_upper[[1L]],
      coupled_lower  = row_coup$amspaf_lower[[1L]],
      coupled_upper  = row_coup$amspaf_upper[[1L]],
      indep_covered  = truth >= row_indep$amspaf_lower[[1L]] &&
                       truth <= row_indep$amspaf_upper[[1L]],
      coupled_covered = truth >= row_coup$amspaf_lower[[1L]] &&
                        truth <= row_coup$amspaf_upper[[1L]],
      indep_width    = row_indep$amspaf_upper[[1L]] - row_indep$amspaf_lower[[1L]],
      coupled_width  = row_coup$amspaf_upper[[1L]]  - row_coup$amspaf_lower[[1L]]
    )
  }) |>
    purrr::compact() |>
    dplyr::bind_rows()

  qs2::qs_save(loo_rows, LOO_OUT)
  cat("Saved LOO results to", LOO_OUT, "\n")
}

## ── Coverage summary ──────────────────────────────────────────────────────────
n <- nrow(loo_rows)
cat(sprintf("\n== LOO coverage (n=%d hold-out dates, interval=%.0f%%) ==\n", n, 100 * INTERVAL))
cat(sprintf("                         Independent   Coupled\n"))
cat(sprintf("Empirical coverage:      %5.1f%%        %5.1f%%    (target %.0f%%)\n",
            100 * mean(loo_rows$indep_covered),
            100 * mean(loo_rows$coupled_covered),
            100 * INTERVAL))
cat(sprintf("Mean CI width:           %7.2f       %7.2f\n",
            mean(loo_rows$indep_width, na.rm = TRUE),
            mean(loo_rows$coupled_width, na.rm = TRUE)))
cat(sprintf("Width ratio (cpl/indep): %7.3f\n",
            mean(loo_rows$coupled_width, na.rm = TRUE) /
            mean(loo_rows$indep_width, na.rm = TRUE)))

## ── λ sensitivity ─────────────────────────────────────────────────────────────
## Fit the target model once to extract anchor residuals for the sensitivity check.
cat("\n== Ridge λ sensitivity (mean |off-diag R|, min eigenvalue before nearPD) ==\n")

tm <- suppressMessages(
  fit_target_model(da$df, da$reference_model)
)
analytes <- names(tm$models)

## Base case: use .anchor_residual_cor() with the hardcoded λ=0.10.
cor_base <- hydroSense:::.anchor_residual_cor(tm, analytes)

## Build the wide anchor-residual matrix once (reused across λ values).
extract_S <- function(nm) {
  m    <- tm$models[[nm]]
  da_m <- m$d_anchors
  anch <- if (!is.null(da_m) && nrow(da_m) >= 2L) da_m else m$anchors
  if (is.null(anch) || nrow(anch) == 0L) return(NULL)
  data.frame(date = as.Date(anch$date), S = anch$S)
}
anch_list <- stats::setNames(lapply(analytes, extract_S), analytes)

make_df <- function(nm) {
  df  <- anch_list[[nm]]
  out <- if (is.null(df)) data.frame(date = as.Date(character(0)), stringsAsFactors = FALSE) else
    data.frame(date = as.Date(df$date), stringsAsFactors = FALSE)
  out[[nm]] <- if (is.null(df)) numeric(0) else df$S
  out
}
wide_list <- lapply(analytes, make_df)
S_wide <- Reduce(function(a, b) merge(a, b, by = "date", all = TRUE), wide_list)
S_mat  <- as.matrix(S_wide[, analytes, drop = FALSE])
R_hat  <- suppressWarnings(stats::cor(S_mat, use = "pairwise.complete.obs"))
R_hat[is.na(R_hat)] <- 0
diag(R_hat) <- 1
p <- length(analytes)

lambdas <- c(0.05, 0.10, 0.20)
cat(sprintf("%-6s  mean|off-diag|  min_eig_before_nearPD\n", "lambda"))
for (lam in lambdas) {
  R_ridge <- (1 - lam) * R_hat + lam * diag(p)
  min_eig <- min(eigen(R_ridge, symmetric = TRUE, only.values = TRUE)$values)
  off     <- R_ridge[upper.tri(R_ridge)]
  cat(sprintf("%.2f    %.4f          %.6f\n", lam, mean(abs(off)), min_eig))
}
