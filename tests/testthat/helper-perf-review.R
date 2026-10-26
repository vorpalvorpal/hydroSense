## ============================================================================
## Synthetic fixtures for the review findings pinned by test-perf-review.R
## ============================================================================
##
## Used by test-perf-review.R and the golden generator
## data-raw/make_review_goldens.R.  Builds on helper-perf-fixtures.R (the
## synthetic scenario and reference model).  Everything here is SYNTHETIC:
## fixed seeds, invented site codes, no monitoring data.
##
## R1  Batch prediction failure.  One Inf in the reference model's hydrology
##     makes predict.gam() error for every date whose hydro features are
##     non-finite (it returns NA for an NA input, but errors on Inf).  The
##     pre-refactor resolver predicted one date at a time, so only those
##     dates fell back to the static tier.
## R2  Prefix-sum precision.  time_weighted_aggregate(summary = "arith_mean")
##     on series whose earlier values dwarf the values in the current window.

## R1 -------------------------------------------------------------------------

.pfr_inf_date <- as.Date("2023-06-01")

## Daily grid of target dates spanning the whole hydro record and beyond it.
.pfr_inf_grid <- function() {
  seq(as.Date("2020-06-01"), as.Date("2024-03-01"), by = "day")
}

.pfr_inf_target <- function() {
  g <- .pfr_inf_grid()
  tibble::tibble(sample_id = as.character(g), datetime = g)
}

## The fixture reference model with hydro$value = Inf on one day.
##   * "rainfall": the API recursion carries the Inf forward, so every date
##     from 2023-06-01 on has non-finite features (predict errors there).
##   * "stage": hydro features are antecedent means, so only dates whose
##     tau_long window contains 2023-06-01 are non-finite; dates after the
##     window are finite again (model -> static -> model).
.pfr_inf_model <- function(hydro_type = c("rainfall", "stage")) {
  hydro_type <- match.arg(hydro_type)
  rm <- .pf_reference_model()
  rm$hydro$value[as.Date(rm$hydro$date) == .pfr_inf_date] <- Inf
  rm$hydro_type <- hydro_type
  rm
}

.pfr_resolver_inf_cases <- function() {
  tdf <- .pfr_inf_target()
  list(
    rainfall = hydroSense:::.resolve_ref_norm_instant(.pfr_inf_model("rainfall"), tdf),
    stage    = hydroSense:::.resolve_ref_norm_instant(.pfr_inf_model("stage"), tdf)
  )
}

## R2 -------------------------------------------------------------------------

## One point, one analyte, one sample per day.
.pfr_daily_series <- function(v, start = as.Date("2018-01-01")) {
  d <- start + seq_along(v) - 1L
  tibble::tibble(
    sample_id = sprintf("s%05d", seq_along(v)), site_id = "P1", datetime = d,
    analyte = "msPAF", value = v, detected = TRUE
  )
}

## The four reviewer repros.  Each: series, focal dates (every 3 days over the
## record plus a month past its end), tau and window in days.
.pfr_prefix_inputs <- function() {
  focal_for <- function(df) {
    r <- range(df$datetime)
    unique(c(seq(r[1], r[2] + 30L, by = "3 days"), r[2]))
  }
  mk <- function(v, tau, window) {
    df <- .pfr_daily_series(v)
    list(df = df, focal = focal_for(df), tau = tau, window = window)
  }

  # (a) 50 for two years, then 1e-9 for two years.
  a <- mk(c(rep(50, 730), rep(1e-9, 730)), tau = 365, window = 365)
  # (b) 1e3 for two years, then 1e-12 for one year.
  b <- mk(c(rep(1e3, 730), rep(1e-12, 365)), tau = 30, window = 30)
  # (c) U(0, 100) for two years, then 1e-10 * U(0, 1) for two years.
  set.seed(301L)
  c_ <- mk(c(stats::runif(730, 0, 100), 1e-10 * stats::runif(730)),
           tau = 90, window = 365)
  # (d) msPAF-like: a 1e-7 log-normal background, sparse decaying pulses
  #     (peak 0.05-0.9) between days 400 and 700, then two quiet years.
  set.seed(302L)
  n <- 4L * 365L
  v <- 1e-7 * exp(stats::rnorm(n, 0, 0.3))
  days <- 400:700
  hit <- days[stats::runif(length(days)) < 0.06]
  for (q in hit) {
    L <- 10L
    v[q:(q + L - 1L)] <- v[q:(q + L - 1L)] +
      stats::runif(1, 0.05, 0.9) * exp(-(0:(L - 1L)) / 3)
  }
  v <- pmin(v, 1)
  d <- mk(v, tau = 365, window = 365)

  list(a = a, b = b, c = c_, d = d)
}

## time_weighted_aggregate() on each repro, arithmetic AND geometric mean.
.pfr_prefix_cases <- function() {
  inp <- .pfr_prefix_inputs()
  out <- list()
  for (nm in names(inp)) {
    x <- inp[[nm]]
    for (s in c("arith_mean", "geom_mean")) {
      out[[paste0(nm, "_", s)]] <- suppressMessages(time_weighted_aggregate(
        x$df, x$focal, tau = x$tau, tau_units = "d", window = x$window,
        window_units = "d", summary = s
      ))
    }
  }
  out
}
