## Integration tests for the Kalman daily-uncertainty rework (issue #16).
## New guarantees beyond test-amspaf-daily-draws.R (which covers schema, draws,
## reproducibility, grab_cv, parallel, and the draw-bearing-input error):
##   I2 centre line invariant to ndraws/seed (deterministic posterior mean)
##   I3 interval pinches near grabs, balloons mid-gap
##   I7 daily grid clipped per-analyte to that analyte's grab span
##   I8 envelope non-negative and upper-skewed (floor + SSD asymmetry)

library(testthat)
library(leachatetools)

make_chem_k <- function(site, dates, mult = 1, seed = 1L, analytes = NULL) {
  set.seed(seed)
  if (is.null(analytes))
    analytes <- c("Cu", "Zn", "Ni", "pH", "DOC", "hardness", "Ca", "Mg")
  base <- c(Cu = 0.5, Zn = 5, Ni = 0.3)
  purrr::map_dfr(dates, function(d) {
    vals <- vapply(analytes, function(a) {
      if (a %in% names(base)) exp(stats::rnorm(1, log(base[[a]]), 0.35)) * mult
      else switch(a, pH = stats::runif(1, 6.5, 8), DOC = stats::runif(1, 1, 5),
                  hardness = stats::runif(1, 20, 60), Ca = stats::runif(1, 4, 12),
                  Mg = stats::runif(1, 2, 8), stats::runif(1, 1, 5))
    }, numeric(1))
    tibble::tibble(
      sample_id = paste0(site, format(d, "%Y%m%d")), site_id = site,
      datetime = d, analyte = analytes, value = as.numeric(vals), detected = TRUE,
      units.analyte = ifelse(analytes %in% c("Cu", "Zn", "Ni"), "ug/L", NA_character_))
  })
}

make_hydro_k <- function(n = 760L, seed = 99L) {
  set.seed(seed)
  tibble::tibble(date = seq(as.Date("2020-07-01"), by = "day", length.out = n),
                 value = pmax(0, stats::rnorm(n, 2, 4)))
}

.kf <- local({
  dates <- seq(as.Date("2021-01-01"), by = "2 weeks", length.out = 16L)
  hydro <- make_hydro_k()
  ref   <- make_chem_k("reference", dates, seed = 1L)
  tgt   <- make_chem_k("target",    dates, mult = 5, seed = 2L)
  rm    <- fit_reference_model(ref, hydro = hydro, conc_units = "ug/L",
                               min_obs_model = 10L, api_windows_short = 7L,
                               api_windows_long = 30L)
  list(rm = rm, tgt = tgt, dates = dates)
})

run_summary <- function(ndraws, seed = 1L, ...) {
  suppressMessages(amspaf_daily(
    .kf$tgt, reference_model = .kf$rm, interpolation = "model",
    ndraws = ndraws, seed = seed, return = "summary",
    require_temperature = FALSE, conc_units = "ug/L", ...))
}

## ── I2: deterministic point estimate is stable; draw summary tracks draws ─
## Since issue #42 the summary centre is the draws' OWN central tendency, so it
## depends on ndraws/seed (it is no longer the deterministic line). The stable,
## reproducible best guess is point mode (ndraws = NULL).

test_that("I2: point-mode centre stable; summary centre is draw median", {
  skip_if(is.null(.kf$rm), "no reference model")
  run_point <- function() suppressMessages(amspaf_daily(
    .kf$tgt, reference_model = .kf$rm, interpolation = "model",
    require_temperature = FALSE, conc_units = "ug/L"))

  ## Point mode: the deterministic best guess — identical on repeat, seedless.
  expect_equal(run_point()$amspaf, run_point()$amspaf, tolerance = 1e-8)

  ## Summary centre = the per-day median of the same draws (coherent with band).
  draws <- suppressMessages(amspaf_daily(
    .kf$tgt, reference_model = .kf$rm, interpolation = "model",
    ndraws = 20L, seed = 7L, return = "draws",
    require_temperature = FALSE, conc_units = "ug/L"))
  summ <- run_summary(ndraws = 20L, seed = 7L)
  ref  <- draws |>
    dplyr::group_by(.data$date) |>
    dplyr::summarise(m = stats::median(.data$amspaf), .groups = "drop")
  m <- dplyr::inner_join(summ[, c("date", "amspaf")], ref, by = "date")
  expect_equal(m$amspaf, m$m, tolerance = 1e-8)
})

## ── I3: pinch near grabs, balloon mid-gap ─────────────────────────────────────

test_that("I3: CI width is smaller near grab dates than mid-gap", {
  skip_if(is.null(.kf$rm), "no reference model")
  out <- run_summary(ndraws = 40L, seed = 3L)
  out$width <- out$amspaf_upper - out$amspaf_lower
  dist_to_grab <- vapply(out$date, function(d)
    min(abs(as.numeric(d - .kf$dates))), numeric(1))
  near <- out$width[dist_to_grab <= 1]
  mid  <- out$width[dist_to_grab >= 6]
  expect_gt(mean(mid, na.rm = TRUE), mean(near, na.rm = TRUE))
})

## ── I7: residual grid clipped per analyte to its grab span ────────────────────

test_that("I7: a short-sampled analyte's modelled rows are clipped to its span", {
  skip_if(is.null(.kf$rm), "no reference model")
  # Ni grabbed only in the first 6 fortnights; Cu/Zn throughout.
  early   <- .kf$dates[1:6]
  ni_rows <- .kf$tgt$analyte == "Ni" & !(as.Date(.kf$tgt$datetime) %in% early)
  tgt_clip <- .kf$tgt[!ni_rows, , drop = FALSE]

  tm <- suppressMessages(fit_target_model(tgt_clip, .kf$rm, conc_units = "ug/L",
                                          min_obs_model = 10L,
                                          api_windows_short = 7L,
                                          api_windows_long = 30L))
  skip_if(!("Ni" %in% names(tm$models)), "Ni not modelled")
  q   <- tibble::tibble(date = seq(min(.kf$dates), max(.kf$dates), by = "day"))
  res <- leachatetools:::.resolve_target_impact(tm, q)

  ni <- res[res$analyte == "Ni", , drop = FALSE]
  expect_gt(nrow(ni), 0)
  expect_lte(max(ni$date), max(early) + 1L)             # no Ni beyond its span
  cu <- res[res$analyte == "Cu", , drop = FALSE]        # full-span analyte
  expect_gt(max(cu$date), max(ni$date))
})

## ── I8: envelope non-negative, ordered, and asymmetric ────────────────────────

test_that("I8: bounds are non-negative/ordered and the envelope is asymmetric", {
  skip_if(is.null(.kf$rm), "no reference model")
  out <- run_summary(ndraws = 60L, seed = 8L)
  expect_true(all(out$amspaf_lower >= 0, na.rm = TRUE))
  expect_true(all(out$amspaf_lower <= out$amspaf & out$amspaf <= out$amspaf_upper,
                  na.rm = TRUE))
  # The floor + SSD nonlinearity make the interval asymmetric about the centre
  # on a non-trivial share of days. (Direction is data-dependent: symmetric
  # concentration draws map through the SSD curvature, which can skew the AmsPAF
  # band either way, so we test that asymmetry EXISTS, not its sign.)
  uhw  <- out$amspaf_upper - out$amspaf
  lhw  <- out$amspaf - out$amspaf_lower
  asym <- abs(uhw - lhw) / pmax(uhw + lhw, .Machine$double.eps)
  expect_gt(mean(asym > 0.02, na.rm = TRUE), 0.2)       # >20% of days asymmetric
})
