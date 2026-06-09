## Tests for amspaf_daily() draws mode (issue #16, Chunk E).
##
## Properties tested:
##   F1. ndraws without interpolation = "model" aborts with informative error
##   F2. return = "draws" yields ndraws rows per (date, site)
##   F3. return = "summary" adds amspaf_lower / amspaf_upper columns
##   F4. seed makes draws reproducible
##   F5. Point mode (ndraws = NULL) unchanged — no regression
##   F6. Summary ordering: amspaf_lower <= amspaf <= amspaf_upper

library(testthat)
library(leachatetools)


## ── Shared setup ─────────────────────────────────────────────────────────────

make_chem_f <- function(site, dates, mult = 1, seed = 1L) {
  set.seed(seed)
  analytes <- c("Cu", "Zn", "Ni", "pH", "DOC", "hardness", "Ca", "Mg")
  purrr::map_dfr(dates, function(d) tibble::tibble(
    sample_id     = paste0(site, format(d, "%Y%m%d")),
    site_id       = site,
    datetime      = d,
    analyte       = analytes,
    value         = c(
      exp(stats::rnorm(1, log(0.5), 0.3)) * mult,
      exp(stats::rnorm(1, log(5),   0.4)) * mult,
      exp(stats::rnorm(1, log(0.3), 0.3)) * mult,
      stats::runif(1, 6.5, 8), stats::runif(1, 1, 5), stats::runif(1, 20, 60),
      stats::runif(1, 4, 12),  stats::runif(1, 2, 8)
    ),
    detected      = TRUE,
    units.analyte = dplyr::case_when(
      analyte %in% c("Cu", "Zn", "Ni") ~ "ug/L",
      TRUE                              ~ NA_character_
    )
  ))
}

make_hydro_f <- function(n = 700L, seed = 99L) {
  set.seed(seed)
  tibble::tibble(
    date  = seq(as.Date("2020-07-01"), by = "day", length.out = n),
    value = pmax(0, stats::rnorm(n, 2, 4))
  )
}

## One-time expensive setup: reference model fitted on synthetic data.
.tf <- local({
  ## Biweekly grabs over ~6 months
  dates <- seq(as.Date("2021-01-01"), by = "2 weeks", length.out = 14L)
  hydro <- make_hydro_f()
  ref   <- make_chem_f("reference", dates, seed = 1L)
  tgt   <- make_chem_f("target",    dates, mult = 5, seed = 2L)
  rm    <- fit_reference_model(ref, hydro = hydro, conc_units = "ug/L",
                                min_obs_model = 10L,
                                api_windows_short = 7L,
                                api_windows_long  = 30L)
  list(rm = rm, tgt = tgt, dates = dates)
})


## ── F1: ndraws without interpolation = "model" errors ────────────────────────

test_that("F1: ndraws with forward_fill interpolation aborts", {
  df <- make_chem_f("target", seq(as.Date("2021-01-01"), by = "2 weeks",
                                   length.out = 10L), seed = 2L)
  expect_error(
    amspaf_daily(df, ndraws = 3L, interpolation = "forward_fill",
                 require_temperature = FALSE),
    regexp = "ndraws.*model|model.*ndraws",
    perl   = TRUE
  )
})


## ── F2: return = "draws" → ndraws rows per (date, site) ──────────────────────

test_that("F2: draws output has ndraws rows per unique (date, site_id)", {
  skip_if(is.null(.tf$rm), "Reference model not fitted")
  n  <- 5L
  out <- suppressMessages(
    amspaf_daily(
      .tf$tgt,
      reference_model = .tf$rm,
      interpolation   = "model",
      ndraws          = n,
      seed            = 1L,
      return          = "draws",
      require_temperature = FALSE,
      conc_units      = "ug/L"
    )
  )
  expect_true("draw_id" %in% names(out))
  ## Each (date, site_id) combination should have exactly ndraws rows
  counts <- out |>
    dplyr::group_by(.data$date, .data$site_id) |>
    dplyr::summarise(n_rows = dplyr::n(), .groups = "drop")
  expect_true(all(counts$n_rows == n),
              label = "every (date, site_id) has exactly ndraws rows")
  ## draw_id values should be 1..n
  expect_setequal(unique(out$draw_id), seq_len(n))
})


## ── F3: return = "summary" adds CI columns ───────────────────────────────────

test_that("F3: summary output has amspaf_lower and amspaf_upper columns", {
  skip_if(is.null(.tf$rm), "Reference model not fitted")
  out <- suppressMessages(
    amspaf_daily(
      .tf$tgt,
      reference_model = .tf$rm,
      interpolation   = "model",
      ndraws          = 5L,
      seed            = 2L,
      return          = "summary",
      interval        = 0.8,
      require_temperature = FALSE,
      conc_units      = "ug/L"
    )
  )
  expect_true(all(c("amspaf", "amspaf_lower", "amspaf_upper") %in% names(out)))
  ## Standard point-mode columns still present
  expect_true(all(c("date", "site_id", "n_analytes_used",
                    "n_measured_analytes", "days_since_last_sample") %in% names(out)))
  ## No draw_id column in summary output
  expect_false("draw_id" %in% names(out))
  ## One row per (date, site_id)
  expect_equal(nrow(out),
               nrow(dplyr::distinct(dplyr::select(out, "date", "site_id"))))
})


## ── F4: seed makes draws reproducible ────────────────────────────────────────

test_that("F4: same seed gives identical draws output", {
  skip_if(is.null(.tf$rm), "Reference model not fitted")
  run <- function() suppressMessages(
    amspaf_daily(
      .tf$tgt,
      reference_model = .tf$rm,
      interpolation   = "model",
      ndraws          = 4L,
      seed            = 42L,
      return          = "draws",
      require_temperature = FALSE,
      conc_units      = "ug/L"
    )
  )
  out1 <- run()
  out2 <- run()
  expect_equal(out1$amspaf, out2$amspaf,
               label = "same seed produces identical amspaf draws")
})


## ── F5: point mode regression — standard schema returned ─────────────────────

test_that("F5: ndraws = NULL returns standard point-mode schema", {
  skip_if(is.null(.tf$rm), "Reference model not fitted")
  out <- suppressMessages(
    amspaf_daily(
      .tf$tgt,
      reference_model = .tf$rm,
      interpolation   = "model",
      ndraws          = NULL,
      require_temperature = FALSE,
      conc_units      = "ug/L"
    )
  )
  expected_cols <- c("date", "site_id", "amspaf", "n_analytes_used",
                     "dominant_analyte", "max_paf",
                     "n_measured_analytes", "days_since_last_sample",
                     "analyte_pafs")
  expect_true(all(expected_cols %in% names(out)))
  expect_false("amspaf_lower" %in% names(out))
  expect_false("draw_id"      %in% names(out))
  expect_s3_class(out$date, "Date")
  expect_true(all(is.finite(out$amspaf)))
})


## ── F6: ordering amspaf_lower <= amspaf <= amspaf_upper ──────────────────────

test_that("F6: summary CI bounds straddle the central estimate", {
  skip_if(is.null(.tf$rm), "Reference model not fitted")
  out <- suppressMessages(
    amspaf_daily(
      .tf$tgt,
      reference_model = .tf$rm,
      interpolation   = "model",
      ndraws          = 20L,
      seed            = 7L,
      return          = "summary",
      interval        = 0.9,
      central         = "median",
      require_temperature = FALSE,
      conc_units      = "ug/L"
    )
  )
  expect_true(all(out$amspaf_lower <= out$amspaf + 1e-9),
              label = "lower bound <= central")
  expect_true(all(out$amspaf <= out$amspaf_upper + 1e-9),
              label = "central <= upper bound")
})
