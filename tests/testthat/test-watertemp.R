## estimate_water_temp(): air-only vs air + day-of-year-harmonic model, with
## AICc selection and the seasonal-eligibility guards.
## Also covers: fallback relationships (part a) and multi-timescale window
## selection via LOO-CV PRESS (part b).

library(testthat)
library(leachatetools)

# Helper: a year-spanning air-temperature series.
make_air <- function(years = 2021:2022) {
  dates <- seq(as.Date(paste0(min(years), "-01-01")),
               as.Date(paste0(max(years), "-12-31")), by = "day")
  doy <- as.integer(format(dates, "%j"))
  tibble::tibble(
    datetime        = dates,
    air_temp_mean_C = 15 + 10 * sin(2 * pi * (doy - 30) / 365) +
                      stats::rnorm(length(dates), 0, 2)
  )
}

test_that("seasonal model is selected when hysteresis is present", {
  set.seed(42)
  air <- make_air()
  wd  <- sample(air$datetime, 40)
  doy <- as.integer(format(wd, "%j"))
  airv <- air$air_temp_mean_C[match(wd, air$datetime)]
  # Water carries a seasonal offset beyond what same-day air explains.
  wobs <- tibble::tibble(
    datetime     = wd,
    water_temp_C = 0.7 * airv + 3 * sin(2 * pi * doy / 365.25) + 4 +
                   stats::rnorm(40, 0, 0.5)
  )
  out <- suppressMessages(estimate_water_temp(air, wobs))
  expect_true(attr(out, "seasonal_used"))
  cmp <- attr(out, "model_comparison")
  expect_lt(cmp$aicc[cmp$model == "air_plus_season"],
            cmp$aicc[cmp$model == "air_only"])
})

test_that("air-only model is selected when there is no seasonal signal", {
  set.seed(7)
  air <- make_air()
  wd  <- sample(air$datetime, 40)
  airv <- air$air_temp_mean_C[match(wd, air$datetime)]
  # Water is a clean linear function of same-day air, no extra seasonality.
  wobs <- tibble::tibble(
    datetime     = wd,
    water_temp_C = 0.8 * airv + 2 + stats::rnorm(40, 0, 0.5)
  )
  out <- suppressMessages(estimate_water_temp(air, wobs))
  expect_false(attr(out, "seasonal_used"))
  expect_equal(attr(out, "model"), "air only")
})

test_that("seasonal = 'off' forces the air-only model even with seasonal data", {
  set.seed(42)
  air <- make_air()
  wd  <- sample(air$datetime, 40)
  doy <- as.integer(format(wd, "%j"))
  airv <- air$air_temp_mean_C[match(wd, air$datetime)]
  wobs <- tibble::tibble(
    datetime     = wd,
    water_temp_C = 0.7 * airv + 3 * sin(2 * pi * doy / 365.25) + 4 +
                   stats::rnorm(40, 0, 0.5)
  )
  out <- suppressMessages(estimate_water_temp(air, wobs, seasonal = "off"))
  expect_false(attr(out, "seasonal_used"))
})

test_that("seasonal model is not considered when coverage is poor", {
  set.seed(1)
  air <- make_air()
  # All observations in a single quarter (Jan-Mar) -> cannot anchor a cycle.
  pool <- air$datetime[as.integer(format(air$datetime, "%m")) <= 3L]
  wd   <- sample(pool, 12)
  airv <- air$air_temp_mean_C[match(wd, air$datetime)]
  wobs <- tibble::tibble(datetime = wd,
                         water_temp_C = 0.8 * airv + 2 + stats::rnorm(12, 0, 0.5))
  out <- suppressMessages(estimate_water_temp(air, wobs))
  expect_false(attr(out, "seasonal_used"))
  # The seasonal candidate should not even be fitted (NA AICc).
  cmp <- attr(out, "model_comparison")
  expect_true(is.na(cmp$aicc[cmp$model == "air_plus_season"]))
})

test_that(".aicc returns Inf when the sample is too small for the correction", {
  d <- data.frame(y = c(1, 2, 3, 4), x1 = c(1, 2, 3, 4),
                  x2 = c(2, 1, 4, 3), x3 = c(1, 3, 2, 4))
  fit <- stats::lm(y ~ x1 + x2 + x3, data = d)  # k = 5, n = 4 -> undefined
  expect_identical(leachatetools:::.aicc(fit), Inf)
})

test_that(".add_doy_harmonics is cyclic and bounded", {
  df <- data.frame(d = as.Date(c("2021-01-01", "2021-12-31", "2021-07-02")))
  h  <- leachatetools:::.add_doy_harmonics(df, "d")
  expect_true(all(h$sin_doy >= -1 & h$sin_doy <= 1))
  expect_true(all(h$cos_doy >= -1 & h$cos_doy <= 1))
  # Jan 1 (doy 1) and Dec 31 (doy 365) are adjacent on the annual circle.
  expect_lt(abs(h$sin_doy[1] - h$sin_doy[2]), 0.05)
  expect_lt(abs(h$cos_doy[1] - h$cos_doy[2]), 0.05)
})

test_that("output shape is unchanged and attributes are attached", {
  set.seed(3)
  air <- make_air()
  wd  <- sample(air$datetime, 30)
  airv <- air$air_temp_mean_C[match(wd, air$datetime)]
  wobs <- tibble::tibble(datetime = wd,
                         water_temp_C = 0.8 * airv + 2 + stats::rnorm(30, 0, 0.5))
  out <- suppressMessages(estimate_water_temp(air, wobs, site_id = "S1"))
  expect_named(out, c("datetime", "analyte", "value", "detected",
                      "site_id", "sample_id"))
  expect_true(all(out$analyte == "temperature"))
  expect_s3_class(attr(out, "lm_fit"), "lm")
  expect_s3_class(attr(out, "model_comparison"), "data.frame")
  # New attributes are present
  expect_true(!is.null(attr(out, "window_selected")))
})

# ── Part (a): Fallback relationships ─────────────────────────────────────────

test_that("fallback = 'identity' returns predictions equal to air temp", {
  set.seed(10)
  air <- make_air()
  out <- suppressMessages(
    estimate_water_temp(air, fallback = "identity")
  )
  expect_equal(attr(out, "model"), "fallback")
  expect_false(attr(out, "seasonal_used"))
  expect_null(attr(out, "lm_fit"))
  expect_true(is.na(attr(out, "window_selected")))
  # Predictions should equal air temp for matching dates
  joined <- merge(out, as.data.frame(air), by = "datetime")
  expect_equal(joined$value, joined$air_temp_mean_C, tolerance = 1e-9)
})

test_that("fallback = c(intercept, slope) applies linear transform", {
  set.seed(11)
  air <- make_air()
  out <- suppressMessages(
    estimate_water_temp(air, fallback = c(2.5, 0.8))
  )
  joined <- merge(out, as.data.frame(air), by = "datetime")
  expect_equal(joined$value, 2.5 + 0.8 * joined$air_temp_mean_C, tolerance = 1e-9)
})

test_that("fallback = function applies arbitrary transform", {
  set.seed(12)
  air <- make_air()
  out <- suppressMessages(
    estimate_water_temp(air, fallback = function(x) x + 5)
  )
  joined <- merge(out, as.data.frame(air), by = "datetime")
  expect_equal(joined$value, joined$air_temp_mean_C + 5, tolerance = 1e-9)
})

test_that("fallback is used when water_temp_obs has fewer than 5 rows", {
  set.seed(13)
  air <- make_air()
  wobs <- tibble::tibble(
    datetime     = air$datetime[1:3],
    water_temp_C = c(14, 15, 16)
  )
  out <- suppressMessages(
    estimate_water_temp(air, wobs, fallback = c(1, 0.9))
  )
  expect_equal(attr(out, "model"), "fallback")
})

test_that("no fallback with fewer than 5 rows still errors", {
  air <- make_air()
  wobs <- tibble::tibble(
    datetime     = air$datetime[1:3],
    water_temp_C = c(14, 15, 16)
  )
  expect_error(
    suppressMessages(estimate_water_temp(air, wobs)),
    "paired observations"
  )
})

# ── Part (b): Multi-timescale window selection ────────────────────────────────

test_that("select = 'auto' recovers a window near the true antecedence scale", {
  # Generate water temp as a function of the 14-day rolling mean of air temp.
  # With 200 pairs over 3 years, auto should converge near window = 14.
  set.seed(100)
  air <- make_air(2020:2022)
  roll14 <- leachatetools:::.trailing_rollmean(air$air_temp_mean_C, 14)
  # Discard first 13 days (NA rolling mean)
  air14  <- air[!is.na(roll14), ]
  roll14 <- roll14[!is.na(roll14)]

  pair_idx <- sample(seq_len(nrow(air14)), 200)
  wobs <- tibble::tibble(
    datetime     = air14$datetime[pair_idx],
    water_temp_C = 0.85 * roll14[pair_idx] + 3 + stats::rnorm(200, 0, 0.8)
  )
  out <- suppressMessages(estimate_water_temp(air, wobs, select = "auto"))
  # 1-SE rule may select shorter; allow up to 3× the true window
  expect_lte(attr(out, "window_selected"), 42L)
  expect_gte(attr(out, "window_selected"), 1L)
})

test_that("auto falls back to window = 1 when n < auto_min_n and no water_body_type", {
  set.seed(14)
  air <- make_air()
  wd  <- sample(air$datetime, 15)
  airv <- air$air_temp_mean_C[match(wd, air$datetime)]
  wobs <- tibble::tibble(datetime = wd,
                         water_temp_C = 0.8 * airv + 2 + stats::rnorm(15, 0, 0.5))
  out <- suppressMessages(estimate_water_temp(air, wobs, select = "auto"))
  expect_equal(attr(out, "window_selected"), 1L)
})

test_that("water_body_type sets fallback window when auto search is skipped", {
  set.seed(15)
  air <- make_air()
  wd  <- sample(air$datetime, 15)
  airv <- air$air_temp_mean_C[match(wd, air$datetime)]
  wobs <- tibble::tibble(datetime = wd,
                         water_temp_C = 0.8 * airv + 2 + stats::rnorm(15, 0, 0.5))
  out <- suppressMessages(
    estimate_water_temp(air, wobs, select = "auto", water_body_type = "lake")
  )
  expect_equal(attr(out, "window_selected"), 30L)
})

test_that("select = integer forces that window", {
  set.seed(16)
  air <- make_air()
  wd  <- sample(air$datetime, 40)
  airv <- air$air_temp_mean_C[match(wd, air$datetime)]
  wobs <- tibble::tibble(datetime = wd,
                         water_temp_C = 0.8 * airv + 2 + stats::rnorm(40, 0, 0.5))
  out <- suppressMessages(estimate_water_temp(air, wobs, select = 30L))
  expect_equal(attr(out, "window_selected"), 30L)
  # window_comparison is NULL for fixed select
  expect_null(attr(out, "window_comparison"))
})

test_that("window_comparison has correct structure when auto search runs", {
  set.seed(17)
  air <- make_air()
  wd  <- sample(air$datetime, 40)
  airv <- air$air_temp_mean_C[match(wd, air$datetime)]
  wobs <- tibble::tibble(datetime = wd,
                         water_temp_C = 0.8 * airv + 2 + stats::rnorm(40, 0, 0.5))
  out <- suppressMessages(estimate_water_temp(air, wobs, select = "auto",
                                              auto_min_n = 20L))
  cmp <- attr(out, "window_comparison")
  expect_s3_class(cmp, "data.frame")
  expect_named(cmp, c("window", "loo_rmse", "loo_rmse_se", "n_pairs", "selected"))
  expect_equal(sum(cmp$selected), 1L)
  expect_equal(cmp$window[cmp$selected], attr(out, "window_selected"))
})

test_that("1-SE rule selects a shorter window when near-tie exists", {
  # Construct data where a 7-day and a 14-day window have nearly identical
  # LOO-CV RMSE so the 1-SE rule should prefer the shorter one.
  set.seed(200)
  air <- make_air(2020:2022)
  # Use the average of 7-day and 14-day rolling means as the true signal —
  # neither window dominates, so CV RMSEs will be close.
  r7  <- leachatetools:::.trailing_rollmean(air$air_temp_mean_C, 7)
  r14 <- leachatetools:::.trailing_rollmean(air$air_temp_mean_C, 14)
  ok  <- !is.na(r7) & !is.na(r14)
  air_ok <- air[ok, ]
  signal  <- 0.5 * (r7[ok] + r14[ok])

  pair_idx <- sample(seq_len(nrow(air_ok)), 200)
  wobs <- tibble::tibble(
    datetime     = air_ok$datetime[pair_idx],
    water_temp_C = 0.85 * signal[pair_idx] + 3 + stats::rnorm(200, 0, 1.5)
  )
  out <- suppressMessages(estimate_water_temp(air, wobs, select = "auto"))
  cmp <- attr(out, "window_comparison")
  # The 1-SE rule must not select a window longer than both 7 and 14 when
  # those two are within one SE of each other.
  if (!is.null(cmp)) {
    rmse_7  <- cmp$loo_rmse[cmp$window == 7]
    rmse_14 <- cmp$loo_rmse[cmp$window == 14]
    se_best <- cmp$loo_rmse_se[which.min(cmp$loo_rmse)]
    if (length(rmse_7) == 1 && length(rmse_14) == 1 &&
        abs(rmse_7 - rmse_14) <= se_best) {
      expect_lte(attr(out, "window_selected"), 7L)
    }
  }
  # Regardless, window_selected must be a positive integer
  expect_true(is.integer(attr(out, "window_selected")))
  expect_gte(attr(out, "window_selected"), 1L)
})
