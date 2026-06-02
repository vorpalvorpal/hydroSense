## estimate_water_temp(): air-only vs air + day-of-year-harmonic model, with
## AICc selection and the seasonal-eligibility guards.

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
})
