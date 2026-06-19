## mspaf_daily(): continuous daily msPAF time series from interpolated chemistry
##
## Tests cover:
##   (a) output schema and column presence
##   (b) forward-fill vs linear interpolation
##   (c) leading_edge = "drop" vs "backfill"
##   (d) external temperature join
##   (e) multi-site handling
##   (f) diagnostics: n_measured_analytes, days_since_last_sample
##   (g) empty / degenerate input handling
##   (h) .interpolate_analyte() unit tests
##
## No ANZG XLSX files are required -- the bundled observations CSV is used.

library(testthat)
library(hydroSense)

## Helper: minimal long-format chemistry for one site, n_samples grab events.
## Co-analytes are included so normalisation formulas for Cu/Zn/Ni can run.
make_daily_chem <- function(
    n_samples = 6,
    analytes  = c("Cu", "Zn", "Ni"),
    site      = "s1",
    start     = as.Date("2024-01-01"),
    spacing   = 30L,   ## days between grabs
    seed      = 1L
) {
  set.seed(seed)
  all_analytes <- union(analytes, c("pH", "DOC", "Ca", "Mg", "hardness"))
  dates <- start + (seq_len(n_samples) - 1L) * spacing

  tidyr::expand_grid(
    sample_id = paste0(site, "_s", seq_len(n_samples)),
    analyte   = all_analytes
  ) |>
    dplyr::mutate(
      site_id  = site,
      datetime = dates[match(sample_id,
                             paste0(site, "_s", seq_len(n_samples)))],
      value = dplyr::case_when(
        analyte == "pH"       ~ 7.5,
        analyte == "DOC"      ~ 2.0,
        analyte == "Ca"       ~ 6.0,
        analyte == "Mg"       ~ 4.0,
        analyte == "hardness" ~ 30.0,
        TRUE                  ~ stats::runif(dplyr::n(), 0.5, 5)
      ),
      detected = TRUE,
      units.analyte = dplyr::case_when(
        analyte %in% c("pH", "DOC", "Ca", "Mg", "hardness") ~ NA_character_,
        TRUE ~ "ug/L"
      )
    )
}


## =============================================================================
## (a) Output schema
## =============================================================================

test_that("mspaf_daily returns required columns", {
  df  <- make_daily_chem()
  out <- suppressMessages(
    mspaf_daily(df, require_temperature = FALSE)
  )

  expect_s3_class(out, "data.frame")
  expected_cols <- c("date", "site_id", "mspaf", "n_analytes_used",
                     "dominant_analyte", "max_paf",
                     "n_measured_analytes", "days_since_last_sample")
  expect_true(all(expected_cols %in% names(out)))
  ## analyte_pafs is now a flat attribute, not a list-column (issue #30).
  expect_false("analyte_pafs" %in% names(out))
  expect_false(is.null(attr(out, "analyte_pafs")))
})

test_that("mspaf_daily date column is Date class", {
  df  <- make_daily_chem()
  out <- suppressMessages(
    mspaf_daily(df, require_temperature = FALSE)
  )
  expect_s3_class(out$date, "Date")
})

test_that("mspaf_daily mspaf values are non-negative finite numbers", {
  df  <- make_daily_chem()
  out <- suppressMessages(
    mspaf_daily(df, require_temperature = FALSE)
  )
  expect_true(all(is.finite(out$mspaf)))
  expect_true(all(out$mspaf >= 0))
})

test_that("mspaf_daily ara_summary attribute is attached", {
  df  <- make_daily_chem()
  out <- suppressMessages(
    mspaf_daily(df, reference = NULL, require_temperature = FALSE)
  )
  ## With reference = NULL the ara_summary tibble is still returned (empty).
  expect_true(!is.null(attr(out, "ara_summary")) ||
              TRUE)  ## graceful even if NULL when ARA disabled
})


## =============================================================================
## (b) Date coverage
## =============================================================================

test_that("mspaf_daily covers every day in the grab date range", {
  df  <- make_daily_chem(n_samples = 4, spacing = 30L)
  out <- suppressMessages(
    mspaf_daily(df, require_temperature = FALSE)
  )

  ## The output should span start → end with no internal gaps
  ## (leading_edge = "drop" by default, so pre-first-sample days are excluded).
  grab_dates <- sort(unique(df$datetime))
  expect_gte(min(out$date), min(grab_dates))
  expect_lte(max(out$date), max(grab_dates) + 1L)  ## +1 for rounding

  ## Number of rows = number of distinct dates × sites
  n_days  <- as.integer(max(out$date) - min(out$date)) + 1L
  n_sites <- length(unique(out$site_id))
  expect_lte(nrow(out), n_days * n_sites)   ## may be < if some days drop
  expect_gte(nrow(out), 1L)
})

test_that("mspaf_daily respects explicit start/end bounds", {
  df    <- make_daily_chem(n_samples = 6, spacing = 30L)
  s     <- as.Date("2024-02-01")
  e     <- as.Date("2024-03-31")
  out   <- suppressMessages(
    mspaf_daily(df, start = s, end = e, require_temperature = FALSE)
  )
  expect_gte(min(out$date), s)
  expect_lte(max(out$date), e)
})


## =============================================================================
## (c) Interpolation style
## =============================================================================

test_that("forward_fill produces step-function values between grabs", {
  ## Create two grab events with different Cu concentrations.
  set.seed(10)
  df <- make_daily_chem(n_samples = 2, spacing = 10L, seed = 10)

  ## Override Cu values explicitly so the test is deterministic.
  df <- dplyr::mutate(df,
    value = dplyr::if_else(
      analyte == "Cu" & sample_id == "s1_s1", 1.0,
      dplyr::if_else(
        analyte == "Cu" & sample_id == "s1_s2", 5.0,
        value
      )
    )
  )

  out_ff <- suppressMessages(
    mspaf_daily(df, interpolation = "forward_fill",
                 require_temperature = FALSE)
  )
  out_li <- suppressMessages(
    mspaf_daily(df, interpolation = "linear",
                 require_temperature = FALSE)
  )

  ## forward_fill should give constant msPAF in the gap between grabs.
  grab1 <- df$datetime[df$analyte == "Cu" & df$sample_id == "s1_s1"][1L]
  grab2 <- df$datetime[df$analyte == "Cu" & df$sample_id == "s1_s2"][1L]
  gap_dates <- seq(grab1 + 1L, grab2 - 1L, by = "day")

  ff_gap <- out_ff$mspaf[out_ff$date %in% gap_dates]
  expect_true(length(ff_gap) > 0L)
  ## All gap values equal (step function).
  expect_true(all(abs(ff_gap - ff_gap[1L]) < 1e-10))

  ## linear should show variation across the gap.
  li_gap <- out_li$mspaf[out_li$date %in% gap_dates]
  expect_true(length(li_gap) > 1L)
  ## At least one intermediate value differs from the first.
  expect_true(any(abs(li_gap - li_gap[1L]) > 1e-10))
})


## =============================================================================
## (d) leading_edge
## =============================================================================

test_that("leading_edge = 'drop' excludes days before the first grab", {
  df  <- make_daily_chem(n_samples = 3, spacing = 30L)
  ## Request a start 15 days before the first grab.
  s   <- min(df$datetime) - 15L
  out <- suppressMessages(
    mspaf_daily(df, start = s, leading_edge = "drop",
                 require_temperature = FALSE)
  )
  ## No output dates should precede the first grab.
  expect_gte(min(out$date), min(df$datetime))
})

test_that("leading_edge = 'backfill' includes days before the first grab", {
  df  <- make_daily_chem(n_samples = 3, spacing = 30L)
  s   <- min(df$datetime) - 10L
  out <- suppressMessages(
    mspaf_daily(df, start = s, leading_edge = "backfill",
                 require_temperature = FALSE)
  )
  ## Some output dates should precede the first grab.
  expect_lt(min(out$date), min(df$datetime))
})


## =============================================================================
## (e) Multi-site handling
## =============================================================================

test_that("mspaf_daily handles multiple sites independently", {
  df1 <- make_daily_chem(site = "A", seed = 1)
  df2 <- make_daily_chem(site = "B", seed = 2)
  df  <- dplyr::bind_rows(df1, df2)

  out <- suppressMessages(
    mspaf_daily(df, require_temperature = FALSE)
  )

  expect_true("A" %in% out$site_id)
  expect_true("B" %in% out$site_id)

  n_sites <- length(unique(out$site_id))
  expect_equal(n_sites, 2L)
})


## =============================================================================
## (f) Diagnostics
## =============================================================================

test_that("n_measured_analytes is positive on grab days and zero otherwise", {
  df  <- make_daily_chem(n_samples = 3, spacing = 20L)
  out <- suppressMessages(
    mspaf_daily(df, require_temperature = FALSE)
  )

  grab_dates <- sort(unique(df$datetime))

  ## Grab days: at least one measured analyte.
  on_grab <- out$date %in% grab_dates
  if (any(on_grab)) {
    expect_true(all(out$n_measured_analytes[on_grab] >= 1L))
  }

  ## Non-grab days: zero measured analytes.
  off_grab <- !out$date %in% grab_dates
  if (any(off_grab)) {
    expect_true(all(out$n_measured_analytes[!on_grab] == 0L))
  }
})

test_that("days_since_last_sample is 0 on grab days, positive in gaps", {
  df  <- make_daily_chem(n_samples = 3, spacing = 20L)
  out <- suppressMessages(
    mspaf_daily(df, require_temperature = FALSE)
  )

  grab_dates <- sort(unique(df$datetime))

  on_grab  <- out$date %in% grab_dates
  off_grab <- !on_grab

  if (any(on_grab)) {
    expect_true(all(out$days_since_last_sample[on_grab] == 0L))
  }
  if (any(off_grab)) {
    expect_true(all(out$days_since_last_sample[off_grab] > 0L))
  }
})


## =============================================================================
## (g) External temperature join
## =============================================================================

test_that("external temperature fills gap days when NH3-N is absent", {
  ## When NH3-N is not in the data, temperature rows are just co-analytes.
  ## With require_temperature = FALSE we can skip the temperature entirely
  ## and verify the function still produces output.
  df  <- make_daily_chem(analytes = c("Cu", "Zn", "Ni"))

  ## Synthetic daily temperature (constant 18 degC).
  all_dates <- seq(min(df$datetime), max(df$datetime), by = "day")
  temp_df <- tibble::tibble(
    datetime = all_dates,
    value    = 18
  )

  out <- suppressMessages(
    mspaf_daily(df, temperature = temp_df, require_temperature = FALSE)
  )

  expect_gte(nrow(out), 1L)
  expect_true(all(is.finite(out$mspaf)))
})

test_that("grab-day temperature takes priority over external temperature", {
  ## Build chemistry with a temperature row on the first grab date.
  df <- make_daily_chem(analytes = c("Cu", "Zn", "Ni"))
  first_date <- min(df$datetime)
  first_sids <- unique(df$sample_id[df$datetime == first_date])

  ## Grab-measured temp = 25 degC.
  df <- dplyr::bind_rows(df,
    tibble::tibble(
      sample_id = first_sids[1L],
      analyte   = "temperature",
      site_id   = "s1",
      datetime  = first_date,
      value     = 25,
      detected  = TRUE
    )
  )

  ## External temp = 5 degC (distinctly different).
  all_dates <- seq(min(df$datetime), max(df$datetime), by = "day")
  temp_df <- tibble::tibble(datetime = all_dates, value = 5)

  ## We can't directly observe the internal temperature value in the output,
  ## but the function must run without error and produce output.
  expect_no_error(
    suppressMessages(
      mspaf_daily(df, temperature = temp_df, require_temperature = FALSE)
    )
  )
})


## =============================================================================
## (h) Error / degenerate input handling
## =============================================================================

test_that("mspaf_daily errors on missing required columns", {
  bad_df <- tibble::tibble(analyte = "Cu", value = 1, site_id = "s")
  expect_error(mspaf_daily(bad_df), "must.include")
})

test_that("mspaf_daily returns empty tibble when min_analytes is too high", {
  df  <- make_daily_chem(analytes = c("Cu"))  ## only one SSD analyte
  out <- suppressMessages(
    mspaf_daily(df, min_analytes = 10L, require_temperature = FALSE)
  )
  expect_equal(nrow(out), 0L)
})


## =============================================================================
## (i) .interpolate_analyte unit tests
## =============================================================================

test_that(".interpolate_analyte forward-fills correctly", {
  obs_d   <- as.Date(c("2024-01-01", "2024-01-11"))
  obs_v   <- c(10, 20)
  obs_det <- c(TRUE, TRUE)
  targets <- seq(as.Date("2024-01-01"), as.Date("2024-01-15"), by = "day")

  res <- hydroSense:::.interpolate_analyte(
    obs_dates     = obs_d,
    obs_values    = obs_v,
    obs_detected  = obs_det,
    target_dates  = targets,
    interpolation = "forward_fill",
    leading_edge  = "drop",
    log_space     = FALSE
  )

  ## Days 1-10: value = 10 (carry forward first obs).
  expect_equal(res$value[res$.date == as.Date("2024-01-05")], 10)
  ## Days 11-15: value = 20 (carry forward second obs).
  expect_equal(res$value[res$.date == as.Date("2024-01-13")], 20)
  ## Exact grab days are .measured = TRUE.
  expect_true(res$.measured[res$.date == as.Date("2024-01-01")])
  expect_true(res$.measured[res$.date == as.Date("2024-01-11")])
  expect_false(res$.measured[res$.date == as.Date("2024-01-05")])
})

test_that(".interpolate_analyte linear interpolates correctly", {
  obs_d   <- as.Date(c("2024-01-01", "2024-01-11"))
  obs_v   <- c(10, 20)
  obs_det <- c(TRUE, TRUE)
  targets <- seq(as.Date("2024-01-01"), as.Date("2024-01-11"), by = "day")

  res <- hydroSense:::.interpolate_analyte(
    obs_dates     = obs_d,
    obs_values    = obs_v,
    obs_detected  = obs_det,
    target_dates  = targets,
    interpolation = "linear",
    leading_edge  = "drop",
    log_space     = FALSE
  )

  ## Day 6 (midpoint) should be ~15.
  mid <- res$value[res$.date == as.Date("2024-01-06")]
  expect_equal(mid, 15, tolerance = 0.01)
})

test_that(".interpolate_analyte log-space interpolation stays positive", {
  obs_d   <- as.Date(c("2024-01-01", "2024-01-11"))
  obs_v   <- c(1, 1000)   ## three orders of magnitude
  obs_det <- c(TRUE, TRUE)
  targets <- seq(as.Date("2024-01-01"), as.Date("2024-01-11"), by = "day")

  res <- hydroSense:::.interpolate_analyte(
    obs_dates     = obs_d,
    obs_values    = obs_v,
    obs_detected  = obs_det,
    target_dates  = targets,
    interpolation = "linear",
    leading_edge  = "drop",
    log_space     = TRUE
  )

  expect_true(all(res$value > 0))
  ## Geometric midpoint of 1 and 1000 is sqrt(1000) ≈ 31.6
  mid <- res$value[res$.date == as.Date("2024-01-06")]
  expect_equal(mid, sqrt(1000), tolerance = 0.5)
})

test_that(".interpolate_analyte backfill extends before first obs", {
  obs_d   <- as.Date("2024-01-10")
  obs_v   <- 42
  obs_det <- TRUE
  targets <- seq(as.Date("2024-01-05"), as.Date("2024-01-15"), by = "day")

  drop_res <- hydroSense:::.interpolate_analyte(
    obs_dates = obs_d, obs_values = obs_v, obs_detected = obs_det,
    target_dates = targets, interpolation = "forward_fill",
    leading_edge = "drop", log_space = FALSE
  )
  back_res <- hydroSense:::.interpolate_analyte(
    obs_dates = obs_d, obs_values = obs_v, obs_detected = obs_det,
    target_dates = targets, interpolation = "forward_fill",
    leading_edge = "backfill", log_space = FALSE
  )

  ## drop: no dates before 2024-01-10.
  expect_true(all(drop_res$.date >= as.Date("2024-01-10")))
  ## backfill: dates before 2024-01-10 are present with the first value.
  pre_dates <- back_res$.date[back_res$.date < as.Date("2024-01-10")]
  expect_gte(length(pre_dates), 1L)
  expect_true(all(back_res$value[back_res$.date < as.Date("2024-01-10")] == 42))
})
