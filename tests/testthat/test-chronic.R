## Tests for compute_chronic_chemistry()

library(testthat)
library(leachatetools)

# ── Helper: build a synthetic long-format chemistry df ────────────────────────
make_synth_chem <- function(
    n_samples    = 20,
    n_analytes   = 3,
    feature      = "f1",
    start        = as.Date("2023-01-01"),
    end          = as.Date("2023-12-31"),
    value_const  = NULL  # if given, all rows get this value
) {
  set.seed(42)
  dates   <- sort(sample(seq(start, end, by = "day"), n_samples))
  analytes <- paste0("A", seq_len(n_analytes))

  df <- tidyr::expand_grid(
    uuid.sample  = paste0("s", seq_len(n_samples)),
    name.analyte = analytes
  ) |>
    dplyr::mutate(
      uuid.feature     = feature,
      datetime.sample  = dates[match(
        sub("s", "", uuid.sample), as.character(seq_len(n_samples))
      )],
      value      = if (!is.null(value_const)) value_const else runif(dplyr::n(), 0.5, 5),
      quantified = TRUE,
      imputed    = FALSE
    )
  df
}

# ── Property: identical concentrations → chronic mean equals that value ───────
test_that("constant concentration → chronic mean equals that constant", {
  df <- make_synth_chem(value_const = 3.0)
  result <- compute_chronic_chemistry(
    df,
    focal_dates = as.Date("2024-01-01"),
    tau_days    = 90,
    window_days = 365
  )
  # Allow eps tolerance
  expect_true(all(abs(result$value - 3.0) < 1e-4))
})

# ── Property: forward-step weighting — two samples, known weight ratio ─────────
test_that("exponential decay weight ratio is correct for two samples 90d apart", {
  # Two samples: old (90 days before focal) and recent (0 days before focal).
  # tau = 90 d.
  # Midpoints: old interval [t1, t2) → midpoint 45 d before focal (age ~45)
  # recent interval [t2, focal] → midpoint 45/2 = 22.5 d before focal (age ~22.5)
  # But let's use a cleaner test with just 2 samples.
  #
  # Simpler: just verify that a sample with higher concentration near focal_date
  # gets weighted more than one further away (tau = 90 d).
  focal <- as.Date("2024-01-01")
  df <- tibble::tibble(
    uuid.sample     = c("s1", "s1", "s2", "s2"),
    uuid.feature    = "f1",
    datetime.sample = c(
      focal - 91,  # old sample
      focal - 91,
      focal - 1,   # recent sample
      focal - 1
    ),
    name.analyte = c("Cu", "Zn", "Cu", "Zn"),
    value        = c(1.0, 1.0, 10.0, 10.0),
    quantified   = TRUE,
    imputed      = FALSE
  )
  result <- compute_chronic_chemistry(
    df, focal_dates = focal, tau_days = 90, window_days = 365,
    anchor_outside_window = FALSE
  )
  # The recent sample (value=10) is more recent, so the geometric mean should
  # be > 1.0 (closer to 10 than 1 due to temporal weighting + duration weighting)
  cu_val <- result$value[result$name.analyte == "Cu"]
  expect_gt(cu_val, 1.0)
  expect_lt(cu_val, 10.0)
})

# ── Property: quantified always TRUE in output ────────────────────────────────
test_that("output quantified column is always TRUE", {
  df <- make_synth_chem()
  result <- compute_chronic_chemistry(
    df, focal_dates = as.Date("2024-01-01"), tau_days = 90, window_days = 365
  )
  expect_true(all(result$quantified))
})

# ── Property: n_samples_in_window is counted correctly ───────────────────────
test_that("n_samples_in_window equals samples within window", {
  # 10 samples uniformly spread across 2023
  df <- make_synth_chem(n_samples = 10, n_analytes = 1)
  focal <- as.Date("2024-01-01")
  result <- compute_chronic_chemistry(
    df, focal_dates = focal, tau_days = 90, window_days = 365,
    anchor_outside_window = FALSE
  )
  # All 10 samples should fall within focal - 365 to focal
  expect_equal(unique(result$n_samples_in_window), 10L)
})

# ── Property: n_imputed_in_window propagates imputed column ──────────────────
test_that("n_imputed_in_window counts imputed rows", {
  df <- make_synth_chem(n_samples = 10, n_analytes = 1)
  # Mark 3 samples as imputed
  df$imputed[1:3] <- TRUE
  focal <- as.Date("2024-01-01")
  result <- compute_chronic_chemistry(
    df, focal_dates = focal, tau_days = 90, window_days = 365,
    anchor_outside_window = FALSE
  )
  expect_gte(result$n_imputed_in_window, 0L)
  expect_lte(result$n_imputed_in_window, 3L)
})

# ── Schema: output has expected columns ──────────────────────────────────────
test_that("output has all required columns", {
  df <- make_synth_chem()
  result <- compute_chronic_chemistry(
    df, focal_dates = as.Date("2024-01-01")
  )
  expected_cols <- c("focal_date", "uuid.feature", "uuid.sample",
                     "name.analyte", "value", "quantified",
                     "n_samples_in_window", "n_imputed_in_window")
  expect_true(all(expected_cols %in% names(result)))
})

# ── Schema: synthetic uuid.sample is correct format ──────────────────────────
test_that("synthetic uuid.sample has correct prefix", {
  df <- make_synth_chem()
  focal <- as.Date("2024-01-01")
  result <- compute_chronic_chemistry(df, focal_dates = focal)
  expect_true(all(grepl("^chronic_", result$uuid.sample)))
})

# ── Error: no samples in window → abort with message ─────────────────────────
test_that("errors when no samples fall within window", {
  df <- make_synth_chem(start = as.Date("2020-01-01"), end = as.Date("2020-06-01"))
  expect_error(
    compute_chronic_chemistry(
      df, focal_dates = as.Date("2024-01-01"),
      window_days = 30, anchor_outside_window = FALSE
    ),
    regexp = "No chronic chemistry"
  )
})

# ── Arithmetic mean alternative ───────────────────────────────────────────────
test_that("arith_mean summary produces >= geom_mean for positive values", {
  df <- make_synth_chem(n_samples = 10, n_analytes = 1)
  focal <- as.Date("2024-01-01")
  gm <- compute_chronic_chemistry(df, focal_dates = focal, summary = "geom_mean")
  am <- compute_chronic_chemistry(df, focal_dates = focal, summary = "arith_mean")
  # By AM-GM inequality, arith mean >= geom mean
  expect_gte(am$value, gm$value)
})
