## Tests for draw_measurement_error() (Chunk 5).
## Stan-free. Verifies the lognormal expansion, carrier contract, and guards.
##
## Twelve properties tested:
##   1.  absent error_col     → clear error
##   2.  all-NA error_col     → identity (no draw_id added, schema unchanged)
##   3.  draws shape          → draw_id 1..N per eligible cell; N rows each
##   4.  geometric mean       → median of draws ≈ v; arithmetic mean > v
##   5.  CV vs SD equivalence → cv = sd/v gives identical results (same seed)
##   6.  spread monotonic     → larger CV → wider spread
##   7.  BDL untouched        → detected=FALSE rows not expanded
##   8.  imputed pass-through → imputed=TRUE rows not expanded (keep own draw_id)
##   9.  N from domain        → N inferred from existing draw domain
##   10. N conflict error     → mismatched ndraws vs domain aborts
##   11. ndraws required      → point frame + return="draws" without ndraws → error
##   12. seed reproducibility → same seed → identical draws; different → different

library(testthat)
library(leachatetools)

## ── Shared helpers ────────────────────────────────────────────────────────────

## Simple point frame: one analyte per sample, with optional error_col.
make_point_df <- function(vals       = c(1, 5, 10, 50, 100),
                           cv         = 0.10,
                           add_error  = TRUE,
                           analyte    = "Cu") {
  n <- length(vals)
  df <- tibble::tibble(
    sample_id = paste0("s", seq_len(n)),
    site_id   = "A",
    datetime  = as.Date("2024-01-01"),
    analyte   = analyte,
    value     = vals,
    detected  = TRUE
  )
  if (add_error) df[["value_cv"]] <- cv
  df
}

## Frame that already carries metals draws (simulates impute_chemistry output).
make_draws_df <- function(vals = c(10, 50), n_draws = 3L, cv = 0.10) {
  n <- length(vals)
  dplyr::bind_rows(
    # Drawn cells (imputed)
    tibble::tibble(
      sample_id    = rep(paste0("s", seq_len(n)), each = n_draws),
      site_id      = "A",
      datetime     = as.Date("2024-01-01"),
      analyte      = "Cu",
      value        = rep(vals, each = n_draws),
      detected     = TRUE,
      imputed      = TRUE,
      imputed_kind = "missing",
      draw_id      = rep(seq_len(n_draws), times = n),
      value_cv     = cv
    ),
    # Exact observed cells (co-analytes, draw_id = NA, not imputed)
    tibble::tibble(
      sample_id    = paste0("s", seq_len(n)),
      site_id      = "A",
      datetime     = as.Date("2024-01-01"),
      analyte      = "pH",
      value        = 7.5,
      detected     = TRUE,
      imputed      = FALSE,
      imputed_kind = "observed",
      draw_id      = NA_integer_,
      value_cv     = cv
    )
  )
}


## ── 1. Absent error_col → error ───────────────────────────────────────────────

test_that("absent error_col aborts with a clear message", {
  df <- make_point_df(add_error = FALSE)
  expect_error(
    draw_measurement_error(df, error_col = "value_cv", ndraws = 5L),
    "value_cv"
  )
})


## ── 2. All-NA error_col → identity ───────────────────────────────────────────

test_that("all-NA error_col returns df unchanged (no draw_id added)", {
  df <- make_point_df()
  df[["value_cv"]] <- NA_real_

  out <- draw_measurement_error(df, error_col = "value_cv", ndraws = 5L)

  expect_identical(out, df)
  expect_false("draw_id" %in% names(out))
})


## ── 3. Draws shape ────────────────────────────────────────────────────────────

test_that("eligible cells expand to N rows each with draw_id 1..N", {
  N  <- 5L
  df <- make_point_df(vals = c(1, 10, 100), cv = 0.10)
  out <- draw_measurement_error(df, error_col = "value_cv",
                                 ndraws = N, seed = 1L)

  # All rows are now expanded (all detected, all non-imputed, all have error)
  expect_true("draw_id" %in% names(out))
  expect_false(anyNA(out$draw_id))
  expect_equal(sort(unique(out$draw_id)), seq_len(N))
  expect_equal(nrow(out), 3L * N)
})

test_that("draw values are positive and finite", {
  df  <- make_point_df(cv = 0.20)
  out <- draw_measurement_error(df, error_col = "value_cv",
                                 ndraws = 20L, seed = 2L)
  expect_true(all(is.finite(out$value)))
  expect_true(all(out$value > 0))
})


## ── 4. Geometric-mean convention ─────────────────────────────────────────────

test_that("median of draws ≈ reported value (geometric-mean convention)", {
  # Single cell with value = 10, CV = 0.10; at large N median ≈ 10
  df <- tibble::tibble(
    sample_id = "s1", site_id = "A", datetime = as.Date("2024-01-01"),
    analyte = "Cu", value = 10, detected = TRUE, value_cv = 0.10
  )
  out <- draw_measurement_error(df, error_col = "value_cv",
                                 ndraws = 5000L, seed = 42L)
  expect_equal(median(out$value), 10, tolerance = 0.05)
})

test_that("arithmetic mean of draws slightly exceeds reported value", {
  df <- tibble::tibble(
    sample_id = "s1", site_id = "A", datetime = as.Date("2024-01-01"),
    analyte = "Cu", value = 10, detected = TRUE, value_cv = 0.20
  )
  out <- draw_measurement_error(df, error_col = "value_cv",
                                 ndraws = 5000L, seed = 42L)
  # mean = v * exp(sigma_log^2 / 2); sigma_log^2 = log(1 + 0.04) ≈ 0.039
  sigma_log2 <- log(1 + 0.20^2)
  expected_mean <- 10 * exp(sigma_log2 / 2)
  expect_equal(mean(out$value), expected_mean, tolerance = 0.05)
})


## ── 5. CV vs SD equivalence ───────────────────────────────────────────────────

test_that("error_type='sd' with sd=cv*v gives identical draws as error_type='cv'", {
  df_cv <- tibble::tibble(
    sample_id = "s1", site_id = "A", datetime = as.Date("2024-01-01"),
    analyte = "Cu", value = 10, detected = TRUE, value_cv = 0.15
  )
  df_sd <- dplyr::mutate(df_cv, value_sd = 0.15 * 10)

  out_cv <- draw_measurement_error(df_cv, error_col = "value_cv",
                                    error_type = "cv", ndraws = 50L, seed = 7L)
  out_sd <- draw_measurement_error(df_sd, error_col = "value_sd",
                                    error_type = "sd", ndraws = 50L, seed = 7L)

  expect_equal(out_cv$value, out_sd$value)
})


## ── 6. Spread monotonic in CV ─────────────────────────────────────────────────

test_that("larger CV produces wider spread of draws", {
  make_single <- function(cv) {
    df <- tibble::tibble(
      sample_id = "s1", site_id = "A", datetime = as.Date("2024-01-01"),
      analyte = "Cu", value = 10, detected = TRUE, value_cv = cv
    )
    draw_measurement_error(df, error_col = "value_cv",
                            ndraws = 1000L, seed = 99L)$value
  }
  d_small <- make_single(0.05)
  d_large <- make_single(0.30)

  expect_gt(stats::IQR(d_large), stats::IQR(d_small))
})


## ── 7. BDL cells untouched ────────────────────────────────────────────────────

test_that("BDL cells (detected=FALSE) are not expanded", {
  df <- tibble::tibble(
    sample_id = c("s1", "s2"),
    site_id   = "A",
    datetime  = as.Date("2024-01-01"),
    analyte   = "Cu",
    value     = c(0.5, 10),
    detected  = c(FALSE, TRUE),
    value_cv  = 0.10
  )
  out <- draw_measurement_error(df, error_col = "value_cv",
                                 ndraws = 4L, seed = 1L)

  # BDL row (s1) stays exact: one row with draw_id=NA
  s1 <- dplyr::filter(out, sample_id == "s1")
  expect_equal(nrow(s1), 1L)
  expect_true(is.na(s1$draw_id))

  # Detected row (s2) is expanded: 4 rows
  s2 <- dplyr::filter(out, sample_id == "s2")
  expect_equal(nrow(s2), 4L)
})


## ── 8. Imputed cells pass through ────────────────────────────────────────────

test_that("imputed=TRUE rows are not re-expanded by draw_measurement_error", {
  df <- make_draws_df(vals = c(10, 50), n_draws = 3L, cv = 0.10)

  # Cu rows are imputed (draw_id 1..3); pH rows are exact (draw_id NA)
  out <- draw_measurement_error(df, error_col = "value_cv", seed = 1L)

  # Cu rows unchanged (still 3 rows per sample, draw_id 1..3)
  cu <- dplyr::filter(out, analyte == "Cu")
  expect_equal(nrow(cu), 6L)  # 2 samples × 3 draws
  expect_equal(sort(unique(cu$draw_id)), 1:3)

  # pH rows expanded (exact + error → draws); 2 samples × 3 draws = 6
  ph <- dplyr::filter(out, analyte == "pH")
  expect_equal(nrow(ph), 6L)
  expect_false(anyNA(ph$draw_id))
})


## ── 9. N inferred from existing draw domain ───────────────────────────────────

test_that("N is inferred from existing draw domain; ndraws can be NULL", {
  N  <- 4L
  df <- make_draws_df(vals = c(10), n_draws = N, cv = 0.10)

  # pH rows (exact) should be expanded to 4 draws each
  out <- draw_measurement_error(df, error_col = "value_cv")

  ph <- dplyr::filter(out, analyte == "pH")
  expect_equal(sort(unique(ph$draw_id)), seq_len(N))
})


## ── 10. N conflict error ──────────────────────────────────────────────────────

test_that("ndraws conflicting with existing domain aborts", {
  N  <- 3L
  df <- make_draws_df(vals = c(10), n_draws = N, cv = 0.10)

  expect_error(
    draw_measurement_error(df, error_col = "value_cv", ndraws = N + 1L),
    "conflicts"
  )
})


## ── 11. ndraws required for point input ───────────────────────────────────────

test_that("point input without ndraws aborts with clear message", {
  df <- make_point_df()
  expect_error(
    draw_measurement_error(df, error_col = "value_cv"),
    "ndraws"
  )
})


## ── 12. Seed reproducibility ──────────────────────────────────────────────────

test_that("same seed produces identical draws", {
  df <- make_point_df(vals = c(10, 50))
  out1 <- draw_measurement_error(df, error_col = "value_cv",
                                   ndraws = 20L, seed = 42L)
  out2 <- draw_measurement_error(df, error_col = "value_cv",
                                   ndraws = 20L, seed = 42L)
  expect_identical(out1$value, out2$value)
})

test_that("different seeds produce different draws", {
  df <- make_point_df(vals = c(10, 50))
  out1 <- draw_measurement_error(df, error_col = "value_cv",
                                   ndraws = 20L, seed = 1L)
  out2 <- draw_measurement_error(df, error_col = "value_cv",
                                   ndraws = 20L, seed = 2L)
  expect_false(isTRUE(all.equal(out1$value, out2$value)))
})

test_that("draw output passes .draw_domain() ragged-N validation", {
  df  <- make_point_df(vals = c(1, 10))
  out <- draw_measurement_error(df, error_col = "value_cv",
                                 ndraws = 7L, seed = 5L)
  expect_no_error(leachatetools:::.draw_domain(out))
  expect_equal(leachatetools:::.draw_domain(out), 1:7)
})
