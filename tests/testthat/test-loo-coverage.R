## Tests for the leave-one-anchor-out coverage harness (issue #16 rework).
## Stan-free, simulated data.
##
## Internal API under test (R/kalman_bridge.R):
##   .loo_coverage_series(anchor_dates, anchor_S, interval, block, n_fit_min,
##                        scale) -> list(coverage, mean_width, n, n_held)
##   .loo_anchor_coverage(target_model, interval, block, ...)
##       -> tibble (analyte, tier, n, coverage, mean_width) + a pooled row.

library(testthat)
library(hydroSense)

sim_ou <- function(n_days = 420L, theta = 0.05, gamma = 4, anchor_every = 10L,
                   seed = 1L, r = 1e-2) {
  set.seed(seed)
  phi <- exp(-theta); q <- gamma * (1 - phi^2)
  x <- numeric(n_days); x[1L] <- stats::rnorm(1, 0, sqrt(gamma))
  for (t in 2:n_days) x[t] <- phi * x[t - 1L] + stats::rnorm(1, 0, sqrt(q))
  dates   <- as.Date("2021-01-01") + (seq_len(n_days) - 1L)
  anc_idx <- seq(1L, n_days, by = anchor_every)
  list(anchor_dates = dates[anc_idx],
       anchor_S     = x[anc_idx] + stats::rnorm(length(anc_idx), 0, sqrt(r)))
}

## ── L1: calibration + sensitivity to mis-calibration ─────────────────────────

test_that("L1: coverage is reasonable, and too-narrow scale lowers it", {
  s <- sim_ou(seed = 11L)
  cov_ok <- hydroSense:::.loo_coverage_series(s$anchor_dates, s$anchor_S,
                                                 interval = 0.9, scale = 1)
  expect_gte(cov_ok$coverage, 0.7)
  expect_lte(cov_ok$coverage, 1.0)
  expect_gt(cov_ok$n, 0)

  cov_narrow <- hydroSense:::.loo_coverage_series(s$anchor_dates, s$anchor_S,
                                                     interval = 0.9, scale = 0.02)
  expect_lt(cov_narrow$coverage, cov_ok$coverage)        # harness detects under-spread
})

## ── L2: block-holdout ─────────────────────────────────────────────────────────

test_that("L2: block holdout runs and returns coverage/width/n", {
  s <- sim_ou(seed = 12L)
  out <- hydroSense:::.loo_coverage_series(s$anchor_dates, s$anchor_S,
                                              interval = 0.9, block = 3L)
  expect_true(all(c("coverage", "mean_width", "n", "n_held") %in% names(out)))
  expect_gt(out$n_held, 0)
  expect_true(is.finite(out$mean_width) && out$mean_width > 0)
})

## ── L3: wider scale -> wider mean interval width ─────────────────────────────

test_that("L3: mean interval width grows with scale", {
  s <- sim_ou(seed = 13L)
  w1 <- hydroSense:::.loo_coverage_series(s$anchor_dates, s$anchor_S,
                                             scale = 1)$mean_width
  w4 <- hydroSense:::.loo_coverage_series(s$anchor_dates, s$anchor_S,
                                             scale = 4)$mean_width
  expect_gt(w4, w1)
})

## ── L4: sparse / degenerate series handled without error ──────────────────────

test_that("L4: too-few anchors returns NA coverage, no crash", {
  out <- hydroSense:::.loo_coverage_series(as.Date("2021-01-01") + c(0, 14),
                                              c(1, 2), interval = 0.9)
  expect_true(is.na(out$coverage) || out$n == 0L)
})

## ── L5: target_model wrapper schema ───────────────────────────────────────────

test_that("L5: .loo_anchor_coverage returns per-analyte + pooled rows", {
  mk_anchors <- function(seed) {
    s <- sim_ou(anchor_every = 12L, seed = seed)
    tibble::tibble(date = s$anchor_dates, S = s$anchor_S)
  }
  fake_tm <- list(models = list(
    Cu = list(tier = "model",  anchors = mk_anchors(21L)),
    Zn = list(tier = "bridge", anchors = mk_anchors(22L))
  ))
  res <- hydroSense:::.loo_anchor_coverage(fake_tm, interval = 0.9)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("analyte", "tier", "n", "coverage", "mean_width") %in%
                    names(res)))
  expect_true("Cu" %in% res$analyte && "Zn" %in% res$analyte)
  expect_true(any(res$analyte == "(pooled)"))
})


## ── L6: .resolve_tier_scale (tier-aware ou_scale) ────────────────────────────

test_that("L6: .resolve_tier_scale picks per-tier scales or falls back to 1", {
  rs <- hydroSense:::.resolve_tier_scale
  expect_equal(rs(2, "model"), 2)                              # scalar -> all
  expect_equal(rs(c(model = 1, bridge = 3), "bridge"), 3)      # named -> tier
  expect_equal(rs(c(model = 1, bridge = 3), "model"), 1)
  expect_equal(rs(c(bridge = 3), "model"), 1)                  # absent tier -> 1
  expect_equal(rs(c(model = 2), NA_character_), 1)             # NA tier -> 1
})


## ── L7: public loo_anchor_coverage() ─────────────────────────────────────────

test_that("L7: loo_anchor_coverage() validates input and returns the table", {
  expect_error(loo_anchor_coverage(list()), "target_model")

  mk_anch <- function(seed) {
    s <- sim_ou(anchor_every = 12L, seed = seed)
    tibble::tibble(date = s$anchor_dates, S = s$anchor_S)
  }
  tm <- structure(list(models = list(
    Cu = list(tier = "model",  anchors = mk_anch(21L)),
    Zn = list(tier = "bridge", anchors = mk_anch(22L))
  )), class = "target_model")

  res <- loo_anchor_coverage(tm, interval = 0.9)
  expect_s3_class(res, "tbl_df")
  expect_true(all(c("analyte", "tier", "n", "coverage", "mean_width") %in%
                    names(res)))
  expect_true(any(res$analyte == "(pooled)"))
  expect_setequal(stats::na.omit(res$tier), c("model", "bridge"))
  expect_error(loo_anchor_coverage(tm, interval = 2), "interval")
})
