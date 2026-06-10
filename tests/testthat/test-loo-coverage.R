## Tests for the leave-one-anchor-out coverage harness (issue #16 rework).
## Stan-free, simulated data.
##
## Internal API under test (R/kalman_bridge.R):
##   .loo_coverage_series(anchor_dates, anchor_S, interval, block, n_fit_min,
##                        scale) -> list(coverage, mean_width, n, n_held)
##   .loo_anchor_coverage(target_model, interval, block, ...)
##       -> tibble (analyte, tier, n, coverage, mean_width) + a pooled row.

library(testthat)
library(leachatetools)

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
  cov_ok <- leachatetools:::.loo_coverage_series(s$anchor_dates, s$anchor_S,
                                                 interval = 0.9, scale = 1)
  expect_gte(cov_ok$coverage, 0.7)
  expect_lte(cov_ok$coverage, 1.0)
  expect_gt(cov_ok$n, 0)

  cov_narrow <- leachatetools:::.loo_coverage_series(s$anchor_dates, s$anchor_S,
                                                     interval = 0.9, scale = 0.02)
  expect_lt(cov_narrow$coverage, cov_ok$coverage)        # harness detects under-spread
})

## ── L2: block-holdout ─────────────────────────────────────────────────────────

test_that("L2: block holdout runs and returns coverage/width/n", {
  s <- sim_ou(seed = 12L)
  out <- leachatetools:::.loo_coverage_series(s$anchor_dates, s$anchor_S,
                                              interval = 0.9, block = 3L)
  expect_true(all(c("coverage", "mean_width", "n", "n_held") %in% names(out)))
  expect_gt(out$n_held, 0)
  expect_true(is.finite(out$mean_width) && out$mean_width > 0)
})

## ── L3: wider scale -> wider mean interval width ─────────────────────────────

test_that("L3: mean interval width grows with scale", {
  s <- sim_ou(seed = 13L)
  w1 <- leachatetools:::.loo_coverage_series(s$anchor_dates, s$anchor_S,
                                             scale = 1)$mean_width
  w4 <- leachatetools:::.loo_coverage_series(s$anchor_dates, s$anchor_S,
                                             scale = 4)$mean_width
  expect_gt(w4, w1)
})

## ── L4: sparse / degenerate series handled without error ──────────────────────

test_that("L4: too-few anchors returns NA coverage, no crash", {
  out <- leachatetools:::.loo_coverage_series(as.Date("2021-01-01") + c(0, 14),
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
  res <- leachatetools:::.loo_anchor_coverage(fake_tm, interval = 0.9)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("analyte", "tier", "n", "coverage", "mean_width") %in%
                    names(res)))
  expect_true("Cu" %in% res$analyte && "Zn" %in% res$analyte)
  expect_true(any(res$analyte == "(pooled)"))
})
