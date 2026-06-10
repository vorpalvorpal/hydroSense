## Tests for the Kalman/OU residual smoother core (issue #16 rework).
## Stan-free, simulated data. Replaces test-ou-bridge.R.
##
## Internal API under test (R/kalman_bridge.R):
##   .estimate_ou_kalman_params(anchor_dates, anchor_S, n_fit_min, scale)
##       -> list(theta, gamma, tier)   tier in c("mle","fallback","degenerate")
##   .build_kalman_model(grid_dates, anchor_dates, anchor_S, theta, gamma,
##                       r_vec = NULL, q_mult = NULL)  -> SSModel or NULL
##   .kalman_smooth(model) -> list(mean, var)
##   .kalman_draw(model, nsim) -> matrix [n_grid x nsim]
##   .residual_smoother(anchor_dates, anchor_S, target_dates, z_hydro, kappa,
##                      n_fit_min, scale, r_vec) -> list(grid_dates, params,
##                      model, mean, var)

library(testthat)
library(leachatetools)

## ── helpers ──────────────────────────────────────────────────────────────────

sim_ou <- function(n_days = 210L, theta = 0.05, gamma = 4, anchor_every = 14L,
                   seed = 1L, r = 1e-3) {
  set.seed(seed)
  phi <- exp(-theta); q <- gamma * (1 - phi^2)
  x <- numeric(n_days); x[1L] <- stats::rnorm(1, 0, sqrt(gamma))
  for (t in 2:n_days) x[t] <- phi * x[t - 1L] + stats::rnorm(1, 0, sqrt(q))
  dates   <- as.Date("2021-01-01") + (seq_len(n_days) - 1L)
  anc_idx <- seq(1L, n_days, by = anchor_every)
  list(
    dates        = dates,
    x            = x,
    target_dates = dates,
    anchor_dates = dates[anc_idx],
    anchor_S     = x[anc_idx] + stats::rnorm(length(anc_idx), 0, sqrt(r)),
    anc_idx      = anc_idx
  )
}

## ── K1: MLE parameter recovery ────────────────────────────────────────────────

test_that("K1: estimate recovers gamma well and a finite positive theta (mle tier)", {
  s <- sim_ou(n_days = 420L, theta = 0.05, gamma = 4, anchor_every = 10L, seed = 7L)
  p <- leachatetools:::.estimate_ou_kalman_params(s$anchor_dates, s$anchor_S,
                                                  n_fit_min = 8L)
  expect_identical(p$tier, "mle")
  # gamma (marginal variance) is the well-identified moment
  expect_gt(p$gamma, 1.5)
  expect_lt(p$gamma, 9)
  # theta finite, positive, bounded by the spacing-derived correlation length
  expect_true(is.finite(p$theta) && p$theta > 0)
})

## ── K2: fallback ladder ───────────────────────────────────────────────────────

test_that("K2: <n_fit_min anchors -> fallback tier; <2 or constant -> degenerate", {
  s <- sim_ou(anchor_every = 40L, seed = 2L)           # ~6 anchors over 210d
  expect_lt(length(s$anchor_S), 8L)
  p_fb <- leachatetools:::.estimate_ou_kalman_params(s$anchor_dates, s$anchor_S,
                                                     n_fit_min = 8L)
  expect_identical(p_fb$tier, "fallback")
  expect_gt(p_fb$gamma, 0)

  one <- leachatetools:::.estimate_ou_kalman_params(as.Date("2021-01-01"), 3)
  expect_identical(one$tier, "degenerate")

  const <- leachatetools:::.estimate_ou_kalman_params(
    as.Date("2021-01-01") + c(0, 14, 28, 42, 56, 70, 84, 98), rep(2, 8))
  expect_identical(const$tier, "degenerate")
})

## ── K3: pinch at anchors ──────────────────────────────────────────────────────

test_that("K3: with tiny r, smoother pins mean to obs and var ~0 at anchors", {
  s <- sim_ou(seed = 3L)
  p <- leachatetools:::.estimate_ou_kalman_params(s$anchor_dates, s$anchor_S)
  m <- leachatetools:::.build_kalman_model(s$target_dates, s$anchor_dates,
                                           s$anchor_S, p$theta, p$gamma,
                                           r_vec = rep(1e-8, length(s$anchor_S)))
  sm <- leachatetools:::.kalman_smooth(m)
  anc_pos <- match(s$anchor_dates, s$target_dates)
  expect_equal(sm$mean[anc_pos], s$anchor_S, tolerance = 1e-3)
  expect_true(all(sm$var[anc_pos] < 1e-3))
})

## ── K4: balloon mid-gap, wider for longer gaps ────────────────────────────────

test_that("K4: mid-gap variance exceeds anchor variance and grows with gap length", {
  s <- sim_ou(seed = 4L)
  p <- leachatetools:::.estimate_ou_kalman_params(s$anchor_dates, s$anchor_S)
  m <- leachatetools:::.build_kalman_model(s$target_dates, s$anchor_dates,
                                           s$anchor_S, p$theta, p$gamma,
                                           r_vec = rep(1e-6, length(s$anchor_S)))
  sm <- leachatetools:::.kalman_smooth(m)
  anc_pos <- match(s$anchor_dates, s$target_dates)
  mid_pos <- anc_pos[1L] + 7L                       # mid first 14-day gap
  expect_gt(sm$var[mid_pos], max(sm$var[anc_pos]))

  # Longer gap -> larger mid-gap variance (compare a 14d gap to a 28d gap)
  s2 <- sim_ou(anchor_every = 28L, seed = 4L)
  p2 <- leachatetools:::.estimate_ou_kalman_params(s2$anchor_dates, s2$anchor_S)
  m2 <- leachatetools:::.build_kalman_model(s2$target_dates, s2$anchor_dates,
                                            s2$anchor_S, p2$theta, p2$gamma,
                                            r_vec = rep(1e-6, length(s2$anchor_S)))
  sm2 <- leachatetools:::.kalman_smooth(m2)
  expect_gt(max(sm2$var), max(sm$var))
})

## ── K5: draws ────────────────────────────────────────────────────────────────

test_that("K5: simulateSSM draw mean tracks the smoother mean and is reproducible", {
  s <- sim_ou(seed = 5L)
  p <- leachatetools:::.estimate_ou_kalman_params(s$anchor_dates, s$anchor_S)
  m <- leachatetools:::.build_kalman_model(s$target_dates, s$anchor_dates,
                                           s$anchor_S, p$theta, p$gamma,
                                           r_vec = rep(1e-6, length(s$anchor_S)))
  sm <- leachatetools:::.kalman_smooth(m)

  set.seed(99L); d1 <- leachatetools:::.kalman_draw(m, nsim = 400L)
  expect_equal(nrow(d1), length(s$target_dates))
  expect_equal(ncol(d1), 400L)
  expect_gt(stats::cor(rowMeans(d1), sm$mean), 0.97)

  # draws pinch at anchors (tiny spread there)
  anc_pos <- match(s$anchor_dates, s$target_dates)
  expect_true(mean(apply(d1[anc_pos, , drop = FALSE], 1, stats::sd)) <
                mean(apply(d1[anc_pos + 7L, , drop = FALSE], 1, stats::sd)))

  set.seed(99L); d2 <- leachatetools:::.kalman_draw(m, nsim = 400L)
  expect_identical(d1, d2)                              # seed reproducibility
})

## ── K6: hydrology-modulated process variance ──────────────────────────────────

test_that("K6: q-modulation widens storm gaps; kappa=0 is stationary", {
  s <- sim_ou(seed = 6L)
  p <- leachatetools:::.estimate_ou_kalman_params(s$anchor_dates, s$anchor_S)

  n  <- length(s$target_dates)
  z  <- rep(0, n)
  storm_gap <- (s$anc_idx[3L] + 1L):(s$anc_idx[4L] - 1L)   # one whole gap
  calm_gap  <- (s$anc_idx[5L] + 1L):(s$anc_idx[6L] - 1L)
  z[storm_gap] <- 3                                        # high flow in storm gap

  q_mult <- exp(0.5 * z)                                   # kappa = 0.5
  m_mod <- leachatetools:::.build_kalman_model(
    s$target_dates, s$anchor_dates, s$anchor_S, p$theta, p$gamma,
    r_vec = rep(1e-6, length(s$anchor_S)), q_mult = q_mult)
  sm_mod <- leachatetools:::.kalman_smooth(m_mod)

  m_flat <- leachatetools:::.build_kalman_model(
    s$target_dates, s$anchor_dates, s$anchor_S, p$theta, p$gamma,
    r_vec = rep(1e-6, length(s$anchor_S)), q_mult = rep(1, n))
  sm_flat <- leachatetools:::.kalman_smooth(m_flat)

  # storm gap widens relative to flat; calm gap unchanged
  expect_gt(max(sm_mod$var[storm_gap]), max(sm_flat$var[storm_gap]))
  expect_equal(max(sm_mod$var[calm_gap]), max(sm_flat$var[calm_gap]),
               tolerance = 1e-6)
})

## ── K7: small-theta approaches the Brownian-bridge variance shape ──────────────

test_that("K7: small theta gives a Brownian-bridge-shaped within-gap variance", {
  # Single long gap between two anchors; tiny theta -> var ~ sigma2 * t(L-t)/L
  grid <- as.Date("2021-01-01") + 0:40
  anc  <- grid[c(1L, 41L)]
  aS   <- c(0, 0)
  gamma <- 100; theta <- 1e-6
  m <- leachatetools:::.build_kalman_model(grid, anc, aS, theta, gamma,
                                           r_vec = rep(1e-8, 2))
  sm <- leachatetools:::.kalman_smooth(m)
  L <- 40; tt <- 0:40
  shape <- tt * (L - tt) / L
  shape_var <- sm$var
  # correlation between the realised variance profile and the Brownian shape ~1
  expect_gt(stats::cor(shape[2:40], shape_var[2:40]), 0.999)
})

## ── K8: clip to grab span ─────────────────────────────────────────────────────

test_that("K8: residual_smoother clips the grid to the anchor span", {
  s <- sim_ou(seed = 8L)
  # target range extends beyond the anchors on both sides
  tgt <- c(s$anchor_dates[1L] - 30L, s$target_dates,
           s$anchor_dates[length(s$anchor_dates)] + 30L)
  rs <- leachatetools:::.residual_smoother(s$anchor_dates, s$anchor_S, tgt)
  expect_gte(min(rs$grid_dates), min(s$anchor_dates))
  expect_lte(max(rs$grid_dates), max(s$anchor_dates))
})

## ── K9: degenerate -> zero-uncertainty path ───────────────────────────────────

test_that("K9: degenerate params produce a finite mean with ~zero variance", {
  dates <- as.Date("2021-01-01") + 0:100
  anc   <- as.Date("2021-01-01")
  rs <- leachatetools:::.residual_smoother(anc, 5, dates)
  expect_true(all(is.finite(rs$mean)))
  expect_true(all(rs$var < 1e-8))
})

## ── K10: ou_scale (scale) knob widens the band, preserves correlation length ───

test_that("K10: scale multiplies gamma (wider var) without changing theta", {
  s <- sim_ou(seed = 10L)
  p1 <- leachatetools:::.estimate_ou_kalman_params(s$anchor_dates, s$anchor_S,
                                                   scale = 1)
  p2 <- leachatetools:::.estimate_ou_kalman_params(s$anchor_dates, s$anchor_S,
                                                   scale = 4)
  expect_equal(p2$theta, p1$theta, tolerance = 1e-8)
  expect_equal(p2$gamma, 4 * p1$gamma, tolerance = 1e-6)
})
