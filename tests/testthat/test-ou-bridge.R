## Tests for OU bridge core (issue #16, Chunk A).
## Stan-free, no external data required.
##
## Properties tested:
##   A. .estimate_ou_params
##      A1. degenerate: < 2 anchors → degenerate flag
##      A2. degenerate: constant S → degenerate flag
##      A3. returns theta/sigma2/gamma for a real series
##      A4. scale multiplies gamma and sigma2 but not theta
##      A5. theta = sigma2 / (2 * gamma) identity
##
##   B. .ou_bridge_factors + .ou_bridge_draw
##      B1. anchor positions carry epsilon = 0
##      B2. mid-gap variance > 0 (empirical from many draws)
##      B3. longer gap → wider mid-gap variance
##      B4. theta → 0 limit matches Brownian bridge variance (closed-form)
##      B5. variance saturates (trailing zone variance <= gamma)
##      B6. degenerate params (gamma=0) → zero-variance path
##      B7. single anchor → no bridge gaps, leading/trailing only
##      B8. seed reproducibility
##      B9. different seeds → different draws
##      B10. anchor pinch: empirical SD at anchor position ≈ 0

library(testthat)
library(leachatetools)

## ── helpers ──────────────────────────────────────────────────────────────────

# Simulate a regular anchor series with known OU parameters
make_anchors <- function(n = 20, theta = 0.05, gamma = 4, seed = 1L,
                          spacing_days = 14) {
  set.seed(seed)
  dates <- as.Date("2020-01-01") + (seq_len(n) - 1) * spacing_days
  # Discrete OU: S_{i+1} = phi * S_i + N(0, sigma2_disc)
  phi        <- exp(-theta * spacing_days)
  sigma2_disc <- gamma * (1 - phi^2)
  S <- numeric(n); S[1L] <- stats::rnorm(1, 0, sqrt(gamma))
  for (i in seq_len(n - 1)) S[i + 1] <- phi * S[i] + stats::rnorm(1, 0, sqrt(sigma2_disc))
  list(dates = dates, S = S, theta = theta, gamma = gamma)
}

# Empirical per-date variance from many .ou_bridge_draw() calls
empirical_var_at_idx <- function(factors, n_rep = 4000, seed = 42L) {
  set.seed(seed)
  mat <- replicate(n_rep, leachatetools:::.ou_bridge_draw(factors))
  apply(mat, 1, stats::var)
}


## ── A. .estimate_ou_params ────────────────────────────────────────────────────

test_that("A1: < 2 anchors → degenerate", {
  p <- leachatetools:::.estimate_ou_params(as.Date("2020-01-01"), 5)
  expect_true(p$degenerate)
})

test_that("A1: 0 anchors → degenerate", {
  p <- leachatetools:::.estimate_ou_params(as.Date(character()), numeric())
  expect_true(p$degenerate)
})

test_that("A2: constant S → degenerate (no diffusion detected)", {
  dates <- as.Date("2020-01-01") + 0:9 * 14
  S     <- rep(3, 10)
  p     <- leachatetools:::.estimate_ou_params(dates, S)
  expect_true(p$degenerate)
})

test_that("A3: real series returns positive finite params", {
  a <- make_anchors()
  p <- leachatetools:::.estimate_ou_params(a$dates, a$S)
  expect_false(p$degenerate)
  expect_true(p$theta  > 0 && is.finite(p$theta))
  expect_true(p$sigma2 > 0 && is.finite(p$sigma2))
  expect_true(p$gamma  > 0 && is.finite(p$gamma))
})

test_that("A4: scale multiplies gamma and sigma2, leaves theta unchanged", {
  a  <- make_anchors()
  p1 <- leachatetools:::.estimate_ou_params(a$dates, a$S, scale = 1)
  p2 <- leachatetools:::.estimate_ou_params(a$dates, a$S, scale = 2)
  expect_equal(p2$gamma,  p1$gamma  * 2)
  expect_equal(p2$sigma2, p1$sigma2 * 2)
  expect_equal(p2$theta,  p1$theta)       # theta unchanged
})

test_that("A5: MoM identity theta = sigma2 / (2 * gamma)", {
  a <- make_anchors()
  p <- leachatetools:::.estimate_ou_params(a$dates, a$S)
  expect_equal(p$theta, p$sigma2 / (2 * p$gamma), tolerance = 1e-10)
})


## ── B. .ou_bridge_factors + .ou_bridge_draw ───────────────────────────────────

test_that("B1: anchor positions carry epsilon = 0", {
  a       <- make_anchors(n = 5, spacing_days = 10)
  t_dates <- seq(a$dates[1], a$dates[5], by = "day")
  p       <- leachatetools:::.estimate_ou_params(a$dates, a$S)
  fac     <- leachatetools:::.ou_bridge_factors(a$dates, t_dates,
                                                p$theta, p$sigma2, p$gamma)
  set.seed(1L)
  eps <- leachatetools:::.ou_bridge_draw(fac)
  # Anchor positions are in target dates; epsilon must be exactly 0
  anchor_pos <- which(fac$anchor_in_target)
  expect_true(all(eps[anchor_pos] == 0))
})

test_that("B2: mid-gap empirical variance > 0 for non-degenerate params", {
  a       <- make_anchors(n = 3, spacing_days = 30)
  t_dates <- seq(a$dates[1], a$dates[3], by = "day")
  p       <- leachatetools:::.estimate_ou_params(a$dates, a$S)
  fac     <- leachatetools:::.ou_bridge_factors(a$dates, t_dates,
                                                p$theta, p$sigma2, p$gamma)
  e_var   <- empirical_var_at_idx(fac, n_rep = 3000)
  # Mid-gap index (approximately)
  mid_idx <- which(!fac$anchor_in_target)[floor(sum(!fac$anchor_in_target) / 2)]
  expect_gt(e_var[mid_idx], 0)
})

test_that("B3: longer gap produces wider mid-gap variance", {
  a_short  <- list(dates = as.Date("2020-01-01") + c(0, 14),
                   S = c(0, 0))
  a_long   <- list(dates = as.Date("2020-01-01") + c(0, 60),
                   S = c(0, 0))

  # Use same fixed OU params for both
  theta <- 0.02; sigma2 <- 0.5; gamma <- sigma2 / (2 * theta)

  make_fac <- function(anchors) {
    t_dates <- seq(anchors$dates[1], anchors$dates[2], by = "day")
    leachatetools:::.ou_bridge_factors(anchors$dates, t_dates, theta, sigma2, gamma)
  }
  fac_s <- make_fac(a_short)
  fac_l <- make_fac(a_long)

  var_s <- empirical_var_at_idx(fac_s, n_rep = 3000)
  var_l <- empirical_var_at_idx(fac_l, n_rep = 3000)

  # Maximum mid-gap variance is larger for the longer gap
  expect_gt(max(var_l), max(var_s))
})

test_that("B4: theta→0 limit matches Brownian bridge variance (closed-form)", {
  # With very small theta, the OU bridge variance at the midpoint should equal
  # sigma2 * L / 4 (Brownian bridge midpoint variance).
  L       <- 30                      # 30-day gap
  sigma2  <- 1.0
  # Use theta so small that theta*L << 1 (Brownian regime)
  theta   <- 1e-12
  gamma   <- sigma2 / (2 * theta)    # large but finite

  a_dates <- as.Date("2020-01-01") + c(0, L)
  t_dates <- as.Date("2020-01-01") + seq(0, L)

  fac <- leachatetools:::.ou_bridge_factors(a_dates, t_dates, theta, sigma2, gamma)

  # Extract the bridge factor for the single gap
  bridge_gap <- fac$gaps[[which(sapply(fac$gaps, `[[`, "type") == "bridge")]]
  chol_L     <- bridge_gap$chol_L
  # Var(eps_i) = sum of squared entries in row i of L (since Sigma = L L^T)
  emp_var    <- rowSums(chol_L^2)

  # Midpoint index within the interior (day 15, i.e. index 14 in 1-based interior)
  mid_i <- floor(nrow(chol_L) / 2)
  expected_var <- sigma2 * (mid_i) * (L - mid_i) / L   # Brownian bridge formula
  expect_equal(emp_var[mid_i], expected_var, tolerance = 1e-6)
})

test_that("B5: trailing zone variance ≤ gamma (saturates at marginal)", {
  a       <- make_anchors(n = 2, spacing_days = 14)
  # Extend target dates 90 days past last anchor
  t_dates <- seq(a$dates[1], a$dates[2] + 90, by = "day")
  p       <- leachatetools:::.estimate_ou_params(a$dates, a$S)
  fac     <- leachatetools:::.ou_bridge_factors(a$dates, t_dates,
                                                p$theta, p$sigma2, p$gamma)
  trail   <- fac$gaps[["trailing"]]
  expect_false(is.null(trail))
  # All one-sided SDs should be ≤ sqrt(gamma)
  expect_true(all(trail$sd_vec <= sqrt(p$gamma) + 1e-10))
})

test_that("B6: degenerate params (gamma=0) → zero epsilon path", {
  a       <- make_anchors(n = 5)
  t_dates <- seq(a$dates[1], a$dates[5], by = "day")
  # gamma = 0 → degenerate
  fac     <- leachatetools:::.ou_bridge_factors(a$dates, t_dates, 0, 0, 0)
  eps     <- leachatetools:::.ou_bridge_draw(fac)
  expect_true(all(eps == 0))
})

test_that("B7: single anchor → no bridge gaps, only leading/trailing", {
  single_anchor <- as.Date("2020-06-01")
  t_dates       <- seq(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day")
  theta <- 0.03; sigma2 <- 0.4; gamma <- sigma2 / (2 * theta)
  fac <- leachatetools:::.ou_bridge_factors(single_anchor, t_dates,
                                             theta, sigma2, gamma)
  types <- sapply(fac$gaps, `[[`, "type")
  expect_false("bridge" %in% types)
  expect_true("one_sided" %in% types)
})

test_that("B8: same seed → identical draws", {
  a       <- make_anchors()
  t_dates <- seq(a$dates[1], a$dates[length(a$dates)], by = "day")
  p       <- leachatetools:::.estimate_ou_params(a$dates, a$S)
  fac     <- leachatetools:::.ou_bridge_factors(a$dates, t_dates,
                                                p$theta, p$sigma2, p$gamma)
  set.seed(7L); eps1 <- leachatetools:::.ou_bridge_draw(fac)
  set.seed(7L); eps2 <- leachatetools:::.ou_bridge_draw(fac)
  expect_identical(eps1, eps2)
})

test_that("B9: different seeds → different draws", {
  a       <- make_anchors()
  t_dates <- seq(a$dates[1], a$dates[length(a$dates)], by = "day")
  p       <- leachatetools:::.estimate_ou_params(a$dates, a$S)
  fac     <- leachatetools:::.ou_bridge_factors(a$dates, t_dates,
                                                p$theta, p$sigma2, p$gamma)
  set.seed(1L); eps1 <- leachatetools:::.ou_bridge_draw(fac)
  set.seed(2L); eps2 <- leachatetools:::.ou_bridge_draw(fac)
  expect_false(isTRUE(all.equal(eps1, eps2)))
})

test_that("B10: empirical SD at anchor positions ≈ 0 across many draws", {
  a       <- make_anchors(n = 5, spacing_days = 10)
  t_dates <- seq(a$dates[1], a$dates[5], by = "day")
  p       <- leachatetools:::.estimate_ou_params(a$dates, a$S)
  fac     <- leachatetools:::.ou_bridge_factors(a$dates, t_dates,
                                                p$theta, p$sigma2, p$gamma)
  e_var   <- empirical_var_at_idx(fac, n_rep = 3000)
  anchor_pos <- which(fac$anchor_in_target)
  # Variance at anchor positions should be (numerically) 0
  expect_true(all(e_var[anchor_pos] < 1e-20))
})


## ── G5: additional calibration tests ─────────────────────────────────────────

test_that("G5a: MoM gamma estimate within factor-4 of simulated truth", {
  ## Dense-enough grab series lets the quadratic-variation estimator converge
  ## to a rough multiple of the truth.  Factor-4 tolerance accounts for
  ## finite-sample and Euler-Maruyama discretisation bias.
  gamma_true <- 2.0
  theta_true <- 0.10

  ## Simulate discrete OU at 14-day spacing, 40 grabs.
  set.seed(2025L)
  n     <- 40L
  dates <- as.Date("2020-01-01") + (seq_len(n) - 1L) * 14L
  phi   <- exp(-theta_true * 14)
  s2d   <- gamma_true * (1 - phi^2)
  S     <- numeric(n); S[1L] <- stats::rnorm(1, 0, sqrt(gamma_true))
  for (i in seq_len(n - 1L))
    S[i + 1L] <- phi * S[i] + stats::rnorm(1, 0, sqrt(s2d))

  p <- leachatetools:::.estimate_ou_params(dates, S)
  expect_false(p$degenerate)
  expect_true(p$gamma >= gamma_true / 4 && p$gamma <= gamma_true * 4,
              label = sprintf("gamma_hat=%.3f; true=%.3f", p$gamma, gamma_true))
})


test_that("G5c: empirical OU bridge midpoint variance within factor-2 of theory", {
  ## Use known parameters so the theoretical conditional variance is exact.
  theta  <- 0.15
  gamma  <- 3.0
  sigma2 <- 2 * theta * gamma

  ## 20-day gap; midpoint at t=10.
  L       <- 20L
  a_dates <- as.Date(c("2021-01-01", "2021-01-21"))
  t_dates <- seq(a_dates[1L], a_dates[2L], by = "day")
  fac     <- leachatetools:::.ou_bridge_factors(a_dates, t_dates, theta, sigma2, gamma)

  set.seed(42L)
  draws <- replicate(1500L, leachatetools:::.ou_bridge_draw(fac))

  mid_idx  <- which(t_dates == as.Date("2021-01-11"))   # t = 10
  emp_var  <- stats::var(draws[mid_idx, ])

  ## Theoretical: Var = γ(1−u²)(1−v²)/(1−s²), u=v=exp(-θ·10), s=exp(-θ·20)
  u        <- exp(-theta * 10)
  s        <- exp(-theta * L)
  theo_var <- gamma * (1 - u^2)^2 / (1 - s^2)

  expect_true(emp_var >= theo_var / 2 && emp_var <= theo_var * 2,
              label = sprintf("emp_var=%.4f; theo_var=%.4f", emp_var, theo_var))
})
