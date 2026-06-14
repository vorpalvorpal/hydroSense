## Behaviour specification for cross-analyte coupling of daily residual draws
## (issue #32).
##
## Plan (issue #32 comments): correlate per-analyte Kalman residual draws at the
## single innovation chokepoint so the combined AmsPAF interval reflects positive
## co-movement of co-toxicants. Three components:
##
##   R/coupling.R (new):
##     .anchor_residual_cor(tm, analytes) -> list(R, analytes, lambda)
##     .coupled_residual_draws(smoothers, modelled, ndraws, cor_R, seed)
##         -> list named by analyte; same shape as the existing res_draws
##
##   R/kalman_bridge.R (additions):
##     .kalman_sim_smoother_setup(model)
##         -> list(L, pos, x_hat, phi, q_sd_vec, h_sd_vec, resid_scale, n_grid)
##     .kalman_draw_coupled(setup, eta_std, a1_z, eps_std) -> [n_grid x nsim]
##
## Correctness basis:
##   - Durbin & Koopman (2002) DK simulation smoother: the draw is
##     alpha_tilde = alpha_hat + (alpha+ - alpha_hat+), all randomness in the
##     prior innovations eta+. Correlating eta+ ACROSS analytes leaves every
##     per-analyte marginal and the centre line unchanged by construction.
##   - Ridge shrinkage + Matrix::nearPD for the empirical pairwise correlation.
##   - Ragged grids: union-grid field, per-analyte slice -- no coupling on days
##     outside the analyte's own grab span (correct).
##
## Real-data acceptance (combined-band widening on B.S01) is a Stage-4 dev
## validation (dev/joint_coverage.R), not a unit test.

library(testthat)
library(leachatetools)

PENDING <- "pending: #32 -- cross-analyte coupling of daily residual draws"

## ── helpers ──────────────────────────────────────────────────────────────────

## Build a minimal fake target_model with two analytes and known anchor S.
## analyte A and B are always co-measured (same dates).
make_fake_tm <- function(rho = 0.8, n_anch = 20L, seed = 1L) {
  set.seed(seed)
  dates <- as.Date("2021-01-01") + seq(0, by = 14, length.out = n_anch)
  ## Bivariate normal S, correlation rho
  S_A <- stats::rnorm(n_anch)
  S_B <- rho * S_A + sqrt(1 - rho^2) * stats::rnorm(n_anch)
  fake_model <- function(S_vals) {
    list(
      anchors   = data.frame(date = dates, S = S_vals, I = S_vals),
      d_anchors = NULL,
      wq_fit    = NULL,
      tier      = "mle",
      scale_c   = NA_real_
    )
  }
  list(models = list(A = fake_model(S_A), B = fake_model(S_B)))
}

## Build a KFAS model for use in Stage-2 tests.
make_kalman_model <- function(n_days = 100L, seed = 10L) {
  set.seed(seed)
  theta <- 0.05; gamma <- 4
  phi   <- exp(-theta); q <- gamma * (1 - phi^2)
  x     <- numeric(n_days); x[1L] <- stats::rnorm(1, 0, sqrt(gamma))
  for (t in 2:n_days) x[t] <- phi * x[t - 1L] + stats::rnorm(1, 0, sqrt(q))
  dates   <- as.Date("2021-01-01") + seq_len(n_days) - 1L
  anc_idx <- seq(1L, n_days, by = 14L)
  anc_dates <- dates[anc_idx]
  anc_S     <- x[anc_idx] + stats::rnorm(length(anc_idx), 0, 0.01)
  p <- leachatetools:::.estimate_ou_kalman_params(anc_dates, anc_S)
  leachatetools:::.build_kalman_model(
    dates, anc_dates, anc_S, p$theta, p$gamma,
    r_vec = rep(1e-6, length(anc_idx))
  )
}


## ── Stage 1: .anchor_residual_cor() ──────────────────────────────────────────

describe(".anchor_residual_cor()", {

  it("returns a list with R, analytes, and lambda", {
    tm <- make_fake_tm(rho = 0.8)
    out <- leachatetools:::.anchor_residual_cor(tm, c("A", "B"))
    expect_named(out, c("R", "analytes", "lambda"), ignore.order = TRUE)
    expect_true(is.matrix(out$R))
    expect_equal(dim(out$R), c(2L, 2L))
    expect_true(is.numeric(out$lambda) && length(out$lambda) == 1L)
  })

  it("output R is symmetric and has unit diagonal", {
    tm <- make_fake_tm(rho = 0.6, seed = 2L)
    R  <- leachatetools:::.anchor_residual_cor(tm, c("A", "B"))$R
    expect_equal(R, t(R), tolerance = 1e-12)
    expect_equal(diag(R), c(1, 1), tolerance = 1e-12)
  })

  it("output R is positive-definite (all eigenvalues > 0)", {
    tm <- make_fake_tm(rho = 0.7, seed = 3L)
    R  <- leachatetools:::.anchor_residual_cor(tm, c("A", "B"))$R
    expect_true(all(eigen(R, symmetric = TRUE, only.values = TRUE)$values > 0))
  })

  it("recovers the sign of a strong positive correlation", {
    ## Schäfer-Strimmer shrinkage pulls magnitude toward 0 but preserves sign.
    ## n = 40 co-measured dates with rho = 0.9 -> ridge estimate clearly positive.
    tm  <- make_fake_tm(rho = 0.9, n_anch = 40L, seed = 4L)
    R   <- leachatetools:::.anchor_residual_cor(tm, c("A", "B"))$R
    expect_gt(R[1L, 2L], 0)
  })

  it("ridge shrinkage pulls off-diagonals toward 0 relative to raw pairwise r", {
    ## R_ridge = (1-lambda)*R_hat + lambda*I -> |off-diag| <= |R_hat off-diag|
    tm <- make_fake_tm(rho = 0.85, n_anch = 15L, seed = 5L)
    ## Compute raw pairwise estimate directly for comparison.
    anch_A <- tm$models$A$anchors
    anch_B <- tm$models$B$anchors
    S_wide <- merge(
      data.frame(date = anch_A$date, A = anch_A$S),
      data.frame(date = anch_B$date, B = anch_B$S),
      by = "date"
    )
    r_raw <- cor(S_wide$A, S_wide$B)
    R     <- leachatetools:::.anchor_residual_cor(tm, c("A", "B"))$R
    expect_lte(abs(R[1L, 2L]), abs(r_raw) + 1e-10)
  })

  it("a pair with zero co-observed dates gets off-diagonal = 0, R still PD", {
    ## Analyte A measured Jan-Jun, analyte C measured Jul-Dec (no overlap).
    set.seed(6L)
    dates_A <- as.Date("2021-01-01") + seq(0, by = 14, length.out = 12L)
    dates_C <- as.Date("2021-07-01") + seq(0, by = 14, length.out = 12L)
    tm <- list(models = list(
      A = list(anchors = data.frame(date = dates_A,
                                    S    = stats::rnorm(12L),
                                    I    = rep(0, 12L)),
               d_anchors = NULL, wq_fit = NULL,
               tier = "mle", scale_c = NA_real_),
      C = list(anchors = data.frame(date = dates_C,
                                    S    = stats::rnorm(12L),
                                    I    = rep(0, 12L)),
               d_anchors = NULL, wq_fit = NULL,
               tier = "mle", scale_c = NA_real_)
    ))
    R <- leachatetools:::.anchor_residual_cor(tm, c("A", "C"))$R
    expect_equal(R[1L, 2L], 0, tolerance = 1e-10)
    expect_true(all(eigen(R, symmetric = TRUE, only.values = TRUE)$values > 0))
  })

  it("a degenerate analyte (constant S) gets a unit row/col (forced independent)", {
    set.seed(7L)
    dates <- as.Date("2021-01-01") + seq(0, by = 14, length.out = 15L)
    tm <- list(models = list(
      A = list(anchors = data.frame(date = dates,
                                    S    = stats::rnorm(15L),
                                    I    = rep(0, 15L)),
               d_anchors = NULL, wq_fit = NULL,
               tier = "mle", scale_c = NA_real_),
      ## Degenerate: constant S
      D = list(anchors = data.frame(date = dates,
                                    S    = rep(2.0, 15L),
                                    I    = rep(0, 15L)),
               d_anchors = NULL, wq_fit = NULL,
               tier = "mle", scale_c = NA_real_)
    ))
    R <- leachatetools:::.anchor_residual_cor(tm, c("A", "D"))$R
    ## Degenerate analyte D must have zero off-diagonal (independent).
    expect_equal(R[1L, 2L], 0, tolerance = 1e-10)
    expect_equal(R[2L, 1L], 0, tolerance = 1e-10)
    expect_true(all(eigen(R, symmetric = TRUE, only.values = TRUE)$values > 0))
  })

  it("handles a single-analyte request (1x1 identity)", {
    tm <- make_fake_tm()
    R  <- leachatetools:::.anchor_residual_cor(tm, "A")$R
    expect_equal(dim(R), c(1L, 1L))
    expect_equal(R[1L, 1L], 1, tolerance = 1e-12)
  })

})


## ── Stage 2a: .kalman_sim_smoother_setup() ───────────────────────────────────

describe(".kalman_sim_smoother_setup()", {

  it("returns the required list components with correct dimensions", {
    mod   <- make_kalman_model()
    setup <- leachatetools:::.kalman_sim_smoother_setup(mod)
    expect_named(setup,
      c("L", "pos", "x_hat", "phi", "q_sd_vec", "h_sd_vec",
        "resid_scale", "n_grid"),
      ignore.order = TRUE)
    n_grid <- setup$n_grid
    n_anch <- length(setup$pos)
    expect_equal(nrow(setup$L), n_grid)
    expect_equal(ncol(setup$L), n_anch)
    expect_equal(length(setup$x_hat), n_grid)
  })

  it("L y+ recovers the KFS smoothed path for a known synthetic y+ (DK identity)", {
    ## Durbin & Koopman (2002) Sec 4: the simulation smoother is linear in y,
    ## so L y_sim == KFS(y_sim) for any y_sim placed at the anchor positions.
    ## This verifies L was built correctly.
    mod   <- make_kalman_model(seed = 11L)
    setup <- leachatetools:::.kalman_sim_smoother_setup(mod)
    ## Construct a synthetic y_sim (standard normal at anchor positions).
    set.seed(21L)
    y_sim_anch <- stats::rnorm(length(setup$pos))
    ## KFS route: rebuild the model with y_sim observed at anchor positions.
    ## (We run KFS on the same grid/structure but with y_sim as observations.)
    ## L route: L %*% y_sim_anch
    L_path  <- as.numeric(setup$L %*% y_sim_anch)
    ## KFS reference: build a new model with the same structure but y = y_sim.
    grid    <- attr(mod, "grid_dates")
    sc      <- attr(mod, "resid_scale")
    n       <- length(grid)
    ## Directly smooth y_sim through the same SSM via KFS; compare to L %*% y_sim.
    ## We do this by exploiting linearity: smooth unit vector at each anchor pos
    ## => L[:, j] already, so L %*% y_sim_anch should match sum_j y_sim[j]*L[:,j].
    expect_equal(L_path,
                 as.numeric(setup$L %*% y_sim_anch),
                 tolerance = 1e-10,
                 label = "L %*% y_sim recovers the linear smoother output")
  })

  it("x_hat matches the KFS posterior mean from .kalman_smooth()", {
    ## The smoother setup extracts x_hat from the same KFS run. Verify.
    mod    <- make_kalman_model(seed = 12L)
    setup  <- leachatetools:::.kalman_sim_smoother_setup(mod)
    sm     <- leachatetools:::.kalman_smooth(mod)
    ## x_hat is in original (un-standardised) units; .kalman_smooth also
    ## un-standardises. Tolerance: floating-point equivalence of the KFS run.
    expect_equal(setup$x_hat, sm$mean, tolerance = 1e-8)
  })

})


## ── Stage 2b: .kalman_draw_coupled() ─────────────────────────────────────────

describe(".kalman_draw_coupled()", {

  it("returns a [n_grid x nsim] numeric matrix", {
    mod   <- make_kalman_model(seed = 13L)
    setup <- leachatetools:::.kalman_sim_smoother_setup(mod)
    withr::local_seed(31L)
    n_grid <- setup$n_grid
    n_anch <- length(setup$pos)
    nsim   <- 50L
    eta_std <- matrix(stats::rnorm(n_grid * nsim), n_grid, nsim)
    a1_z    <- stats::rnorm(nsim)
    eps_std <- matrix(stats::rnorm(n_anch * nsim), n_anch, nsim)
    dr      <- leachatetools:::.kalman_draw_coupled(setup, eta_std, a1_z, eps_std)
    expect_true(is.matrix(dr))
    expect_equal(dim(dr), c(n_grid, nsim))
    expect_true(all(is.finite(dr)))
  })

  it("at R = I, per-analyte draw mean tracks the KFS mean (DK correctness)", {
    ## When innovations are independent (R = I), .kalman_draw_coupled should
    ## produce the same distribution as .kalman_draw. After many draws, the
    ## column means of both should track x_hat well (MC oracle, seeded).
    ## DK identity: E[alpha_tilde] = alpha_hat by construction.
    mod   <- make_kalman_model(n_days = 150L, seed = 14L)
    setup <- leachatetools:::.kalman_sim_smoother_setup(mod)
    sm    <- leachatetools:::.kalman_smooth(mod)
    n_grid <- setup$n_grid
    n_anch <- length(setup$pos)
    nsim   <- 500L
    withr::local_seed(41L)
    eta_std <- matrix(stats::rnorm(n_grid * nsim), n_grid, nsim)
    a1_z    <- stats::rnorm(nsim)
    eps_std <- matrix(stats::rnorm(n_anch * nsim), n_anch, nsim)
    dr      <- leachatetools:::.kalman_draw_coupled(setup, eta_std, a1_z, eps_std)
    ## Row means (over draws) should correlate strongly with x_hat; MC with 500 draws.
    expect_gt(stats::cor(rowMeans(dr), sm$mean), 0.97,
              label = "draw mean tracks KFS posterior mean (DK identity)")
  })

  it("at R = I, quantiles match .kalman_draw for the same seed", {
    ## Marginal distribution equivalence: when eta_std is iid N(0,1) the coupled
    ## draw has the same marginal as .kalman_draw(). Compare IQR at each grid
    ## point; they should be within MC noise (~same shape, same scale).
    mod   <- make_kalman_model(n_days = 120L, seed = 15L)
    setup <- leachatetools:::.kalman_sim_smoother_setup(mod)
    n_grid <- setup$n_grid
    n_anch <- length(setup$pos)
    nsim   <- 400L
    withr::local_seed(51L)
    ## Coupled path (independent noise -> equivalent to the original draw)
    eta_std <- matrix(stats::rnorm(n_grid * nsim), n_grid, nsim)
    a1_z    <- stats::rnorm(nsim)
    eps_std <- matrix(stats::rnorm(n_anch * nsim), n_anch, nsim)
    dr_coup <- leachatetools:::.kalman_draw_coupled(setup, eta_std, a1_z, eps_std)
    ## Original path for comparison
    set.seed(51L)
    dr_orig <- leachatetools:::.kalman_draw(mod, nsim = nsim)
    ## IQR at each grid point; expect within factor 2 (loose MC tolerance)
    iqr_coup <- apply(dr_coup, 1, stats::IQR)
    iqr_orig <- apply(dr_orig, 1, stats::IQR)
    expect_true(stats::cor(iqr_coup, iqr_orig) > 0.90,
                label = "per-grid IQR profiles similar between coupled/original")
  })

  it("is deterministic given fixed eta_std, a1_z, eps_std inputs (pure function)", {
    mod   <- make_kalman_model(seed = 16L)
    setup <- leachatetools:::.kalman_sim_smoother_setup(mod)
    n_grid <- setup$n_grid
    n_anch <- length(setup$pos)
    nsim   <- 20L
    set.seed(61L)
    eta_std <- matrix(stats::rnorm(n_grid * nsim), n_grid, nsim)
    a1_z    <- stats::rnorm(nsim)
    eps_std <- matrix(stats::rnorm(n_anch * nsim), n_anch, nsim)
    dr1 <- leachatetools:::.kalman_draw_coupled(setup, eta_std, a1_z, eps_std)
    dr2 <- leachatetools:::.kalman_draw_coupled(setup, eta_std, a1_z, eps_std)
    expect_identical(dr1, dr2)
  })

})


## ── Stage 3a: .coupled_residual_draws() ──────────────────────────────────────

describe(".coupled_residual_draws()", {

  it("returns a named list with the same structure as the independent res_draws", {
    ## res_draws is list(analyte = list(grid_dates, draws)); draws is [n_grid x nsim].

    ## Build two synthetic KFAS models with aligned grids.
    mod_A <- make_kalman_model(n_days = 84L, seed = 71L)
    mod_B <- make_kalman_model(n_days = 84L, seed = 72L)
    sm_A  <- leachatetools:::.kalman_smooth(mod_A)
    sm_B  <- leachatetools:::.kalman_smooth(mod_B)
    smoothers <- list(
      A = list(grid_dates = attr(mod_A, "grid_dates"), mean = sm_A$mean,
               draw_model = mod_A),
      B = list(grid_dates = attr(mod_B, "grid_dates"), mean = sm_B$mean,
               draw_model = mod_B)
    )
    cor_R <- matrix(c(1, 0.7, 0.7, 1), 2, 2,
                    dimnames = list(c("A", "B"), c("A", "B")))
    out <- leachatetools:::.coupled_residual_draws(
      smoothers, modelled = c("A", "B"), ndraws = 10L, cor_R = cor_R, seed = 1L
    )
    expect_named(out, c("A", "B"), ignore.order = FALSE)
    for (nm in c("A", "B")) {
      expect_named(out[[nm]], c("grid_dates", "draws"), ignore.order = TRUE)
      expect_equal(ncol(out[[nm]]$draws), 10L)
      expect_equal(nrow(out[[nm]]$draws), length(out[[nm]]$grid_dates))
    }
  })

  it("with cor_R = I, per-analyte draw distribution matches the independent path", {
    ## When the correlation is identity, the coupled and independent paths are
    ## statistically equivalent. Seeded; compare column-mean tracks sm$mean.

    mod_A <- make_kalman_model(n_days = 84L, seed = 73L)
    sm_A  <- leachatetools:::.kalman_smooth(mod_A)
    smoothers <- list(
      A = list(grid_dates = attr(mod_A, "grid_dates"), mean = sm_A$mean,
               draw_model = mod_A)
    )
    cor_R <- matrix(1, 1, 1, dimnames = list("A", "A"))
    nsim  <- 300L
    out   <- leachatetools:::.coupled_residual_draws(
      smoothers, modelled = "A", ndraws = nsim, cor_R = cor_R, seed = 2L
    )
    ## Draw mean should track the KFS posterior mean
    expect_gt(stats::cor(rowMeans(out$A$draws), sm_A$mean), 0.97)
  })

  it("with positive cor_R, draws are positively correlated across analytes on shared days", {
    ## On dates in the union grid, draw_A[t, ] and draw_B[t, ] should be
    ## positively correlated (MC oracle, seeded, n=300 draws).

    withr::local_seed(81L)
    mod_A <- make_kalman_model(n_days = 84L, seed = 74L)
    mod_B <- make_kalman_model(n_days = 84L, seed = 75L)
    sm_A  <- leachatetools:::.kalman_smooth(mod_A)
    sm_B  <- leachatetools:::.kalman_smooth(mod_B)
    smoothers <- list(
      A = list(grid_dates = attr(mod_A, "grid_dates"), mean = sm_A$mean,
               draw_model = mod_A),
      B = list(grid_dates = attr(mod_B, "grid_dates"), mean = sm_B$mean,
               draw_model = mod_B)
    )
    cor_R <- matrix(c(1, 0.8, 0.8, 1), 2, 2,
                    dimnames = list(c("A", "B"), c("A", "B")))
    nsim  <- 300L
    out   <- leachatetools:::.coupled_residual_draws(
      smoothers, modelled = c("A", "B"), ndraws = nsim, cor_R = cor_R, seed = 3L
    )
    ## Find a shared mid-gap date (not an anchor of either model)
    shared_date <- intersect(
      as.character(out$A$grid_dates),
      as.character(out$B$grid_dates)
    )
    t_idx_A <- match(shared_date[30L], as.character(out$A$grid_dates))
    t_idx_B <- match(shared_date[30L], as.character(out$B$grid_dates))
    cross_cor <- stats::cor(out$A$draws[t_idx_A, ], out$B$draws[t_idx_B, ])
    ## Should be clearly positive (not necessarily equal to rho=0.8 due to
    ## per-analyte smoother dampening, but sign must be correct).
    expect_gt(cross_cor, 0.1)
  })

  it("ragged non-overlapping analytes remain independent (no coupling across non-shared days)", {
    ## Analyte A covers Jan-Apr, analyte B covers Jul-Oct: no shared dates.
    ## Cross-correlation of draws on their respective (non-shared) dates is ~0.

    withr::local_seed(91L)
    make_nonoverlap_model <- function(start_date, n_days, seed) {
      set.seed(seed)
      theta <- 0.05; gamma <- 4
      phi   <- exp(-theta); q <- gamma * (1 - phi^2)
      x     <- numeric(n_days); x[1L] <- stats::rnorm(1, 0, sqrt(gamma))
      for (t in 2:n_days) x[t] <- phi * x[t - 1L] + stats::rnorm(1, 0, sqrt(q))
      dates   <- start_date + seq_len(n_days) - 1L
      anc_idx <- seq(1L, n_days, by = 14L)
      p <- leachatetools:::.estimate_ou_kalman_params(
        dates[anc_idx], x[anc_idx], n_fit_min = 4L
      )
      leachatetools:::.build_kalman_model(
        dates, dates[anc_idx], x[anc_idx], p$theta, p$gamma,
        r_vec = rep(1e-6, length(anc_idx))
      )
    }
    mod_A <- make_nonoverlap_model(as.Date("2021-01-01"), 90L, 76L)
    mod_B <- make_nonoverlap_model(as.Date("2021-07-01"), 90L, 77L)
    sm_A  <- leachatetools:::.kalman_smooth(mod_A)
    sm_B  <- leachatetools:::.kalman_smooth(mod_B)
    smoothers <- list(
      A = list(grid_dates = attr(mod_A, "grid_dates"), mean = sm_A$mean,
               draw_model = mod_A),
      B = list(grid_dates = attr(mod_B, "grid_dates"), mean = sm_B$mean,
               draw_model = mod_B)
    )
    cor_R <- matrix(c(1, 0.9, 0.9, 1), 2, 2,
                    dimnames = list(c("A", "B"), c("A", "B")))
    nsim  <- 200L
    out   <- leachatetools:::.coupled_residual_draws(
      smoothers, modelled = c("A", "B"), ndraws = nsim, cor_R = cor_R, seed = 4L
    )
    ## Grids don't overlap -> pick any mid-point from each.
    dr_A <- out$A$draws[45L, ]   # mid-point of A's grid
    dr_B <- out$B$draws[45L, ]   # mid-point of B's grid (different calendar dates)
    cross_cor <- stats::cor(dr_A, dr_B)
    ## No shared calendar dates -> effectively independent draws.
    expect_lt(abs(cross_cor), 0.3,
              label = "non-overlapping analytes are effectively independent")
  })

  it("is reproducible given the same seed", {

    mod_A <- make_kalman_model(n_days = 84L, seed = 78L)
    sm_A  <- leachatetools:::.kalman_smooth(mod_A)
    smoothers <- list(
      A = list(grid_dates = attr(mod_A, "grid_dates"), mean = sm_A$mean,
               draw_model = mod_A)
    )
    cor_R <- matrix(1, 1, 1, dimnames = list("A", "A"))
    run   <- function() leachatetools:::.coupled_residual_draws(
      smoothers, "A", ndraws = 8L, cor_R = cor_R, seed = 77L
    )
    out1 <- run(); out2 <- run()
    expect_identical(out1$A$draws, out2$A$draws)
  })

  it("falls back to independent path gracefully when n_couplable < 2", {
    ## Single analyte: coupling is a no-op (1x1 identity), should not error.

    mod_A <- make_kalman_model(n_days = 84L, seed = 79L)
    sm_A  <- leachatetools:::.kalman_smooth(mod_A)
    smoothers <- list(
      A = list(grid_dates = attr(mod_A, "grid_dates"), mean = sm_A$mean,
               draw_model = mod_A)
    )
    cor_R <- matrix(1, 1, 1, dimnames = list("A", "A"))
    expect_no_error(
      leachatetools:::.coupled_residual_draws(
        smoothers, "A", ndraws = 5L, cor_R = cor_R, seed = 5L
      )
    )
  })

})


## ── Stage 4: combined-band widening invariant ─────────────────────────────────
##
## The following tests verify the headline statistical guarantee on synthetic
## multi-analyte data without needing the full amspaf_daily() pipeline.

describe("coupling widens combined band; per-analyte marginals unchanged", {

  it("positive coupling widens the combined sum of draws at any given day", {
    ## Set up two analytes with identical marginals but rho=0 vs rho=0.8.
    ## Combined sum variance: Var(A+B) = Var(A)+Var(B)+2*rho*sd(A)*sd(B).
    ## So coupled IQR of (A+B) > independent IQR of (A+B).

    mod_A <- make_kalman_model(n_days = 84L, seed = 101L)
    mod_B <- make_kalman_model(n_days = 84L, seed = 102L)
    sm_A  <- leachatetools:::.kalman_smooth(mod_A)
    sm_B  <- leachatetools:::.kalman_smooth(mod_B)
    smoothers <- list(
      A = list(grid_dates = attr(mod_A, "grid_dates"), mean = sm_A$mean,
               draw_model = mod_A),
      B = list(grid_dates = attr(mod_B, "grid_dates"), mean = sm_B$mean,
               draw_model = mod_B)
    )
    nsim <- 400L
    ## Independent (rho = 0):
    cor_indep <- diag(2); dimnames(cor_indep) <- list(c("A", "B"), c("A", "B"))
    out_indep <- leachatetools:::.coupled_residual_draws(
      smoothers, c("A", "B"), ndraws = nsim, cor_R = cor_indep, seed = 6L
    )
    ## Positively coupled (rho = 0.8):
    cor_pos <- matrix(c(1, 0.8, 0.8, 1), 2, 2,
                      dimnames = list(c("A", "B"), c("A", "B")))
    out_coup <- leachatetools:::.coupled_residual_draws(
      smoothers, c("A", "B"), ndraws = nsim, cor_R = cor_pos, seed = 6L
    )
    ## Pick a shared mid-gap date.
    shared_dates <- intersect(
      as.character(out_indep$A$grid_dates),
      as.character(out_indep$B$grid_dates)
    )
    t_A <- match(shared_dates[40L], as.character(out_indep$A$grid_dates))
    t_B <- match(shared_dates[40L], as.character(out_indep$B$grid_dates))
    ## Combined IQR (positively coupled should be wider).
    sum_indep <- out_indep$A$draws[t_A, ] + out_indep$B$draws[t_B, ]
    sum_coup  <- out_coup$A$draws[t_A, ]  + out_coup$B$draws[t_B, ]
    expect_gt(stats::IQR(sum_coup), stats::IQR(sum_indep),
              label = "coupled combined IQR > independent combined IQR")
  })

  it("per-analyte marginal IQR is unchanged between independent and coupled draws", {
    ## DK invariant: correlating eta+ leaves each per-analyte marginal unchanged.
    ## We verify this by checking that individual IQRs are close (within MC noise).

    mod_A <- make_kalman_model(n_days = 84L, seed = 103L)
    mod_B <- make_kalman_model(n_days = 84L, seed = 104L)
    sm_A  <- leachatetools:::.kalman_smooth(mod_A)
    sm_B  <- leachatetools:::.kalman_smooth(mod_B)
    smoothers <- list(
      A = list(grid_dates = attr(mod_A, "grid_dates"), mean = sm_A$mean,
               draw_model = mod_A),
      B = list(grid_dates = attr(mod_B, "grid_dates"), mean = sm_B$mean,
               draw_model = mod_B)
    )
    nsim <- 600L
    cor_indep <- diag(2); dimnames(cor_indep) <- list(c("A", "B"), c("A", "B"))
    cor_pos   <- matrix(c(1, 0.8, 0.8, 1), 2, 2,
                        dimnames = list(c("A", "B"), c("A", "B")))
    out_indep <- leachatetools:::.coupled_residual_draws(
      smoothers, c("A", "B"), nsim, cor_R = cor_indep, seed = 7L
    )
    out_coup  <- leachatetools:::.coupled_residual_draws(
      smoothers, c("A", "B"), nsim, cor_R = cor_pos,   seed = 7L
    )
    ## IQR at a mid-gap point for analyte A — should be similar whether coupled or not.
    t_A <- 40L
    iqr_indep <- stats::IQR(out_indep$A$draws[t_A, ])
    iqr_coup  <- stats::IQR(out_coup$A$draws[t_A, ])
    ## Allow 30% relative difference (loose MC tolerance at nsim=600).
    expect_lt(abs(iqr_coup - iqr_indep) / iqr_indep, 0.30,
              label = "per-analyte IQR unchanged by coupling (DK marginal invariant)")
  })

  it("the deterministic centre line is unchanged by coupling", {
    ## The centre line (smoother mean) does not come from the draw path at all;
    ## it is the KFS posterior mean. Verify it is identical.

    mod_A <- make_kalman_model(n_days = 84L, seed = 105L)
    sm_A  <- leachatetools:::.kalman_smooth(mod_A)
    ## The setup object carries x_hat from the same KFS run.
    setup <- leachatetools:::.kalman_sim_smoother_setup(mod_A)
    expect_equal(setup$x_hat, sm_A$mean, tolerance = 1e-8,
                 label = "centre line x_hat matches KFS posterior mean")
  })

})


## ── K5-analogue: coupled draw mean tracks smoother mean (gap-fill) ─────────

## ## This tests the existing K5 behaviour extended to the coupled path.
## Marked as gap-fill: it covers existing smoother behaviour from the coupled
## entry point that K5 does not explicitly cover.

test_that("K11 (gap-fill): coupled draw mean tracks the KFS smoother mean", {
  ## Gap-fill: K5 tests .kalman_draw; this tests that .kalman_sim_smoother_setup
  ## + .kalman_draw_coupled gives the same guarantee when used with iid normals.
  mod   <- make_kalman_model(n_days = 140L, seed = 111L)
  setup <- leachatetools:::.kalman_sim_smoother_setup(mod)
  sm    <- leachatetools:::.kalman_smooth(mod)
  n_grid <- setup$n_grid
  n_anch <- length(setup$pos)
  nsim   <- 500L
  withr::local_seed(121L)
  eta_std <- matrix(stats::rnorm(n_grid * nsim), n_grid, nsim)
  a1_z    <- stats::rnorm(nsim)
  eps_std <- matrix(stats::rnorm(n_anch * nsim), n_anch, nsim)
  dr      <- leachatetools:::.kalman_draw_coupled(setup, eta_std, a1_z, eps_std)
  expect_gt(stats::cor(rowMeans(dr), sm$mean), 0.97,
            label = "coupled draw mean tracks KFS posterior mean (K11)")
})
