## Tests for GAM perturbation helpers (issue #16, Chunk B).
## Stan-free: all GAMs fitted with mgcv on synthetic data.
##
## Properties tested:
##   C1. .perturb_gam: NULL input → NULL output
##   C2. .perturb_gam: GAM without $Vp → returned unchanged
##   C3. .perturb_gam: perturbed coefficients differ from original
##   C4. .perturb_gam: mean prediction over many draws ≈ original prediction
##   C5. .perturb_gam: original object is not mutated
##   C6. .perturb_target_model: pooled analytes get the same perturbed fit
##   C7. .perturb_target_model: non-pooled analytes get independent fits
##   C8. .perturb_target_model: original model not mutated
##   C9. .perturb_reference_model: each gamm_fit has perturbed coefficients
##   C10. .perturb_target_model(perturb_reference=TRUE) perturbs ref model too

library(testthat)
library(leachatetools)

## ── helpers ──────────────────────────────────────────────────────────────────

make_simple_gam <- function(n = 60, seed = 1L) {
  set.seed(seed)
  x   <- seq(0, 1, length.out = n)
  y   <- sin(2 * pi * x) + stats::rnorm(n, 0, 0.2)
  df  <- data.frame(x = x, y = y)
  mgcv::gam(y ~ s(x), data = df, method = "REML")
}

# Build a mock target_model: two pooled analytes (Cu, Zn) + one bridge (Ni)
make_mock_target_model <- function(gam_fit) {
  pool_fit <- gam_fit                      # shared pooled impact_fit
  models <- list(
    Cu = list(
      impact_fit   = pool_fit,
      pooled       = TRUE,
      analyte      = "Cu",
      pool_levels  = c("Cu", "Zn"),
      pool_center  = 0,
      pool_scale   = 1,
      wq_fit       = NULL,
      tier         = "model",
      n_obs        = 20L,
      anchors      = data.frame(date = Sys.Date(), S = 0)
    ),
    Zn = list(
      impact_fit   = pool_fit,
      pooled       = TRUE,
      analyte      = "Zn",
      pool_levels  = c("Cu", "Zn"),
      pool_center  = 0,
      pool_scale   = 1,
      wq_fit       = NULL,
      tier         = "model",
      n_obs        = 20L,
      anchors      = data.frame(date = Sys.Date(), S = 0)
    ),
    Ni = list(
      impact_fit   = NULL,
      pooled       = FALSE,
      wq_fit       = NULL,
      tier         = "bridge",
      n_obs        = 5L,
      anchors      = data.frame(date = Sys.Date(), S = 0)
    )
  )
  structure(
    list(
      models          = models,
      reference_model = NULL,
      pca             = NULL,
      pc_cols         = character(0),
      hydro           = data.frame(date = Sys.Date(), value = 1),
      hydro_type      = "rainfall",
      fit_date        = Sys.Date()
    ),
    class = "target_model"
  )
}

make_mock_reference_model <- function(gam_fit) {
  models <- list(
    Cu = list(gamm_fit = gam_fit, tier = "model", n_obs = 25L,
              static_ref = 0, tau_short = 7, tau_long = 90,
              best_aic = 0, null_aic = 10,
              obs = data.frame()),
    pH = list(gamm_fit = NULL, tier = "static", n_obs = 5L,
              static_ref = 7.5, tau_short = 7, tau_long = 90,
              best_aic = NA, null_aic = NA,
              obs = data.frame())
  )
  structure(
    list(
      models            = models,
      hydro             = data.frame(date = Sys.Date(), value = 1),
      hydro_type        = "rainfall",
      match_window_days = 3L,
      match_hydro_tol   = Inf,
      static_ref        = 0,
      fit_date          = Sys.Date(),
      summary           = NULL
    ),
    class = "reference_model"
  )
}


## ── C1–C5: .perturb_gam ──────────────────────────────────────────────────────

test_that("C1: NULL input returns NULL", {
  expect_null(leachatetools:::.perturb_gam(NULL))
})

test_that("C2: GAM without $Vp returned unchanged", {
  g     <- make_simple_gam()
  g_nov <- g
  g_nov$Vp <- NULL
  result <- leachatetools:::.perturb_gam(g_nov)
  expect_identical(result$coefficients, g_nov$coefficients)
})

test_that("C3: perturbed coefficients differ from original", {
  g   <- make_simple_gam()
  set.seed(1L)
  g_p <- leachatetools:::.perturb_gam(g)
  expect_false(isTRUE(all.equal(g$coefficients, g_p$coefficients)))
})

test_that("C4: mean prediction over many perturbation draws ≈ original", {
  g    <- make_simple_gam()
  nd   <- data.frame(x = c(0.1, 0.5, 0.9))
  pred_orig <- as.numeric(stats::predict(g, newdata = nd))

  # Draw 2000 perturbed predictions
  set.seed(42L)
  preds <- replicate(2000L, {
    g_p <- leachatetools:::.perturb_gam(g)
    as.numeric(stats::predict(g_p, newdata = nd))
  })
  pred_mean <- rowMeans(preds)

  # Mean over draws should be close to original (tolerance ~3 SE from truth)
  expect_equal(pred_mean, pred_orig, tolerance = 0.05)
})

test_that("C5: .perturb_gam does not mutate the original object", {
  g    <- make_simple_gam()
  orig <- g$coefficients
  set.seed(1L)
  leachatetools:::.perturb_gam(g)
  expect_identical(g$coefficients, orig)
})


## ── C6–C8: .perturb_target_model ─────────────────────────────────────────────

test_that("C6: pooled analytes (Cu, Zn) get the SAME perturbed impact_fit", {
  g   <- make_simple_gam()
  tm  <- make_mock_target_model(g)
  set.seed(3L)
  tm_p <- leachatetools:::.perturb_target_model(tm)
  # Coefficients must be identical between Cu and Zn (same draw)
  expect_identical(
    tm_p$models$Cu$impact_fit$coefficients,
    tm_p$models$Zn$impact_fit$coefficients
  )
})

test_that("C6: pooled analytes' perturbed coefficients differ from original", {
  g   <- make_simple_gam()
  tm  <- make_mock_target_model(g)
  set.seed(3L)
  tm_p <- leachatetools:::.perturb_target_model(tm)
  expect_false(isTRUE(all.equal(
    tm_p$models$Cu$impact_fit$coefficients,
    g$coefficients
  )))
})

test_that("C7: bridge-tier analyte (Ni) impact_fit stays NULL after perturbation", {
  g   <- make_simple_gam()
  tm  <- make_mock_target_model(g)
  tm_p <- leachatetools:::.perturb_target_model(tm)
  expect_null(tm_p$models$Ni$impact_fit)
})

test_that("C8: .perturb_target_model does not mutate the original model", {
  g    <- make_simple_gam()
  tm   <- make_mock_target_model(g)
  orig_coef <- tm$models$Cu$impact_fit$coefficients
  set.seed(1L)
  leachatetools:::.perturb_target_model(tm)
  expect_identical(tm$models$Cu$impact_fit$coefficients, orig_coef)
})


## ── C9–C10: .perturb_reference_model / perturb_reference flag ────────────────

test_that("C9: .perturb_reference_model perturbs gamm_fit coefficients", {
  g   <- make_simple_gam()
  rm  <- make_mock_reference_model(g)
  set.seed(5L)
  rm_p <- leachatetools:::.perturb_reference_model(rm)
  expect_false(isTRUE(all.equal(
    rm_p$models$Cu$gamm_fit$coefficients,
    g$coefficients
  )))
  # Static-tier analyte (pH) has NULL gamm_fit — should stay NULL
  expect_null(rm_p$models$pH$gamm_fit)
})

test_that("C10: perturb_reference=TRUE perturbs embedded reference model", {
  g   <- make_simple_gam()
  tm  <- make_mock_target_model(g)
  # Attach a mock reference model
  rm  <- make_mock_reference_model(g)
  tm$reference_model <- rm
  orig_coef <- rm$models$Cu$gamm_fit$coefficients

  set.seed(6L)
  tm_p <- leachatetools:::.perturb_target_model(tm, perturb_reference = TRUE)
  expect_false(isTRUE(all.equal(
    tm_p$reference_model$models$Cu$gamm_fit$coefficients,
    orig_coef
  )))
})

test_that("C10: perturb_reference=FALSE leaves embedded reference unchanged", {
  g   <- make_simple_gam()
  tm  <- make_mock_target_model(g)
  rm  <- make_mock_reference_model(g)
  tm$reference_model <- rm
  orig_coef <- rm$models$Cu$gamm_fit$coefficients

  tm_p <- leachatetools:::.perturb_target_model(tm, perturb_reference = FALSE)
  expect_identical(
    tm_p$reference_model$models$Cu$gamm_fit$coefficients,
    orig_coef
  )
})
