## Tests for S6 anchor perturbation + S7 co-analyte perturbation
## (issue #16, Chunk D).
##
## Properties tested:
##   E1. fdm$ou[[nm]]$c_norm_obs_anch is non-NULL and positive for model-tier
##   E2. fdm$co_grab_map is a named list with entries for co-analytes
##   E3. .perturb_anchors_in_model changes anchor S when CV > 0
##   E4. .perturb_anchors_in_model leaves anchors unchanged when CV = 0
##   E5. .perturb_co_split returns different co values when CV > 0
##   E6. .perturb_co_split coherence: same-grab days share the same multiplier
##   E7. Anchor + co perturbations change .predict_daily_tox output

library(testthat)
library(leachatetools)


## ── Shared setup (reuses helpers from test-daily-target.R) ───────────────────

make_chem_e <- function(site, dates, mult = 1, seed = 1) {
  set.seed(seed)
  analytes <- c("Cu", "Zn", "Ni", "pH", "DOC", "hardness", "Ca", "Mg")
  purrr::map_dfr(dates, function(d) tibble::tibble(
    sample_id = paste0(site, format(d, "%Y%m%d")),
    site_id   = site,
    datetime  = d,
    analyte   = analytes,
    value     = c(
      exp(stats::rnorm(1, log(0.5), 0.3)) * mult,
      exp(stats::rnorm(1, log(5),   0.4)) * mult,
      exp(stats::rnorm(1, log(0.3), 0.3)) * mult,
      stats::runif(1, 6.5, 8), stats::runif(1, 1, 5), stats::runif(1, 20, 60),
      stats::runif(1, 4, 12),  stats::runif(1, 2, 8)
    ),
    detected = TRUE
  ))
}

make_hydro_e <- function(n = 700, seed = 99) {
  set.seed(seed)
  tibble::tibble(
    date  = seq(as.Date("2020-07-01"), by = "day", length.out = n),
    value = pmax(0, stats::rnorm(n, 2, 4))
  )
}

## One-time setup (shared across all tests in this file)
.te <- local({
  dates <- seq(as.Date("2021-01-01"), by = "2 weeks", length.out = 40)
  hydro <- make_hydro_e()
  ref   <- make_chem_e("reference", dates, seed = 1)
  tgt   <- make_chem_e("target",    dates, mult = 5, seed = 2)
  rm    <- fit_reference_model(ref, hydro = hydro, conc_units = "ug/L",
                                min_obs_model = 10L,
                                api_windows_short = 7L,
                                api_windows_long  = 30L)
  all_dates <- seq(dates[1], dates[length(dates)], by = "day")
  meta      <- leachatetools:::.load_analyte_metadata(NULL)
  tox_nms   <- meta$analyte[!is.na(meta$ssd_available) & meta$ssd_available]
  daily_long <- leachatetools:::.build_daily_chem(
    site_rows     = tgt,
    dates         = all_dates,
    interpolation = "forward_fill",
    leading_edge  = "drop",
    tox_analytes  = tox_nms
  )
  fdm <- leachatetools:::.fit_daily_target(
    site_rows        = tgt,
    reference_model  = rm,
    imputation_model = NULL,
    conc_units       = "ug/L",
    meta             = meta,
    tox_analytes     = tox_nms,
    daily_long       = daily_long
  )
  list(rm = rm, tgt = tgt, all_dates = all_dates, daily_long = daily_long,
       meta = meta, tox_nms = tox_nms, fdm = fdm)
})


## ── E1: smoother scaffold carries a grid/mean and a draw model ───────────────

test_that("E1: fdm$smoothers carry grid/mean and a draw model for >= 1 analyte", {
  fdm <- .te$fdm
  has_dm <- vapply(fdm$modelled, function(nm) {
    sm <- fdm$smoothers[[nm]]
    !is.null(sm) && length(sm$grid_dates) > 0L &&
      length(sm$mean) == length(sm$grid_dates) && !is.null(sm$draw_model)
  }, logical(1L))
  expect_true(any(has_dm), label = "at least one analyte has a draw model")
})


## ── E2: co_grab_map structure ─────────────────────────────────────────────────

test_that("E2: fdm$co_grab_map is a named list covering co-analytes", {
  fdm <- .te$fdm
  expect_false(is.null(fdm$co_grab_map))
  expect_type(fdm$co_grab_map, "list")
  expect_true(length(fdm$co_grab_map) > 0L)
  # Each element should be a named character vector (names = date strings)
  first_co <- fdm$co_grab_map[[1L]]
  expect_type(first_co, "character")
  expect_true(length(first_co) == length(fdm$qdates))
  # At least some entries should be non-NA (grabs exist within date range)
  expect_true(any(!is.na(first_co)),
              label = "at least some dates have a source grab")
})


## ── E3: grab_cv > 0 widens the draw spread at anchors (S6 as obs noise) ───────

test_that("E3: grab_cv > 0 widens draw spread at anchors vs no grab_cv", {
  fdm0 <- .te$fdm                                       # built without grab_cv
  fdm_cv <- leachatetools:::.fit_daily_target(
    site_rows = .te$tgt, reference_model = .te$rm, imputation_model = NULL,
    conc_units = "ug/L", meta = .te$meta, tox_analytes = .te$tox_nms,
    daily_long = .te$daily_long, grab_cv = 0.2)

  nm_test <- NULL
  for (nm in fdm0$modelled) {
    if (!is.null(fdm0$smoothers[[nm]]$draw_model) &&
        !is.null(fdm_cv$smoothers[[nm]]$draw_model)) { nm_test <- nm; break }
  }
  skip_if(is.null(nm_test), "No analyte with a draw model in both fits")

  grid <- fdm0$smoothers[[nm_test]]$grid_dates
  pos  <- match(fdm0$tm$models[[nm_test]]$anchors$date, grid)
  pos  <- pos[!is.na(pos)]
  set.seed(1L); d0  <- leachatetools:::.kalman_draw(fdm0$smoothers[[nm_test]]$draw_model, 200L)
  set.seed(1L); dcv <- leachatetools:::.kalman_draw(fdm_cv$smoothers[[nm_test]]$draw_model, 200L)
  sd0  <- mean(apply(d0[pos, , drop = FALSE],  1L, stats::sd))
  sdcv <- mean(apply(dcv[pos, , drop = FALSE], 1L, stats::sd))
  expect_gt(sdcv, sd0)                                  # S6 adds anchor spread
})


## ── E4: without grab_cv the draws pin (tighter) at anchors than mid-gap ───────

test_that("E4: with no grab_cv, draw spread is smaller at anchors than mid-gap", {
  fdm <- .te$fdm
  nm  <- NULL
  for (a in fdm$modelled) if (!is.null(fdm$smoothers[[a]]$draw_model)) { nm <- a; break }
  skip_if(is.null(nm), "No analyte with a draw model")

  grid <- fdm$smoothers[[nm]]$grid_dates
  pos  <- match(fdm$tm$models[[nm]]$anchors$date, grid); pos <- pos[!is.na(pos)]
  skip_if(length(pos) < 2L, "Need >= 2 anchors on the grid")
  mid  <- pmin(pos[-length(pos)] + 7L, length(grid))     # mid-gap samples
  set.seed(2L); d <- leachatetools:::.kalman_draw(fdm$smoothers[[nm]]$draw_model, 200L)
  sd_anchor <- mean(apply(d[pos, , drop = FALSE], 1L, stats::sd))
  sd_mid    <- mean(apply(d[mid, , drop = FALSE], 1L, stats::sd))
  expect_lt(sd_anchor, sd_mid)
})


## ── E5: .perturb_co_split returns different values when CV > 0 ───────────────

test_that("E5: co_split values differ after perturbation with CV > 0", {
  fdm <- .te$fdm
  set.seed(7L)
  res <- leachatetools:::.perturb_co_split(fdm, grab_cv_co = 0.15)
  co_orig  <- fdm$co_split
  co_pert  <- res$co_split

  # At least some values should differ
  any_diff <- FALSE
  for (d in intersect(names(co_orig), names(co_pert))) {
    if (!isTRUE(all.equal(co_orig[[d]]$value, co_pert[[d]]$value))) {
      any_diff <- TRUE; break
    }
  }
  expect_true(any_diff, label = "at least one date has different co-analyte values")
})


## ── E6: coherence — same-grab days share the same multiplier ────────────────

test_that("E6: forward-filled days from the same grab get the same multiplier", {
  fdm <- .te$fdm

  # Find a co-analyte with multiple consecutive forward-filled days from one grab
  # (i.e., a gap between two consecutive grabs)
  a <- names(fdm$co_grab_map)[1L]
  src_map <- fdm$co_grab_map[[a]]
  src_tab <- table(src_map[!is.na(src_map)])
  # Pick the most common source grab (likely has many forward-filled days)
  common_src <- names(which.max(src_tab))
  day_strs   <- names(src_map)[!is.na(src_map) & src_map == common_src]
  skip_if(length(day_strs) < 2L, "Need >= 2 days sharing a grab for coherence test")

  set.seed(3L)
  res <- leachatetools:::.perturb_co_split(fdm, grab_cv_co = 0.15)

  # All days sharing the same grab should have the same ratio to original
  orig_vals <- vapply(day_strs, function(d) {
    cd <- fdm$co_split[[d]]
    cd$value[cd$analyte == a]
  }, numeric(1L))
  pert_vals <- vapply(day_strs, function(d) {
    cd <- res$co_split[[d]]
    cd$value[cd$analyte == a]
  }, numeric(1L))
  ratios <- unname(pert_vals / orig_vals)
  # All ratios should be the same (one multiplier per grab)
  expect_equal(ratios, rep(ratios[1L], length(ratios)), tolerance = 1e-10,
               label = "all same-grab days have identical multiplier")
})


## ── E7: a draw (S4+S6 residual path + GAM) + S7 shift .predict_daily_tox ─────

test_that("E7: residual draw (S4/S6) + GAM + co perturbations shift the output", {
  fdm <- leachatetools:::.fit_daily_target(
    site_rows = .te$tgt, reference_model = .te$rm, imputation_model = NULL,
    conc_units = "ug/L", meta = .te$meta, tox_analytes = .te$tox_nms,
    daily_long = .te$daily_long, grab_cv = 0.2)

  rows_base <- leachatetools:::.predict_daily_tox(fdm)   # centre

  set.seed(11L)
  tm_p <- leachatetools:::.perturb_target_model(fdm$tm)
  rp <- stats::setNames(lapply(fdm$modelled, function(nm) {
    sm <- fdm$smoothers[[nm]]
    if (is.null(sm$draw_model))
      leachatetools:::.residual_on_qdates(sm$grid_dates, sm$mean, fdm$qdates)
    else {
      dr <- leachatetools:::.kalman_draw(sm$draw_model, 1L)
      leachatetools:::.residual_on_qdates(sm$grid_dates, dr[, 1L], fdm$qdates)
    }
  }), fdm$modelled)
  co_p <- leachatetools:::.perturb_co_split(fdm, grab_cv_co = 0.2)

  rows_pert <- leachatetools:::.predict_daily_tox(
    fdm, tm_p = tm_p, residual_paths = rp,
    co_split = fdm$co_split, wq_long = co_p$wq_long)

  common_analytes <- intersect(unique(rows_base$analyte), unique(rows_pert$analyte))
  any_diff <- FALSE
  for (nm in common_analytes) {
    v_base <- rows_base$value[rows_base$analyte == nm]
    v_pert <- rows_pert$value[rows_pert$analyte == nm]
    if (length(v_base) == length(v_pert) && !isTRUE(all.equal(v_base, v_pert))) {
      any_diff <- TRUE; break
    }
  }
  expect_true(any_diff, label = "values differ after S4/S6 + GAM + S7 perturbation")
})
