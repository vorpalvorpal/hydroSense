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


## ── E1: c_norm_obs_anch precomputed for model-tier ───────────────────────────

test_that("E1: fdm$ou[[nm]]$c_norm_obs_anch is non-NULL and positive for >= 1 analyte", {
  fdm <- .te$fdm
  has_obs <- vapply(fdm$modelled, function(nm) {
    obs <- fdm$ou[[nm]]$c_norm_obs_anch
    !is.null(obs) && length(obs) > 0L && all(is.finite(obs)) && all(obs >= 0)
  }, logical(1L))
  expect_true(any(has_obs), label = "at least one analyte has valid c_norm_obs_anch")
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


## ── E3: .perturb_anchors_in_model changes S when CV > 0 ──────────────────────

test_that("E3: anchor S values change when grab_cv > 0", {
  fdm <- .te$fdm
  tm  <- fdm$tm
  set.seed(42L)
  tm_p <- leachatetools:::.perturb_anchors_in_model(tm, fdm, grab_cv = 0.1)

  # Find an analyte with a valid c_norm_obs_anch
  nm_test <- NULL
  for (nm in fdm$modelled) {
    if (!is.null(fdm$ou[[nm]]$c_norm_obs_anch) && nrow(tm$models[[nm]]$anchors) >= 2L) {
      nm_test <- nm; break
    }
  }
  skip_if(is.null(nm_test), "No suitable analyte with valid c_norm_obs_anch")

  orig_S <- tm$models[[nm_test]]$anchors$S
  pert_S <- tm_p$models[[nm_test]]$anchors$S
  expect_false(isTRUE(all.equal(orig_S, pert_S)),
               label = "perturbed S differs from original")
})


## ── E4: zero CV leaves anchors unchanged ────────────────────────────────────

test_that("E4: grab_cv = 0 leaves anchors identical to original", {
  fdm <- .te$fdm
  tm  <- fdm$tm
  tm_p <- leachatetools:::.perturb_anchors_in_model(tm, fdm, grab_cv = 0)
  for (nm in fdm$modelled) {
    expect_equal(tm_p$models[[nm]]$anchors$S, tm$models[[nm]]$anchors$S,
                 info = nm)
  }
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


## ── E7: S6+S7 perturbations change .predict_daily_tox output ────────────────

test_that("E7: anchor + co perturbations shift .predict_daily_tox output", {
  fdm <- .te$fdm

  # Point-mode baseline
  rows_base <- leachatetools:::.predict_daily_tox(fdm, eps_paths = NULL)

  # Perturb anchors (S6) and co-analytes (S7)
  set.seed(11L)
  tm_p <- leachatetools:::.perturb_anchors_in_model(fdm$tm, fdm, grab_cv = 0.2)
  co_p <- leachatetools:::.perturb_co_split(fdm, grab_cv_co = 0.2)

  rows_pert <- leachatetools:::.predict_daily_tox(
    fdm,
    tm_p     = tm_p,
    co_split = co_p$co_split,
    wq_long  = co_p$wq_long
  )

  # At least some values should differ after perturbation
  common_analytes <- intersect(
    unique(rows_base$analyte), unique(rows_pert$analyte)
  )
  any_diff <- FALSE
  for (nm in common_analytes) {
    v_base <- rows_base$value[rows_base$analyte == nm]
    v_pert <- rows_pert$value[rows_pert$analyte == nm]
    if (!isTRUE(all.equal(v_base, v_pert))) { any_diff <- TRUE; break }
  }
  expect_true(any_diff,
              label = "at least one analyte has different values after S6+S7 perturbation")
})
