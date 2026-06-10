## Tests for the season-blind target impact model (issue #14).
##
## All Stan-free: the imputation (tier-2) front-end is gated and not exercised
## here. Covers fit_target_model(), the season-blind property, the
## perfect-management invariant, the elevated-site contrast, the residual
## bridge interpolation, and the amspaf_daily(interpolation = "model") path.

library(testthat)
library(leachatetools)


## ── Helpers ───────────────────────────────────────────────────────────────────

make_chem <- function(site, dates, mult = 1, seed = 1) {
  set.seed(seed)
  analytes <- c("Cu", "Zn", "Ni", "pH", "DOC", "hardness", "Ca", "Mg")
  purrr::map_dfr(dates, function(d) tibble::tibble(
    sample_id = paste0(site, format(d, "%Y%m%d")),
    site_id   = site,
    datetime  = d,
    analyte   = analytes,
    value     = c(
      exp(rnorm(1, log(0.5), 0.3)) * mult,
      exp(rnorm(1, log(5),   0.4)) * mult,
      exp(rnorm(1, log(0.3), 0.3)) * mult,
      runif(1, 6.5, 8), runif(1, 1, 5), runif(1, 20, 60),
      runif(1, 4, 12), runif(1, 2, 8)
    ),
    detected = TRUE
  ))
}

make_hydro <- function(n = 700, seed = 99) {
  set.seed(seed)
  tibble::tibble(
    date  = seq(as.Date("2020-07-01"), by = "day", length.out = n),
    value = pmax(0, rnorm(n, 2, 4))
  )
}

fit_rm <- function(ref, hydro) {
  fit_reference_model(
    ref, hydro = hydro, conc_units = "ug/L", min_obs_model = 10L,
    api_windows_short = c(7L), api_windows_long = c(30L)
  )
}


## ── fit_target_model() structure ──────────────────────────────────────────────

test_that("fit_target_model() returns a target_model with per-analyte models", {
  dates <- seq(as.Date("2021-01-01"), by = "2 weeks", length.out = 40)
  ref   <- make_chem("reference", dates)
  hydro <- make_hydro()
  rm    <- fit_rm(ref, hydro)
  tgt   <- make_chem("target", dates, mult = 5, seed = 2)

  tm <- fit_target_model(tgt, rm, conc_units = "ug/L", min_obs_model = 8L,
                         api_windows_short = c(7L), api_windows_long = c(30L))

  expect_s3_class(tm, "target_model")
  expect_true(length(tm$models) > 0L)
  expect_s3_class(tm$reference_model, "reference_model")
  # every analyte model has the expected slots
  m1 <- tm$models[[1L]]
  expect_true(all(c("impact_fit", "window_short", "window_long", "tier",
                    "n_obs", "anchors") %in% names(m1)))
  expect_true(m1$tier %in% c("model", "bridge"))
  expect_true(all(c("date", "I", "S") %in% names(m1$anchors)))
})

test_that("fit_target_model() errors without a reference_model", {
  dates <- seq(as.Date("2021-01-01"), by = "2 weeks", length.out = 10)
  tgt   <- make_chem("target", dates)
  expect_error(
    fit_target_model(tgt, reference_model = list(), conc_units = "ug/L"),
    regexp = "reference_model"
  )
})

test_that("print.target_model() runs and reports the season-blind tiers", {
  dates <- seq(as.Date("2021-01-01"), by = "2 weeks", length.out = 40)
  rm    <- fit_rm(make_chem("reference", dates), make_hydro())
  tm    <- fit_target_model(make_chem("target", dates, mult = 5, seed = 2),
                            rm, conc_units = "ug/L", min_obs_model = 8L,
                            api_windows_short = c(7L), api_windows_long = c(30L))
  expect_output(print(tm), regexp = "target_model")
  expect_output(print(tm), regexp = "season-blind")
})


## ── Perfect-management invariant & elevated contrast ──────────────────────────

test_that("clean site -> impact centred near zero; elevated site -> positive", {
  dates <- seq(as.Date("2021-01-01"), by = "2 weeks", length.out = 40)
  ref   <- make_chem("reference", dates)
  hydro <- make_hydro()
  rm    <- fit_rm(ref, hydro)
  q     <- tibble::tibble(date = seq(as.Date("2021-02-01"),
                                     as.Date("2021-11-15"), by = "month"))

  tm_clean <- fit_target_model(make_chem("target", dates, mult = 1, seed = 2),
                               rm, conc_units = "ug/L", min_obs_model = 8L,
                               api_windows_short = c(7L), api_windows_long = c(30L))
  tm_hot   <- fit_target_model(make_chem("target", dates, mult = 20, seed = 2),
                               rm, conc_units = "ug/L", min_obs_model = 8L,
                               api_windows_short = c(7L), api_windows_long = c(30L))

  cu_clean <- leachatetools:::.resolve_target_impact(tm_clean, q)
  cu_hot   <- leachatetools:::.resolve_target_impact(tm_hot,   q)
  cu_clean <- cu_clean[cu_clean$analyte == "Cu", ]
  cu_hot   <- cu_hot[cu_hot$analyte == "Cu", ]

  # Clean: centred near zero (small log-normal mean-median offset allowed)
  expect_lt(abs(median(cu_clean$impact)), 0.25)
  # Elevated: systematically and substantially higher than clean
  expect_gt(median(cu_hot$impact), median(cu_clean$impact) + 0.5)
})

test_that("the impact model carries no day-of-year term (season-blind)", {
  # The fitted GAM formula must not reference doy / season.
  dates <- seq(as.Date("2021-01-01"), by = "2 weeks", length.out = 40)
  rm    <- fit_rm(make_chem("reference", dates), make_hydro())
  tm    <- fit_target_model(make_chem("target", dates, mult = 10, seed = 2),
                            rm, conc_units = "ug/L", min_obs_model = 8L,
                            api_windows_short = c(7L), api_windows_long = c(30L))
  for (m in tm$models) {
    if (m$tier == "model" && !is.null(m$impact_fit)) {
      terms_chr <- as.character(stats::formula(m$impact_fit))
      expect_false(any(grepl("doy", terms_chr)))
    }
  }
  succeed()
})


## (Residual interpolation unit tests removed: `.interp_residual()` is replaced
## by the state-space smoother. Pinch-at-anchor / mid-gap behaviour is covered by
## tests/testthat/test-kalman-bridge.R; the residual mean is now a pure temporal
## smooth — hydrology modulates the *variance*, not the residual mean.)


## ── .resolve_target_impact() output shape ─────────────────────────────────────

test_that(".resolve_target_impact() returns C_norm = max(ref_norm + impact, 0)", {
  dates <- seq(as.Date("2021-01-01"), by = "2 weeks", length.out = 40)
  rm    <- fit_rm(make_chem("reference", dates), make_hydro())
  tm    <- fit_target_model(make_chem("target", dates, mult = 10, seed = 2),
                            rm, conc_units = "ug/L", min_obs_model = 8L,
                            api_windows_short = c(7L), api_windows_long = c(30L))
  res <- leachatetools:::.resolve_target_impact(
    tm, tibble::tibble(date = seq(as.Date("2021-02-01"),
                                  as.Date("2021-06-01"), by = "week"))
  )
  expect_named(res, c("date", "analyte", "ref_norm", "impact", "C_norm",
                      "impact_tier"), ignore.order = TRUE)
  expect_true(all(res$C_norm >= 0))
  expect_equal(res$C_norm, pmax(res$ref_norm + res$impact, 0), tolerance = 1e-9)
})


## ── amspaf_daily(interpolation = "model") end-to-end ──────────────────────────

test_that("amspaf_daily(interpolation='model') requires a reference_model", {
  dates <- seq(as.Date("2021-01-01"), by = "2 weeks", length.out = 10)
  tgt   <- make_chem("target", dates)
  expect_error(
    amspaf_daily(tgt, interpolation = "model", conc_units = "ug/L",
                 require_temperature = FALSE),
    regexp = "reference_model"
  )
})

## ── WQ layer + residual d (issue #14, item B) ────────────────────────────────

# Chemistry where Cu tracks EC (a leachate tracer), so the WQ layer has signal.
make_wq_chem <- function(site, dates, mult = 1, seed = 1) {
  set.seed(seed)
  analytes <- c("Cu", "Zn", "Ni", "pH", "EC", "DOC", "hardness", "Ca", "Mg")
  purrr::map_dfr(dates, function(d) {
    ec <- runif(1, 100, 600)
    tibble::tibble(
      sample_id = paste0(site, format(d, "%Y%m%d")), site_id = site,
      datetime = d, analyte = analytes,
      value = c(exp(rnorm(1, log(0.2 + ec / 1000), 0.2)) * mult,
                exp(rnorm(1, log(5), 0.4)) * mult, exp(rnorm(1, log(0.3), 0.3)) * mult,
                runif(1, 6.5, 8), ec, runif(1, 1, 5),
                runif(1, 20, 60), runif(1, 4, 12), runif(1, 2, 8)),
      detected = TRUE)
  })
}

# A Stan-free PCA-only imputation_model (same shape the bs01 script builds).
make_pca_model <- function(chem, wq = c("pH", "EC", "DOC", "Ca", "Mg", "hardness")) {
  pca <- leachatetools:::.prepare_chem_pca(chem, wq_vars = wq)
  structure(list(pca = pca, pca_vars = wq), class = "imputation_model")
}

test_that("fit_target_model() fits a WQ layer when a PCA model is supplied", {
  dates <- seq(as.Date("2021-01-01"), by = "week", length.out = 60)
  rm    <- fit_rm(make_wq_chem("reference", dates), make_hydro())
  tgt   <- make_wq_chem("target", dates, mult = 8, seed = 2)
  tm <- fit_target_model(tgt, rm, imputation_model = make_pca_model(tgt),
                         conc_units = "ug/L", min_obs_model = 10L,
                         api_windows_short = c(7L), api_windows_long = c(30L))
  expect_false(is.null(tm$pca))
  expect_true(length(tm$pc_cols) > 0L)
  # At least one analyte earned a WQ layer
  expect_true(any(vapply(tm$models, function(m) !is.null(m$wq_fit), logical(1L))))
})

test_that("WQ-only PCA model does not trigger the (brms) impute-first warning", {
  dates <- seq(as.Date("2021-01-01"), by = "week", length.out = 60)
  rm    <- fit_rm(make_wq_chem("reference", dates), make_hydro())
  tgt   <- make_wq_chem("target", dates, mult = 8, seed = 2)
  expect_no_warning(
    fit_target_model(tgt, rm, imputation_model = make_pca_model(tgt),
                     conc_units = "ug/L", min_obs_model = 10L,
                     api_windows_short = c(7L), api_windows_long = c(30L))
  )
})

test_that(".resolve_target_impact() uses the 'wq' tier when wq is supplied", {
  dates <- seq(as.Date("2021-01-01"), by = "week", length.out = 60)
  rm    <- fit_rm(make_wq_chem("reference", dates), make_hydro())
  tgt   <- make_wq_chem("target", dates, mult = 8, seed = 2)
  tm <- fit_target_model(tgt, rm, imputation_model = make_pca_model(tgt),
                         conc_units = "ug/L", min_obs_model = 10L,
                         api_windows_short = c(7L), api_windows_long = c(30L))

  qdates <- as.Date(c("2021-04-07", "2021-04-14"))
  wq <- tibble::tibble(
    sample_id = rep(as.character(qdates), each = 2),
    analyte   = rep(c("EC", "pH"), times = 2),
    value     = c(550, 7.2, 120, 7.2)     # day 1 high EC, day 2 low EC
  )
  res <- leachatetools:::.resolve_target_impact(
    tm, tibble::tibble(date = qdates), analytes = "Cu", wq = wq
  )
  expect_true(all(res$impact_tier == "wq"))
  # Cu tracks EC: the high-EC day should predict a higher normalised Cu
  cu_hi <- res$C_norm[res$date == qdates[1]]
  cu_lo <- res$C_norm[res$date == qdates[2]]
  expect_gt(cu_hi, cu_lo)
})

test_that("without a wq argument the resolver falls back to the impact tiers", {
  dates <- seq(as.Date("2021-01-01"), by = "week", length.out = 60)
  rm    <- fit_rm(make_wq_chem("reference", dates), make_hydro())
  tgt   <- make_wq_chem("target", dates, mult = 8, seed = 2)
  tm <- fit_target_model(tgt, rm, imputation_model = make_pca_model(tgt),
                         conc_units = "ug/L", min_obs_model = 10L,
                         api_windows_short = c(7L), api_windows_long = c(30L))
  res <- leachatetools:::.resolve_target_impact(
    tm, tibble::tibble(date = as.Date("2021-04-07")), analytes = "Cu"  # no wq
  )
  expect_true(all(res$impact_tier %in% c("model", "bridge")))
})


## ── hierarchical pooling of the hydro response (issue #14, item C) ────────────

# Chemistry whose metal impact is genuinely driven by short-window antecedent
# rainfall (first flush) — so the pooled factor-smooth response has signal.
make_hydro_driven <- function(site, dates, hydro, mult = 1, seed = 5) {
  set.seed(seed)
  analytes <- c("Cu", "Zn", "Ni", "pH", "DOC", "hardness", "Ca", "Mg")
  api7 <- leachatetools:::.compute_api(hydro$value, hydro$date, dates, 7L)
  api7 <- (api7 - mean(api7)) / (stats::sd(api7) + 1e-9)
  purrr::pmap_dfr(list(dates, api7), function(d, a) {
    flush <- exp(0.8 * a)   # metals rise with antecedent rainfall
    tibble::tibble(
      sample_id = paste0(site, format(d, "%Y%m%d")), site_id = site,
      datetime = d, analyte = analytes,
      value = c(exp(rnorm(1, log(0.5), 0.15)) * mult * flush,
                exp(rnorm(1, log(5),   0.15)) * mult * flush,
                exp(rnorm(1, log(0.3), 0.15)) * mult * flush,
                runif(1, 6.5, 8), runif(1, 1, 5), runif(1, 20, 60),
                runif(1, 4, 12), runif(1, 2, 8)),
      detected = TRUE)
  })
}

test_that("pool = TRUE produces a valid model and finite predictions", {
  dates <- seq(as.Date("2021-01-01"), by = "week", length.out = 70)
  rm    <- fit_rm(make_chem("reference", dates), make_hydro(800))
  tgt   <- make_chem("target", dates, mult = 10, seed = 2)
  tm <- fit_target_model(tgt, rm, conc_units = "ug/L", min_obs_model = 10L,
                         pool = TRUE, api_windows_short = c(7L, 14L),
                         api_windows_long = c(30L, 90L))
  expect_s3_class(tm, "target_model")
  res <- leachatetools:::.resolve_target_impact(
    tm, tibble::tibble(date = seq(as.Date("2021-03-01"), as.Date("2021-09-01"),
                                  by = "month"))
  )
  expect_true(all(is.finite(res$impact)))
})

test_that("pool = TRUE fires the pooled 'model' tier on hydro-driven impact", {
  hydro <- make_hydro(900)
  dates <- seq(as.Date("2021-01-01"), by = "week", length.out = 80)
  rm    <- fit_rm(make_chem("reference", dates), hydro)
  tgt   <- make_hydro_driven("target", dates, hydro, mult = 10, seed = 5)
  tm <- fit_target_model(tgt, rm, conc_units = "ug/L", min_obs_model = 10L,
                         pool = TRUE, api_windows_short = c(7L, 14L),
                         api_windows_long = c(30L, 90L))
  tiers <- vapply(tm$models, `[[`, character(1), "tier")
  # at least one analyte's response is pooled-modelled
  expect_true(any(tiers == "model"))
  pooled <- Filter(function(m) isTRUE(m$pooled), tm$models)
  expect_true(length(pooled) >= 1L)
  # pooled analytes share one common window (the hallmark of a single joint fit)
  expect_equal(length(unique(vapply(pooled, `[[`, integer(1), "window_short"))), 1L)
})

test_that("pool = TRUE with a single modelled analyte falls back gracefully", {
  dates <- seq(as.Date("2021-01-01"), by = "week", length.out = 40)
  rm    <- fit_rm(make_chem("reference", dates), make_hydro())
  # one toxicant only
  tgt <- make_chem("target", dates, mult = 10, seed = 2) |>
    dplyr::filter(!analyte %in% c("Zn", "Ni"))
  expect_no_error(
    fit_target_model(tgt, rm, conc_units = "ug/L", min_obs_model = 8L,
                     pool = TRUE, api_windows_short = c(7L),
                     api_windows_long = c(30L))
  )
})

# Four large hydro-driven toxicants + one tiny one (Ni).  The large siblings
# form a high-amplitude "population" hydro response; the bug drags the tiny
# analyte toward it.  All share the same first-flush shape so they pool.
make_multi_scale <- function(site, dates, hydro, big_mult, ni_mult, seed = 7) {
  set.seed(seed)
  api7 <- leachatetools:::.compute_api(hydro$value, hydro$date, dates, 7L)
  api7 <- (api7 - mean(api7)) / (stats::sd(api7) + 1e-9)
  bases <- c(Cu = 0.5, Zn = 5, Pb = 0.8, Cr = 1.0)           # large toxicants
  purrr::pmap_dfr(list(dates, api7), function(d, a) {
    flush <- exp(0.8 * a)
    big <- vapply(bases, function(b) exp(rnorm(1, log(b), 0.1)) * big_mult * flush,
                  numeric(1))
    tibble::tibble(
      sample_id = paste0(site, format(d, "%Y%m%d")), site_id = site,
      datetime = d,
      analyte = c(names(bases), "Ni", "pH", "DOC", "hardness", "Ca", "Mg"),
      value = c(unname(big),
                exp(rnorm(1, log(0.3), 0.1)) * ni_mult * flush,  # tiny Ni
                7, 3, 40, 8, 5),                                 # constant co-analytes
      detected = TRUE)
  })
}

test_that("pool = TRUE preserves per-analyte magnitude (no cross-contamination)", {
  # Regression test for the pooling bug where a single bs='fs' fit on the raw
  # impact shrank each analyte's response toward the population: with several
  # large-signal toxicants present, a near-zero analyte (Ni) was dragged up
  # toward their high-amplitude hydro response (on the real B.S01 data Ni's
  # daily AmsPAF share jumped from ~0% to ~40% while Cu collapsed from ~45% to
  # 0%).  The fix pools the standardised (z) SHAPE and restores each analyte's
  # own magnitude.
  #
  # The bug is visible BETWEEN anchors where the pooled hydro response dominates
  # (the daily-grid-over-sparse-grabs case).  The smoother clips to the grab span
  # (no prediction beyond the last anchor), so we query daily WITHIN the 2021
  # span, between the weekly anchors: the day-to-day variation is the pooled
  # hydro response fitted_I plus each analyte's own smoothed residual — the
  # magnitude the bug inflates for Ni.
  hydro <- make_hydro(900)                                   # spans 2020-2022
  dates <- seq(as.Date("2021-01-01"), by = "week", length.out = 52)
  rm    <- fit_rm(make_multi_scale("reference", dates, hydro, 1, 1), hydro)
  tgt   <- make_multi_scale("target", dates, hydro, big_mult = 20, ni_mult = 1.1,
                            seed = 8)
  tm <- fit_target_model(tgt, rm, conc_units = "ug/L", min_obs_model = 10L,
                         pool = TRUE, api_windows_short = c(7L, 14L),
                         api_windows_long = c(30L, 90L))
  # Ni must actually be jointly pooled with the big siblings, else vacuous.
  expect_true(isTRUE(tm$models$Ni$pooled) && isTRUE(tm$models$Cu$pooled))

  q   <- tibble::tibble(date = seq(as.Date("2021-02-03"), as.Date("2021-11-01"),
                                   by = "day"))
  res <- leachatetools:::.resolve_target_impact(tm, q)
  big <- stats::sd(res$impact[res$analyte == "Cu"])
  ni  <- stats::sd(res$impact[res$analyte == "Ni"])
  # Ni's pooled hydro-response swing must stay far below the big toxicants'.
  # Under the old raw-I pool Ni was inflated toward the population amplitude.
  expect_gt(big, ni * 5)
})


## ── impact_tier surfaced in ara_summary() (issue #14, item A) ─────────────────

test_that("amspaf_daily(interpolation='model') surfaces impact_tier in ara_summary()", {
  dates <- seq(as.Date("2021-01-01"), by = "2 weeks", length.out = 40)
  rm    <- fit_rm(make_chem("reference", dates), make_hydro())
  d     <- amspaf_daily(make_chem("target", dates, mult = 10, seed = 2),
                        reference = rm, reference_model = rm,
                        interpolation = "model", conc_units = "ug/L",
                        require_temperature = FALSE, min_analytes = 1)
  s <- ara_summary(d)
  expect_true("impact_tier" %in% names(s))
  expect_true(all(stats::na.omit(s$impact_tier) %in% c("model", "bridge")))
})


## ── impute-first front-end (brms smoke test, Stan-gated) ──────────────────────

test_that("impute-first enriches reference & target anchors (brms smoke test)", {
  skip_if_not(
    requireNamespace("brms", quietly = TRUE) &&
      identical(Sys.getenv("BRMS_SMOKE_TEST"), "1"),
    "Skipping brms smoke test: brms not installed or BRMS_SMOKE_TEST != '1'"
  )

  # Build chemistry where some samples are missing a metal so imputation has
  # something to fill. pH/EC present everywhere (required_vars).
  set.seed(11)
  dates <- seq(as.Date("2021-01-01"), by = "week", length.out = 60)
  base  <- make_chem("reference", dates, mult = 1, seed = 11)
  # add pH/EC duplicates as required_vars (EC not in make_chem) and drop ~30% Cu
  ec <- base |>
    dplyr::distinct(sample_id, site_id, datetime) |>
    dplyr::mutate(analyte = "EC", value = runif(dplyr::n(), 100, 400), detected = TRUE)
  ref <- dplyr::bind_rows(base, ec)
  drop_cu <- sample(unique(ref$sample_id), 18)
  ref <- dplyr::filter(ref, !(sample_id %in% drop_cu & analyte == "Cu"))

  im <- fit_imputation_model(ref, required_vars = c("pH", "EC"),
                             iter = 400, warmup = 200, chains = 1, cores = 1)
  hydro <- make_hydro()

  # Reference model with vs without impute-first: impute-first must not have
  # fewer Cu anchors.
  rm_plain <- fit_reference_model(ref, hydro = hydro, conc_units = "ug/L",
                                  min_obs_model = 10L,
                                  api_windows_short = c(7L), api_windows_long = c(30L))
  rm_imp   <- fit_reference_model(ref, hydro = hydro, conc_units = "ug/L",
                                  imputation_model = im, min_obs_model = 10L,
                                  api_windows_short = c(7L), api_windows_long = c(30L))
  n_plain <- if (!is.null(rm_plain$models[["Cu"]])) rm_plain$models[["Cu"]]$n_obs else 0L
  n_imp   <- if (!is.null(rm_imp$models[["Cu"]]))   rm_imp$models[["Cu"]]$n_obs   else 0L
  expect_gte(n_imp, n_plain)

  # Target model accepts the imputation model and still fits.
  tm <- fit_target_model(dplyr::mutate(ref, site_id = "target"), rm_imp,
                         imputation_model = im, conc_units = "ug/L",
                         min_obs_model = 8L,
                         api_windows_short = c(7L), api_windows_long = c(30L))
  expect_s3_class(tm, "target_model")
  expect_true(length(tm$models) > 0L)
})


test_that("amspaf_daily(interpolation='model'): ARA <= no-ARA, daily series", {
  dates <- seq(as.Date("2021-01-01"), by = "2 weeks", length.out = 40)
  ref   <- make_chem("reference", dates)
  hydro <- make_hydro()
  rm    <- fit_rm(ref, hydro)
  tgt   <- make_chem("target", dates, mult = 10, seed = 2)

  d_ara <- amspaf_daily(tgt, reference = rm, reference_model = rm,
                        interpolation = "model", conc_units = "ug/L",
                        require_temperature = FALSE, min_analytes = 1)
  d_tot <- amspaf_daily(tgt, reference = NULL, reference_model = rm,
                        interpolation = "model", conc_units = "ug/L",
                        require_temperature = FALSE, min_analytes = 1)

  expect_s3_class(d_ara, "tbl_df")
  expect_true(nrow(d_ara) > 100L)           # genuinely daily
  expect_true(all(c("date", "site_id", "amspaf") %in% names(d_ara)))
  # Impact (ARA) cannot exceed total (no ARA), on average
  expect_lte(mean(d_ara$amspaf), mean(d_tot$amspaf) + 1e-6)
  # ara_summary attribute survives
  expect_false(is.null(attr(d_ara, "ara_summary")))
})
