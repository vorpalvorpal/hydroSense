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


## ── Residual bridge interpolation ─────────────────────────────────────────────

test_that(".interp_residual() pinches to the anchor value at an anchor date", {
  anchors <- tibble::tibble(
    date        = as.Date(c("2021-01-01", "2021-02-01", "2021-03-01")),
    S           = c(2, -1, 0.5),
    hydro_short = c(5, 5, 5),
    hydro_long  = c(10, 10, 10)
  )
  # query exactly on the middle anchor -> that anchor's residual
  v <- leachatetools:::.interp_residual(anchors, as.Date("2021-02-01"), 5, 10)
  expect_equal(v, -1, tolerance = 1e-6)
})

test_that(".interp_residual() leans toward the hydrologically similar bracket", {
  anchors <- tibble::tibble(
    date        = as.Date(c("2021-01-01", "2021-02-01")),
    S           = c(2, -1),
    hydro_short = c(1, 9),   # prev = dry, next = wet
    hydro_long  = c(1, 9)
  )
  qd <- as.Date("2021-01-16")  # midway in time
  # Dry query -> closer to the dry (prev) anchor's residual (+2)
  v_dry <- leachatetools:::.interp_residual(anchors, qd, 1, 1)
  # Wet query -> closer to the wet (next) anchor's residual (-1)
  v_wet <- leachatetools:::.interp_residual(anchors, qd, 9, 9)
  expect_gt(v_dry, v_wet)
  expect_gt(v_dry, 0)    # leans positive (toward dry anchor)
  expect_lt(v_wet, 0.5)  # leans toward the wet anchor
})

test_that(".interp_residual() handles single-anchor and empty cases", {
  one <- tibble::tibble(date = as.Date("2021-01-01"), S = 3,
                        hydro_short = 5, hydro_long = 10)
  expect_equal(leachatetools:::.interp_residual(one, as.Date("2021-06-01"), 5, 10), 3)
  none <- one[0, ]
  expect_equal(leachatetools:::.interp_residual(none, as.Date("2021-06-01"), 5, 10), 0)
})


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
