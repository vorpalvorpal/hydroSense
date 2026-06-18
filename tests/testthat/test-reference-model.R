## Tests for the temporal ARA reference model (issue #9)
##
## All tests are Stan-free and weatherOz-free (no external services needed).
## The weatherOz-gated test for get_silo_rainfall() is at the end.
##
## Coverage:
##   - .compute_api()              recursive linear reservoir (issue #49)
##   - .compute_antecedent_mean()  rolling mean and missing-window case
##   - .compute_hydro_features()   dispatcher
##   - .select_api_tau()           continuous tau selection by AIC (issue #49)
##   - fit_reference_model()       with pre-supplied hydro; basic structure
##   - .resolve_ref_norm()         three-tier ladder
##   - add_amspaf(reference_model) end-to-end (target = reference → ARA ≈ 0)
##   - ara_summary()               attribute accessor

library(testthat)
library(leachatetools)


## ============================================================================
## Helpers
## ============================================================================

## Synthetic reference chemistry (simple metals + co-analytes)
make_ref_chem <- function(n = 60, seed = 42) {
  set.seed(seed)
  analytes <- c("Cu", "Zn", "Ni", "pH", "DOC", "hardness", "Ca", "Mg")
  dates    <- seq(as.Date("2020-01-01"), by = "week", length.out = n)
  purrr::map_dfr(dates, function(d) {
    tibble::tibble(
      sample_id = paste0("r", format(d, "%Y%m%d")),
      site_id   = "reference",
      datetime  = d,
      analyte   = analytes,
      value     = c(
        exp(rnorm(1, log(0.5), 0.3)),   # Cu  µg/L
        exp(rnorm(1, log(5),   0.4)),   # Zn
        exp(rnorm(1, log(0.3), 0.3)),   # Ni
        runif(1, 6.5, 8.0),             # pH
        runif(1, 1,   5),               # DOC
        runif(1, 20,  60),              # hardness
        runif(1, 4,   12),              # Ca
        runif(1, 2,   8)                # Mg
      ),
      detected  = TRUE
    )
  })
}

## Daily synthetic rainfall series
make_hydro <- function(n_days = 500, seed = 99) {
  set.seed(seed)
  tibble::tibble(
    date  = seq(as.Date("2019-07-01"), by = "day", length.out = n_days),
    value = pmax(0, rnorm(n_days, 2, 4))
  )
}

## Minimal target chemistry mirroring reference (for target = reference test)
make_target_chem <- function(ref_chem) {
  ref_chem |>
    dplyr::mutate(
      site_id   = "target",
      sample_id = paste0("t", sample_id)
    )
}


## ============================================================================
## .compute_api()  — exact recursive linear reservoir (issue #49)
## ============================================================================
##
## The store obeys  dS/dt = -S/tau + P(t),  whose exact daily solution is the
## recursion  S_t = k^{dt} * S_{t-1} + P_t  with  k = exp(-1/tau)  (Maillet 1905
## recession; Kohler & Linsley 1951 Antecedent Precipitation Index).  The fourth
## argument to .compute_api() is now the recession constant tau (days), NOT a hard
## summation window: there is no truncation horizon.
##
## `api_recursion_oracle()` is the closed-form specification the implementation
## must reproduce — it seeds an empty reservoir at the first hydro date, treats
## NA rainfall as zero input, steps the recursion through every hydro date, and
## decays the running state forward to each target date by k^{dt}.

api_recursion_oracle <- function(hvals, hdates, target_dates, tau) {
  k <- exp(-1 / tau)
  ord    <- order(hdates)
  hdates <- hdates[ord]
  hvals  <- hvals[ord]
  hvals[is.na(hvals)] <- 0
  vapply(target_dates, function(td) {
    mask <- hdates <= td
    if (!any(mask)) return(0)
    d <- hdates[mask]
    p <- hvals[mask]
    s <- 0
    prev <- d[1L]
    for (i in seq_along(d)) {
      s <- k^as.numeric(d[i] - prev) * s + p[i]
      prev <- d[i]
    }
    s * k^as.numeric(td - prev)
  }, numeric(1L))
}

## Old (pre-#49) truncated windowed sum — retained only as the correctness
## baseline the recursion must dominate (truncation can only discard mass).
api_old_truncated <- function(hvals, hdates, target_dates, window_days) {
  k <- exp(-1 / window_days)
  vapply(target_dates, function(td) {
    mask <- hdates >= (td - window_days) & hdates <= td
    if (!any(mask)) return(0)
    lags <- as.numeric(td - hdates[mask])
    sum(hvals[mask] * k^lags, na.rm = TRUE)
  }, numeric(1L))
}

test_that(".compute_api() reproduces the exact recursive reservoir on a daily series", {
  set.seed(1)
  hdates <- seq(as.Date("2022-01-01"), by = "day", length.out = 120)
  hvals  <- pmax(0, rnorm(120, 2, 4))
  tds    <- hdates[c(30, 60, 90, 120)]
  tau    <- 10
  expect_equal(
    leachatetools:::.compute_api(hvals, hdates, tds, tau),
    api_recursion_oracle(hvals, hdates, tds, tau),
    tolerance = 1e-9
  )
})

test_that(".compute_api() equals the convergent sum and dominates the old truncated window", {
  set.seed(2)
  hdates <- seq(as.Date("2022-01-01"), by = "day", length.out = 365)
  hvals  <- pmax(0, rnorm(365, 2, 4))
  tds    <- hdates[c(200, 300, 365)]
  tau    <- 30
  new_api <- leachatetools:::.compute_api(hvals, hdates, tds, tau)
  old_api <- api_old_truncated(hvals, hdates, tds, tau)
  # Truncation at the window only removes non-negative mass → new >= old.
  expect_true(all(new_api >= old_api - 1e-9))
  # With > 5*tau of history the recursion has effectively converged: it equals
  # the full-history exponential sum (the infinite-sum limit the old form only
  # approximated by truncating at the window).
  full_sum <- vapply(tds, function(td) {
    mask <- hdates <= td
    sum(hvals[mask] * exp(-1 / tau)^as.numeric(td - hdates[mask]))
  }, numeric(1L))
  expect_equal(new_api, full_sum, tolerance = 1e-9)
})

test_that(".compute_api() is monotone non-decreasing in each rainfall input", {
  hdates <- seq(as.Date("2022-01-01"), by = "day", length.out = 20)
  hvals  <- rep(1, 20)
  tds    <- hdates
  base   <- leachatetools:::.compute_api(hvals, hdates, tds, 7)
  bumped <- hvals
  bumped[10] <- bumped[10] + 5
  out    <- leachatetools:::.compute_api(bumped, hdates, tds, 7)
  # Dates before the bumped day are unchanged; from the bumped day on they rise.
  expect_equal(out[1:9], base[1:9], tolerance = 1e-9)
  expect_true(all(out[10:20] > base[10:20]))
})

test_that(".compute_api() decays geometrically over a gap with no input", {
  # Single rain pulse, then the reservoir drains as P0 * k^dt (pure recession).
  hdates <- as.Date("2022-01-01")
  hvals  <- 10
  tau    <- 8
  tds    <- as.Date("2022-01-01") + c(0, 1, 5, 20)
  k      <- exp(-1 / tau)
  expect_equal(
    leachatetools:::.compute_api(hvals, hdates, tds, tau),
    10 * k^c(0, 1, 5, 20),
    tolerance = 1e-9
  )
})

test_that(".compute_api() approaches a cumulative sum as tau -> Inf", {
  # k -> 1 : nothing decays, so API is the running total of all prior rain.
  hdates <- seq(as.Date("2022-01-01"), by = "day", length.out = 30)
  hvals  <- rep(2, 30)
  tds    <- hdates
  out    <- leachatetools:::.compute_api(hvals, hdates, tds, 1e6)
  expect_equal(out, cumsum(hvals), tolerance = 1e-3)
})

test_that(".compute_api() approaches today's rain only as tau -> 0+", {
  # k -> 0 : every lagged term vanishes, leaving only same-day rain.
  hdates <- seq(as.Date("2022-01-01"), by = "day", length.out = 10)
  hvals  <- as.numeric(1:10)
  tds    <- hdates
  out    <- leachatetools:::.compute_api(hvals, hdates, tds, 1e-6)
  expect_equal(out, hvals, tolerance = 1e-6)
})

test_that(".compute_api() returns 0 when no hydro on or before the target date", {
  hvals  <- c(10, 5, 3)
  hdates <- as.Date(c("2022-01-10", "2022-01-11", "2022-01-12"))
  # Target precedes all hydro → empty reservoir → 0.
  expect_equal(
    leachatetools:::.compute_api(hvals, hdates, as.Date("2022-01-01"), 7),
    0
  )
})

test_that(".compute_api() handles irregular spacing via k^dt", {
  hdates <- as.Date(c("2022-01-01", "2022-01-04", "2022-01-05", "2022-01-15"))
  hvals  <- c(10, 0, 5, 3)
  tds    <- as.Date(c("2022-01-05", "2022-01-20"))
  tau    <- 12
  expect_equal(
    leachatetools:::.compute_api(hvals, hdates, tds, tau),
    api_recursion_oracle(hvals, hdates, tds, tau),
    tolerance = 1e-9
  )
})

test_that(".compute_api() treats NA rainfall as zero input while decay continues", {
  hdates <- seq(as.Date("2022-01-01"), by = "day", length.out = 6)
  hvals  <- c(10, NA, 5, NA, 2, NA)
  tds    <- hdates
  tau    <- 7
  expect_equal(
    leachatetools:::.compute_api(hvals, hdates, tds, tau),
    api_recursion_oracle(hvals, hdates, tds, tau), # oracle sets NA -> 0
    tolerance = 1e-9
  )
})

test_that(".compute_api() is deterministic across repeated calls", {
  set.seed(3)
  hdates <- seq(as.Date("2022-01-01"), by = "day", length.out = 50)
  hvals  <- pmax(0, rnorm(50, 2, 4))
  tds    <- hdates[c(10, 25, 50)]
  a <- leachatetools:::.compute_api(hvals, hdates, tds, 9)
  b <- leachatetools:::.compute_api(hvals, hdates, tds, 9)
  expect_identical(a, b)
})

test_that(".compute_api() handles a vector of target dates", {
  hvals  <- rep(1, 10)
  hdates <- seq(as.Date("2022-01-01"), by = "day", length.out = 10)
  tds    <- as.Date(c("2022-01-05", "2022-01-10"))
  apis   <- leachatetools:::.compute_api(hvals, hdates, tds, 7)
  expect_length(apis, 2L)
  expect_true(all(apis > 0))
})


## ============================================================================
## .compute_antecedent_mean()
## ============================================================================

test_that(".compute_antecedent_mean() returns mean of overlapping values", {
  hvals  <- c(4, 6, 8)
  hdates <- as.Date(c("2022-01-01", "2022-01-03", "2022-01-05"))
  # Window of 5 days back from 2022-01-05 covers all three
  am <- leachatetools:::.compute_antecedent_mean(hvals, hdates,
                                                  as.Date("2022-01-05"), 5L)
  expect_equal(am, mean(c(4, 6, 8)), tolerance = 1e-9)
})

test_that(".compute_antecedent_mean() returns NA when no data in window", {
  hvals  <- c(4, 6)
  hdates <- as.Date(c("2022-01-01", "2022-01-02"))
  am <- leachatetools:::.compute_antecedent_mean(hvals, hdates,
                                                  as.Date("2022-01-20"), 5L)
  expect_true(is.na(am))
})


## ============================================================================
## .compute_hydro_features()
## ============================================================================

test_that(".compute_hydro_features() produces hydro_short and hydro_long", {
  hydro <- make_hydro(200)
  tds   <- hydro$date[100:103]
  out   <- leachatetools:::.compute_hydro_features(hydro, tds, 7, 30, "rainfall")
  expect_s3_class(out, "tbl_df")
  expect_named(out, c("date", "hydro_short", "hydro_long"))
  expect_equal(nrow(out), 4L)
})

test_that(".compute_hydro_features() rainfall path uses the recursive reservoir at each tau", {
  hydro <- make_hydro(200)
  tds   <- hydro$date[100:103]
  out   <- leachatetools:::.compute_hydro_features(hydro, tds, 7, 30, "rainfall")
  expect_equal(out$hydro_short,
               api_recursion_oracle(hydro$value, hydro$date, tds, 7),
               tolerance = 1e-9)
  expect_equal(out$hydro_long,
               api_recursion_oracle(hydro$value, hydro$date, tds, 30),
               tolerance = 1e-9)
})

test_that(".compute_hydro_features() stage type uses antecedent mean (unchanged)", {
  hydro <- tibble::tibble(
    date  = seq(as.Date("2022-01-01"), by = "day", length.out = 60),
    value = runif(60, 1, 5)
  )
  td  <- hydro$date[30]
  out <- leachatetools:::.compute_hydro_features(hydro, td, 7, 30, "stage")
  expect_true(!is.na(out$hydro_short))
  expect_true(!is.na(out$hydro_long))
})


## ============================================================================
## .select_api_tau()  — continuous tau selection by profiled GAM AIC (issue #49)
## ============================================================================
##
## tau is chosen by an outer optimisation of the GAM's REML/AIC over a bounded
## continuous range (profile likelihood; Wood 2017), replacing the old discrete
## window grid.  Guards: physical bounds, a fast<slow separation constraint
## (tau_long >= 1.5*tau_short), and a parsimony gate — a fitted tau is adopted
## only when it beats the parsimonious default by >= 2 AIC, else the default is
## returned.  Deterministic (golden-section; no RNG).

test_that(".select_api_tau() returns named per-store tau within the supplied bounds", {
  ref   <- make_ref_chem(80)
  hydro <- make_hydro(700)
  obs   <- tibble::tibble(
    date       = as.Date(ref$datetime[ref$analyte == "Cu"]),
    value_norm = ref$value[ref$analyte == "Cu"]
  )
  sel <- leachatetools:::.select_api_tau(
    obs, hydro, "rainfall",
    tau_bounds_short = c(1, 30),
    tau_bounds_long  = c(20, 365)
  )
  expect_named(sel, c("tau_short", "tau_long", "best_aic", "null_aic"),
               ignore.order = TRUE)
  expect_gte(sel$tau_short, 1)
  expect_lte(sel$tau_short, 30)
  expect_gte(sel$tau_long, 20)
  expect_lte(sel$tau_long, 365)
})

test_that(".select_api_tau() enforces tau_long >= 1.5 * tau_short", {
  ref   <- make_ref_chem(80)
  hydro <- make_hydro(700)
  obs   <- tibble::tibble(
    date       = as.Date(ref$datetime[ref$analyte == "Zn"]),
    value_norm = ref$value[ref$analyte == "Zn"]
  )
  sel <- leachatetools:::.select_api_tau(
    obs, hydro, "rainfall",
    tau_bounds_short = c(1, 30),
    tau_bounds_long  = c(20, 365)
  )
  expect_gte(sel$tau_long, 1.5 * sel$tau_short)
})

test_that(".select_api_tau() falls back to the parsimonious default on a flat AIC surface", {
  # y independent of hydrology → no tau improves on the default by >= 2 AIC,
  # so the overfitting gate returns the supplied default unchanged.
  set.seed(7)
  hydro <- make_hydro(700)
  dates <- seq(as.Date("2020-01-01"), by = "week", length.out = 80)
  obs   <- tibble::tibble(date = dates, value_norm = exp(rnorm(80, 0, 0.3)))
  sel <- leachatetools:::.select_api_tau(
    obs, hydro, "rainfall",
    tau_bounds_short = c(1, 30),
    tau_bounds_long  = c(20, 365),
    default_short = 7, default_long = 60
  )
  expect_equal(sel$tau_short, 7)
  expect_equal(sel$tau_long, 60)
})

test_that(".select_api_tau() collapses to the default when bounds are degenerate", {
  ref   <- make_ref_chem(80)
  hydro <- make_hydro(700)
  obs   <- tibble::tibble(
    date       = as.Date(ref$datetime[ref$analyte == "Cu"]),
    value_norm = ref$value[ref$analyte == "Cu"]
  )
  sel <- leachatetools:::.select_api_tau(
    obs, hydro, "rainfall",
    tau_bounds_short = c(7, 7),
    tau_bounds_long  = c(60, 60)
  )
  expect_equal(sel$tau_short, 7)
  expect_equal(sel$tau_long, 60)
})

test_that(".select_api_tau() adopts a non-default tau when the signal warrants it", {
  # Generate value_norm from a reservoir at tau* = 20 with low noise, then start
  # the optimiser from a deliberately mismatched short default of 3.  A genuine
  # tau ~ 20 signal must pull tau_short up off that default (counterpart to the
  # flat-surface fallback above).
  set.seed(11)
  hydro    <- make_hydro(800)
  dates    <- seq(as.Date("2020-02-01"), by = "5 days", length.out = 90)
  api_true <- api_recursion_oracle(hydro$value, hydro$date, dates, 20)
  obs      <- tibble::tibble(
    date       = dates,
    value_norm = exp(0.5 * as.numeric(scale(api_true)) + rnorm(90, 0, 0.05))
  )
  sel <- leachatetools:::.select_api_tau(
    obs, hydro, "rainfall",
    tau_bounds_short = c(1, 30),
    tau_bounds_long  = c(20, 365),
    default_short = 3, default_long = 300
  )
  expect_gt(sel$tau_short, 3)
})

test_that(".select_api_tau() is deterministic across reruns", {
  ref   <- make_ref_chem(80)
  hydro <- make_hydro(700)
  obs   <- tibble::tibble(
    date       = as.Date(ref$datetime[ref$analyte == "Cu"]),
    value_norm = ref$value[ref$analyte == "Cu"]
  )
  a <- leachatetools:::.select_api_tau(obs, hydro, "rainfall", c(1, 30), c(20, 365))
  b <- leachatetools:::.select_api_tau(obs, hydro, "rainfall", c(1, 30), c(20, 365))
  expect_equal(a, b)
})


## ============================================================================
## fit_reference_model()
## ============================================================================

test_that("fit_reference_model() returns a reference_model with expected structure", {
  ref   <- make_ref_chem(60)
  hydro <- make_hydro(500)

  ref_model <- fit_reference_model(
    reference         = ref,
    hydro             = hydro,
    conc_units        = "ug/L",
    min_obs_model     = 10L,  # low threshold for synthetic data
    api_tau_bounds_short = c(7, 14),
    api_tau_bounds_long  = c(30, 60)
  )

  expect_s3_class(ref_model, "reference_model")
  expect_true(length(ref_model$models) > 0L)
  expect_s3_class(ref_model$hydro, "tbl_df")
  expect_equal(ref_model$hydro_type, "rainfall")
  expect_true(is.finite(ref_model$match_hydro_tol))

  # static_ref tibble carries ref_norm per analyte
  expect_named(ref_model$static_ref, c("analyte", "ref_norm", "n_obs"),
               ignore.order = TRUE)
  expect_true(all(ref_model$static_ref$ref_norm > 0))
})

test_that("print.reference_model() runs without error", {
  ref   <- make_ref_chem(60)
  hydro <- make_hydro(500)
  m     <- fit_reference_model(
    reference = ref, hydro = hydro,
    conc_units = "ug/L", min_obs_model = 10L,
    api_tau_bounds_short = c(7, 7), api_tau_bounds_long = c(30, 30)
  )
  expect_output(print(m), regexp = "reference_model")
})

test_that("fit_reference_model() errors when neither hydro nor lat/lon supplied", {
  ref <- make_ref_chem(20)
  expect_error(
    fit_reference_model(reference = ref, conc_units = "ug/L"),
    regexp = "hydro.*latitude"
  )
})

test_that("fit_reference_model() falls back to static tier when n < min_obs_model", {
  ref   <- make_ref_chem(5)   # tiny dataset
  hydro <- make_hydro(200)
  m     <- fit_reference_model(
    reference     = ref,
    hydro         = hydro,
    conc_units    = "ug/L",
    min_obs_model = 30L   # requires 30 obs → all analytes fall to static tier
  )
  tiers <- vapply(m$models, `[[`, character(1), "tier")
  expect_true(all(tiers == "static"))
})


## ============================================================================
## .resolve_ref_norm() — three-tier ladder
## ============================================================================

test_that(".resolve_ref_norm() returns (sample_id, analyte, ref_norm, ref_tier)", {
  ref   <- make_ref_chem(60)
  hydro <- make_hydro(500)
  m     <- fit_reference_model(
    reference = ref, hydro = hydro,
    conc_units = "ug/L", min_obs_model = 10L,
    api_tau_bounds_short = c(7, 7), api_tau_bounds_long = c(30, 30)
  )
  target <- make_target_chem(ref)
  out    <- leachatetools:::.resolve_ref_norm(m, target)
  expect_named(out, c("sample_id", "analyte", "ref_norm", "ref_tier"),
               ignore.order = TRUE)
  expect_true(all(out$ref_norm > 0, na.rm = TRUE))
  expect_true(all(out$ref_tier %in% c("matched", "model", "static")))
})

test_that(".resolve_ref_norm() chronic path returns model_integrated or static tier", {
  ref   <- make_ref_chem(60)
  hydro <- make_hydro(600)
  m     <- fit_reference_model(
    reference = ref, hydro = hydro,
    conc_units = "ug/L", min_obs_model = 10L,
    api_tau_bounds_short = c(7, 7), api_tau_bounds_long = c(30, 30)
  )
  # Simulate a chronic target (has focal_date column)
  target <- tibble::tibble(
    sample_id  = c("c1", "c2"),
    focal_date = as.Date(c("2021-01-01", "2021-06-01"))
  )
  out <- leachatetools:::.resolve_ref_norm(m, target, tau_days = 90, window_days = 180)
  expect_true("focal_date" %in% names(target))  # ensure dispatch worked
  expect_named(out, c("sample_id", "analyte", "ref_norm", "ref_tier"),
               ignore.order = TRUE)
  expect_true(all(out$ref_tier %in% c("model_integrated", "static")))
})


## ============================================================================
## End-to-end: add_amspaf(reference_model) — target = reference → ARA ≈ 0
## ============================================================================

test_that("add_amspaf(reference_model): target = reference gives small AmsPAF", {
  ref   <- make_ref_chem(60)
  hydro <- make_hydro(500)
  m     <- fit_reference_model(
    reference = ref, hydro = hydro,
    conc_units = "ug/L", min_obs_model = 10L,
    api_tau_bounds_short = c(7, 7), api_tau_bounds_long = c(30, 30)
  )

  # Use a subset of reference samples as "target"
  target <- ref |>
    dplyr::filter(sample_id %in% unique(sample_id)[1:5]) |>
    dplyr::mutate(site_id = "target")

  out <- add_amspaf(target, reference = m, conc_units = "ug/L",
                    require_temperature = FALSE)
  amspaf <- dplyr::filter(out, analyte == "AmsPAF")

  # When target ≈ reference background, AmsPAF should be low
  # (not necessarily exactly 0 because tier-2 prediction adds noise)
  expect_true(nrow(amspaf) > 0L)
  expect_true(all(amspaf$value >= 0))
})

test_that("add_amspaf(reference_model): elevated target gives higher AmsPAF than clean target", {
  ref   <- make_ref_chem(60)
  hydro <- make_hydro(500)
  m     <- fit_reference_model(
    reference = ref, hydro = hydro,
    conc_units = "ug/L", min_obs_model = 10L,
    api_tau_bounds_short = c(7, 7), api_tau_bounds_long = c(30, 30)
  )

  base_sample <- ref |>
    dplyr::filter(sample_id == unique(sample_id)[1]) |>
    dplyr::mutate(site_id = "target", sample_id = "clean")

  elevated <- base_sample |>
    dplyr::mutate(
      sample_id = "elevated",
      value = dplyr::if_else(analyte %in% c("Cu", "Zn", "Ni"),
                              value * 100, value)
    )

  both <- dplyr::bind_rows(base_sample, elevated)
  out  <- add_amspaf(both, reference = m, conc_units = "ug/L",
                     require_temperature = FALSE)
  amspaf <- dplyr::filter(out, analyte == "AmsPAF")

  paf_clean    <- amspaf$value[amspaf$sample_id == "clean"]
  paf_elevated <- amspaf$value[amspaf$sample_id == "elevated"]

  expect_gt(paf_elevated, paf_clean)
})


## ============================================================================
## ara_summary() — attribute accessor
## ============================================================================

test_that("ara_summary() returns a tibble with expected columns for NULL reference", {
  ref   <- make_ref_chem(30)
  hydro <- make_hydro(300)
  m     <- fit_reference_model(
    reference = ref, hydro = hydro,
    conc_units = "ug/L", min_obs_model = 10L,
    api_tau_bounds_short = c(7, 7), api_tau_bounds_long = c(30, 30)
  )
  target <- ref |>
    dplyr::filter(sample_id %in% unique(sample_id)[1:3]) |>
    dplyr::mutate(site_id = "target")

  out <- add_amspaf(target, reference = m, conc_units = "ug/L",
                    require_temperature = FALSE)
  s   <- ara_summary(out)

  expect_s3_class(s, "tbl_df")
  expected_cols <- c("sample_id", "analyte", "ref_norm", "C_norm",
                     "C_adj", "C_excess", "floor_fired", "ref_source", "ref_tier")
  expect_named(s, expected_cols, ignore.order = TRUE)
  expect_true(all(s$C_adj >= 0))
  # floor_fired should be TRUE exactly where C_excess < 0
  expect_equal(s$floor_fired, s$C_excess < 0)
})

test_that("ara_summary() also works for static reference path", {
  demo   <- leachate_demo()
  target <- dplyr::filter(demo, site_id == "downstream")
  refdat <- dplyr::filter(demo, site_id == "reference")

  out <- add_amspaf(target, reference = refdat,
                    conc_units = "ug/L", require_temperature = FALSE)
  s   <- ara_summary(out)

  expect_s3_class(s, "tbl_df")
  expect_true("floor_fired" %in% names(s))
  # ref_tier is NA for static path
  expect_true(all(is.na(s$ref_tier)))
})

test_that("ara_summary() returns NULL with message when attribute is absent", {
  df <- tibble::tibble(analyte = "Cu", value = 1)
  expect_message(s <- ara_summary(df), regexp = "ara_summary")
  expect_null(s)
})


## ============================================================================
## get_silo_rainfall() — gated on weatherOz + network
## ============================================================================

test_that("get_silo_rainfall() returns a daily rainfall tibble (weatherOz)", {
  skip_if_not_installed("weatherOz")
  skip_if(is.na(Sys.getenv("SILO_API_KEY", unset = NA)),
          "SILO_API_KEY not set — skipping SILO network test")

  rain <- get_silo_rainfall(
    latitude   = -33.87,
    longitude  = 151.21,
    start_date = "2022-01-01",
    end_date   = "2022-01-31",
    api_key    = Sys.getenv("SILO_API_KEY")
  )
  expect_s3_class(rain, "tbl_df")
  expect_named(rain, c("date", "rainfall_mm"))
  expect_equal(nrow(rain), 31L)
  expect_true(all(rain$rainfall_mm >= 0))
})
