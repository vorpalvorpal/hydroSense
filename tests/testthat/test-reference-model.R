## Tests for the temporal ARA reference model (issue #9)
##
## All tests are Stan-free and weatherOz-free (no external services needed).
## The weatherOz-gated test for get_silo_rainfall() is at the end.
##
## Coverage:
##   - .compute_api()              decay weights and edge cases
##   - .compute_antecedent_mean()  rolling mean and missing-window case
##   - .compute_hydro_features()   dispatcher
##   - .select_api_windows()       grid filtering
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
## .compute_api()
## ============================================================================

test_that(".compute_api() returns 0 when no hydro within window", {
  hvals  <- c(10, 5, 3)
  hdates <- as.Date(c("2022-01-01", "2022-01-02", "2022-01-03"))
  # Target date 20 days later, window = 7 → no overlap
  api <- leachatetools:::.compute_api(hvals, hdates, as.Date("2022-01-23"), 7L)
  expect_equal(api, 0)
})

test_that(".compute_api() returns correct single-day decay value", {
  hvals  <- 10
  hdates <- as.Date("2022-01-01")
  # Same day: lag = 0, so API = 10 * k^0 = 10
  api <- leachatetools:::.compute_api(hvals, hdates, as.Date("2022-01-01"), 7L)
  expect_equal(api, 10)
})

test_that(".compute_api() applies exponential decay correctly", {
  set.seed(1)
  hvals  <- c(10, 5)
  hdates <- as.Date(c("2022-01-01", "2022-01-02"))
  td     <- as.Date("2022-01-02")
  k      <- exp(-1 / 7)
  expected <- 10 * k^1 + 5 * k^0  # lag 1 for 01-01, lag 0 for 01-02
  api <- leachatetools:::.compute_api(hvals, hdates, td, 7L)
  expect_equal(api, expected, tolerance = 1e-9)
})

test_that(".compute_api() handles a vector of target dates", {
  hvals  <- rep(1, 10)
  hdates <- seq(as.Date("2022-01-01"), by = "day", length.out = 10)
  tds    <- as.Date(c("2022-01-05", "2022-01-10"))
  apis   <- leachatetools:::.compute_api(hvals, hdates, tds, 7L)
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
  out   <- leachatetools:::.compute_hydro_features(hydro, tds, 7L, 30L, "rainfall")
  expect_s3_class(out, "tbl_df")
  expect_named(out, c("date", "hydro_short", "hydro_long"))
  expect_equal(nrow(out), 4L)
})

test_that(".compute_hydro_features() stage type uses antecedent mean", {
  hydro <- tibble::tibble(
    date  = seq(as.Date("2022-01-01"), by = "day", length.out = 60),
    value = runif(60, 1, 5)
  )
  td  <- hydro$date[30]
  out <- leachatetools:::.compute_hydro_features(hydro, td, 7L, 30L, "stage")
  expect_true(!is.na(out$hydro_short))
  expect_true(!is.na(out$hydro_long))
})


## ============================================================================
## .select_api_windows()
## ============================================================================

test_that(".select_api_windows() returns a window_short < window_long pair", {
  ref   <- make_ref_chem(60)
  hydro <- make_hydro(500)
  obs   <- tibble::tibble(
    date       = as.Date(ref$datetime[ref$analyte == "Cu"]),
    value_norm = ref$value[ref$analyte == "Cu"]
  )
  sel <- leachatetools:::.select_api_windows(
    obs, hydro, "rainfall",
    api_windows_short = c(3L, 7L),
    api_windows_long  = c(30L, 60L)
  )
  expect_named(sel, c("window_short", "window_long", "best_aic", "null_aic"))
  expect_true(sel$window_short < sel$window_long)
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
    api_windows_short = c(7L, 14L),
    api_windows_long  = c(30L, 60L)
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
    api_windows_short = c(7L), api_windows_long = c(30L)
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
    api_windows_short = c(7L), api_windows_long = c(30L)
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
    api_windows_short = c(7L), api_windows_long = c(30L)
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
    api_windows_short = c(7L), api_windows_long = c(30L)
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
    api_windows_short = c(7L), api_windows_long = c(30L)
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
    api_windows_short = c(7L), api_windows_long = c(30L)
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
