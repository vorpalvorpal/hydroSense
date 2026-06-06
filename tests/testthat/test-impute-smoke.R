## Smoke tests for fit_imputation_model() + impute_chemistry()
##
## Full brms/Stan tests are skipped unless brms + Stan are installed and the
## package-wide BRMS_SMOKE_TEST environment variable is set to "1".  The
## non-brms tests cover input validation and behaviour that does not require
## Stan compilation.

library(testthat)
library(leachatetools)

make_impute_chem <- function(n = 20, n_bdl = 3, n_missing = 2) {
  set.seed(7)
  samples  <- paste0("s", seq_len(n))
  metals   <- c("Cu", "Zn", "Ni")
  drivers  <- c("pH", "EC", "NH3-N", "DOC")

  rows <- tidyr::expand_grid(
    sample_id = samples,
    analyte   = c(metals, drivers)
  ) |>
    dplyr::mutate(
      site_id  = "f1",
      datetime = as.Date("2023-01-01") + (match(sample_id, samples) - 1L),
      value    = dplyr::case_when(
        analyte == "pH"     ~ runif(dplyr::n(), 6.5, 8.5),
        analyte == "EC"     ~ runif(dplyr::n(), 100, 500),
        analyte == "NH3-N"  ~ runif(dplyr::n(), 0.01, 0.5),
        analyte == "DOC"    ~ runif(dplyr::n(), 0.2, 5),
        TRUE                ~ exp(rnorm(dplyr::n(), log(2), 0.5))
      ),
      detected = TRUE,
      imputed  = FALSE
    )

  # Introduce some BDL rows for metals (not drivers)
  bdl_idx <- sample(which(rows$analyte %in% metals), n_bdl)
  rows$detected[bdl_idx] <- FALSE

  # Introduce missing metal rows (remove entirely)
  miss_idx <- sample(which(rows$analyte %in% metals), n_missing)
  rows <- rows[-miss_idx, ]

  rows
}

test_that(".require_brms() passes when brms is installed (brms is Suggests)", {
  # brms moved from Imports to Suggests: only the imputation step needs it.
  # The guard returns TRUE invisibly when present, else aborts with install
  # guidance (exercised in real use where brms is absent).
  skip_if_not_installed("brms")
  expect_true(leachatetools:::.require_brms())
})

test_that("fit_imputation_model errors when required_vars are missing entirely", {
  df <- make_impute_chem()
  df_no_ph <- dplyr::filter(df, analyte != "pH")
  # All samples drop → fewer than min_samples
  expect_error(
    fit_imputation_model(
      df_no_ph,
      required_vars = c("pH", "EC"),
      min_samples   = 10L
    ),
    regexp = "remain after"
  )
})

test_that("fit_imputation_model handles a fully-empty target set gracefully", {
  # Only required_vars present, no metals or organics → warn + empty model
  df <- make_impute_chem() |>
    dplyr::filter(analyte %in% c("pH", "EC", "NH3-N", "DOC"))
  expect_warning(
    m <- fit_imputation_model(df, required_vars = c("pH", "EC"),
                              iter = 100, warmup = 50, chains = 1),
    regexp = "No target analytes"
  )
  expect_null(m$metals)
  expect_null(m$organics)
  expect_s3_class(m, "imputation_model")
})

test_that("fit_imputation_model full run (brms smoke test)", {
  skip_if_not(
    requireNamespace("brms", quietly = TRUE) &&
      identical(Sys.getenv("BRMS_SMOKE_TEST"), "1"),
    "Skipping brms smoke test: brms not installed or BRMS_SMOKE_TEST != '1'"
  )

  df  <- make_impute_chem(n = 30, n_bdl = 4, n_missing = 3)
  m   <- fit_imputation_model(
    df,
    required_vars = c("pH", "EC"),
    iter    = 500,
    warmup  = 250,
    chains  = 1,
    cores   = 1
  )

  expect_s3_class(m, "imputation_model")
  expect_true(!is.null(m$metals))
  expect_true(!is.null(m$pca))

  imp <- impute_chemistry(df, m)

  # All values finite
  expect_true(all(is.finite(imp$value)))

  # imputed/imputed_kind columns present
  expect_true("imputed" %in% names(imp))
  expect_true("imputed_kind" %in% names(imp))
  expect_true(all(imp$imputed_kind %in% c("observed", "censored_left", "missing")))

  # Some BDL/missing rows should now be marked imputed
  expect_true(any(imp$imputed))
})

test_that("fit_imputation_model rejects an unknown impute_method", {
  skip_if_not_installed("brms")
  df <- make_impute_chem()
  # match.arg() runs before any Stan fitting, so this is fast and Stan-free.
  expect_error(
    fit_imputation_model(df, required_vars = c("pH","EC"),
                         impute_method = "bogus"),
    regexp = "should be one of"
  )
})

test_that("cens / cens_factor impute methods fit and impute (brms smoke test)", {
  skip_if_not(
    requireNamespace("brms", quietly = TRUE) &&
      identical(Sys.getenv("BRMS_SMOKE_TEST"), "1"),
    "Skipping brms smoke test: brms not installed or BRMS_SMOKE_TEST != '1'"
  )
  df <- make_impute_chem(n = 30, n_bdl = 4, n_missing = 3)
  for (meth in c("cens", "cens_factor")) {
    m <- fit_imputation_model(df, required_vars = c("pH","EC"),
                              impute_method = meth,
                              iter = 400, warmup = 200, chains = 1, cores = 1)
    expect_equal(m$impute_method, meth)
    expect_equal(m$metals$impute_method, meth)
    imp <- impute_chemistry(df, m)
    expect_true(all(is.finite(imp$value)))
    expect_true(any(imp$imputed))
    expect_true(all(imp$imputed_kind %in% c("observed","censored_left","missing")))

    # cens_factor must carry a genuine shared per-sample latent factor: the
    # long-format model has a `sample_id` group-level effect. (The old wide
    # `(1 |q| sample_id)` form silently produced no such coupling.)
    if (meth == "cens_factor") {
      vc <- brms::VarCorr(m$metals$fit)
      expect_true("sample_id" %in% names(vc))
      expect_gt(vc$sample_id$sd[1, "Estimate"], 0)
    }
  }
})

test_that("impute_coanalytes skips targets not in pca_vars", {
  # Trivial test: construct an imputation_model with no PCA and confirm the
  # function gives a clear error rather than crashing.
  fake_model <- structure(
    list(pca = NULL, pca_vars = c("pH","EC")),
    class = "imputation_model"
  )
  df <- make_impute_chem()
  expect_error(
    impute_coanalytes(df, fake_model),
    regexp = "no fitted PCA"
  )
})
