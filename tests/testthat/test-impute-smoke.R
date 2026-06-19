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
  samples <- paste0("s", seq_len(n))
  metals <- c("Cu", "Zn", "Ni")
  drivers <- c("pH", "EC", "NH3-N", "DOC")

  rows <- tidyr::expand_grid(
    sample_id = samples,
    analyte   = c(metals, drivers)
  ) |>
    dplyr::mutate(
      site_id = "f1",
      datetime = as.Date("2023-01-01") + (match(sample_id, samples) - 1L),
      value = dplyr::case_when(
        analyte == "pH" ~ runif(dplyr::n(), 6.5, 8.5),
        analyte == "EC" ~ runif(dplyr::n(), 100, 500),
        analyte == "NH3-N" ~ runif(dplyr::n(), 0.01, 0.5),
        analyte == "DOC" ~ runif(dplyr::n(), 0.2, 5),
        TRUE ~ exp(rnorm(dplyr::n(), log(2), 0.5))
      ),
      detected = TRUE,
      imputed = FALSE
    )

  # Introduce some BDL rows for metals (not drivers)
  if (n_bdl > 0L) {
    bdl_idx <- sample(which(rows$analyte %in% metals), n_bdl)
    rows$detected[bdl_idx] <- FALSE
  }

  # Introduce missing metal rows (remove entirely). Guard n_missing == 0:
  # negative indexing by integer(0) (`rows[-integer(0), ]`) selects *no* rows.
  if (n_missing > 0L) {
    miss_idx <- sample(which(rows$analyte %in% metals), n_missing)
    rows <- rows[-miss_idx, ]
  }

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
    m <- fit_imputation_model(df,
      required_vars = c("pH", "EC"),
      iter = 100, warmup = 50, chains = 1
    ),
    regexp = "No target analytes"
  )
  expect_length(m$groups, 0L)
  expect_null(m$pca)
  expect_s3_class(m, "imputation_model")
})

test_that("fit_imputation_model full run (brms smoke test)", {
  skip_if_not(
    requireNamespace("brms", quietly = TRUE) &&
      identical(Sys.getenv("BRMS_SMOKE_TEST"), "1"),
    "Skipping brms smoke test: brms not installed or BRMS_SMOKE_TEST != '1'"
  )

  df <- make_impute_chem(n = 30, n_bdl = 4, n_missing = 3)
  m <- fit_imputation_model(
    df,
    required_vars = c("pH", "EC"),
    iter = 500,
    warmup = 250,
    chains = 1,
    cores = 1
  )

  expect_s3_class(m, "imputation_model")
  expect_true(!is.null(m$groups$metals))
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
    fit_imputation_model(df,
      required_vars = c("pH", "EC"),
      impute_method = "bogus"
    ),
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
    m <- fit_imputation_model(df,
      required_vars = c("pH", "EC"),
      impute_method = meth,
      iter = 400, warmup = 200, chains = 1, cores = 1
    )
    expect_equal(m$impute_method, meth)
    expect_equal(m$groups$metals$impute_method, meth)
    imp <- impute_chemistry(df, m)
    expect_true(all(is.finite(imp$value)))
    expect_true(any(imp$imputed))
    expect_true(all(imp$imputed_kind %in% c("observed", "censored_left", "missing")))

    # cens_factor must carry a genuine shared per-sample latent factor: the
    # long-format model has a `sample_id` group-level effect. (The old wide
    # `(1 |q| sample_id)` form silently produced no such coupling.)
    if (meth == "cens_factor") {
      vc <- brms::VarCorr(m$groups$metals$fit)
      expect_true("sample_id" %in% names(vc))
      expect_gt(vc$sample_id$sd[1, "Estimate"], 0)
    }
  }
})

test_that("apply_hurdles = FALSE imputes samples that carry no hurdle analytes (brms smoke test)", {
  skip_if_not(
    requireNamespace("brms", quietly = TRUE) &&
      identical(Sys.getenv("BRMS_SMOKE_TEST"), "1"),
    "Skipping brms smoke test: brms not installed or BRMS_SMOKE_TEST != '1'"
  )

  # Build a data set where one sample is metals-only (has Cu: passes hurdle)
  # and a second sample has NO metals at all (no hurdle analytes).
  df_full <- make_impute_chem(n = 30, n_bdl = 0, n_missing = 0)

  # Remove all metals from sample "s1" so it has no hurdle analytes.
  df_no_metal <- dplyr::filter(df_full, !(sample_id == "s1" & analyte %in% c("Cu", "Zn", "Ni")))

  m <- fit_imputation_model(
    df_full,
    required_vars = c("pH", "EC"),
    iter = 500, warmup = 250, chains = 1, cores = 1
  )

  # apply_hurdles = TRUE (default): s1 has no metals → not eligible → no metal rows added
  imp_with <- impute_chemistry(df_no_metal, m, apply_hurdles = TRUE)
  s1_metal_with <- dplyr::filter(
    imp_with, sample_id == "s1",
    analyte %in% c("Cu", "Zn", "Ni")
  )
  expect_equal(nrow(s1_metal_with), 0L)

  # apply_hurdles = FALSE: s1 gets imputed regardless of missing hurdle analytes
  imp_without <- impute_chemistry(df_no_metal, m, apply_hurdles = FALSE)
  s1_metal_without <- dplyr::filter(
    imp_without, sample_id == "s1",
    analyte %in% c("Cu", "Zn", "Ni")
  )
  expect_gt(nrow(s1_metal_without), 0L)
  expect_true(all(s1_metal_without$imputed))
})

test_that("fabricated rows carry the 'missing' kind, detected flag, and sample metadata (brms smoke test)", {
  # Plan #53 reqs 1-2: entirely-absent target cells for eligible samples gain
  # new rows tagged imputed_kind = "missing", detected = TRUE, with sample-level
  # metadata (site_id, datetime) carried from that sample's existing rows.
  skip_if_not(
    requireNamespace("brms", quietly = TRUE) &&
      identical(Sys.getenv("BRMS_SMOKE_TEST"), "1"),
    "Skipping brms smoke test: brms not installed or BRMS_SMOKE_TEST != '1'"
  )

  df_full <- make_impute_chem(n = 30, n_bdl = 0, n_missing = 0)
  # s1 has all three metals removed entirely -> they are "missing" for s1.
  df_no_metal <- dplyr::filter(
    df_full, !(sample_id == "s1" & analyte %in% c("Cu", "Zn", "Ni"))
  )
  m <- fit_imputation_model(
    df_full,
    required_vars = c("pH", "EC"),
    iter = 500, warmup = 250, chains = 1, cores = 1
  )

  imp <- impute_chemistry(df_no_metal, m, apply_hurdles = FALSE)
  s1_metals <- dplyr::filter(
    imp, sample_id == "s1", analyte %in% c("Cu", "Zn", "Ni")
  )

  # All three absent metals are fabricated, once each.
  expect_setequal(s1_metals$analyte, c("Cu", "Zn", "Ni"))
  # Correct flags on every fabricated row.
  expect_true(all(s1_metals$imputed))
  expect_true(all(s1_metals$imputed_kind == "missing"))
  expect_true(all(s1_metals$detected))
  expect_true(all(is.finite(s1_metals$value) & s1_metals$value > 0))

  # Sample metadata carried from s1's surviving (driver) rows.
  s1_carrier <- dplyr::slice(dplyr::filter(df_no_metal, sample_id == "s1"), 1L)
  expect_true(all(s1_metals$site_id == s1_carrier$site_id))
  expect_true(all(s1_metals$datetime == s1_carrier$datetime))

  # Req 5: observed driver rows for s1 are untouched by fabrication.
  ph_in <- dplyr::filter(df_no_metal, sample_id == "s1", analyte == "pH")$value
  ph_out <- dplyr::filter(
    imp, sample_id == "s1", analyte == "pH", !imputed
  )$value
  expect_equal(sort(ph_out), sort(ph_in))

  # Edge case: a complete frame (every sample has every metal) yields no
  # "missing" fabrication at all.
  imp_complete <- impute_chemistry(df_full, m, apply_hurdles = FALSE)
  expect_false(any(imp_complete$imputed_kind == "missing"))
})

test_that("fabrication in draws mode emits one row per draw for each absent cell (brms smoke test)", {
  # Plan #53 req 4: in draws mode, each fabricated cell expands to one row per
  # draw_id, sourced solely from the posterior draws (no row duplication).
  skip_if_not(
    requireNamespace("brms", quietly = TRUE) &&
      identical(Sys.getenv("BRMS_SMOKE_TEST"), "1"),
    "Skipping brms smoke test: brms not installed or BRMS_SMOKE_TEST != '1'"
  )

  df_full <- make_impute_chem(n = 30, n_bdl = 0, n_missing = 0)
  df_no_metal <- dplyr::filter(
    df_full, !(sample_id == "s1" & analyte %in% c("Cu", "Zn", "Ni"))
  )
  m <- fit_imputation_model(
    df_full,
    required_vars = c("pH", "EC"),
    iter = 500, warmup = 250, chains = 1, cores = 1
  )

  nd <- 20L
  imp <- impute_chemistry(
    df_no_metal, m,
    apply_hurdles = FALSE, return = "draws", ndraws = nd
  )

  cu <- dplyr::filter(imp, sample_id == "s1", analyte == "Cu")
  expect_equal(nrow(cu), nd)
  expect_setequal(cu$draw_id, seq_len(nd))
  expect_true(all(cu$imputed_kind == "missing"))
  expect_true(all(is.finite(cu$value) & cu$value > 0))
})

test_that("impute_coanalytes skips targets not in pca_vars", {
  # Trivial test: construct an imputation_model with no PCA and confirm the
  # function gives a clear error rather than crashing.
  fake_model <- structure(
    list(pca = NULL, pca_vars = c("pH", "EC")),
    class = "imputation_model"
  )
  df <- make_impute_chem()
  expect_error(
    impute_coanalytes(df, fake_model),
    regexp = "no fitted PCA"
  )
})
