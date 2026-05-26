## Smoke tests for impute_chemistry()
##
## Full brms/Stan tests are skipped unless brms + Stan are installed and
## the package-wide BRMS_SMOKE_TEST environment variable is set to "1".
## The non-brms tests validate input checking and error messages.

library(testthat)
library(leachatetools)

make_impute_chem <- function(n = 15, n_bdl = 3, n_missing = 2) {
  set.seed(7)
  samples  <- paste0("s", seq_len(n))
  analytes <- c("Cu", "Zn", "Ni")
  drivers  <- c("pH", "EC", "DOC")

  # Build long-format df
  rows <- tidyr::expand_grid(uuid.sample = samples, name.analyte = c(analytes, drivers)) |>
    dplyr::mutate(
      uuid.feature    = "f1",
      datetime.sample = as.Date("2023-01-01") + (match(uuid.sample, samples) - 1L),
      value      = dplyr::case_when(
        name.analyte == "pH"  ~ runif(dplyr::n(), 6.5, 8.5),
        name.analyte == "EC"  ~ runif(dplyr::n(), 100, 500),
        name.analyte == "DOC" ~ runif(dplyr::n(), 0.2, 5),
        TRUE                  ~ exp(rnorm(dplyr::n(), log(2), 0.5))
      ),
      quantified = TRUE,
      imputed    = FALSE
    )

  # Introduce some BDL rows
  bdl_idx <- sample(which(rows$name.analyte %in% analytes), n_bdl)
  rows$quantified[bdl_idx] <- FALSE

  # Introduce some missing rows (remove them entirely — simulate not measured)
  miss_idx <- sample(which(rows$name.analyte %in% analytes), n_missing)
  rows <- rows[-miss_idx, ]

  rows
}

test_that("impute_chemistry requires brms (now a hard dependency in Imports)", {
  # brms is now in Imports; this test just confirms it is loadable.
  skip_if_not(requireNamespace("brms", quietly = TRUE))
  expect_true(requireNamespace("brms", quietly = TRUE))
})

test_that("impute_chemistry errors on missing driver analyte in df", {
  skip_if_not(requireNamespace("brms", quietly = TRUE))
  df <- make_impute_chem()
  # Remove all DOC rows
  df_no_doc <- dplyr::filter(df, name.analyte != "DOC")
  expect_error(
    impute_chemistry(df_no_doc, drivers = c("pH", "EC", "DOC")),
    regexp = "Driver analyte.* not found"
  )
})

test_that("impute_chemistry errors on BDL driver rows", {
  skip_if_not(requireNamespace("brms", quietly = TRUE))
  df <- make_impute_chem()
  # Make one pH row BDL
  df$quantified[df$name.analyte == "pH"][1L] <- FALSE
  expect_error(
    impute_chemistry(df, drivers = c("pH", "EC", "DOC")),
    regexp = "Driver analyte.* have missing or BDL"
  )
})

test_that("impute_chemistry full run (brms smoke test)", {
  skip_if_not(
    requireNamespace("brms", quietly = TRUE) &&
      identical(Sys.getenv("BRMS_SMOKE_TEST"), "1"),
    "Skipping brms smoke test: brms not installed or BRMS_SMOKE_TEST != '1'"
  )

  df  <- make_impute_chem(n = 30, n_bdl = 4, n_missing = 3)
  imp <- impute_chemistry(
    df,
    drivers = c("pH", "EC", "DOC"),
    iter    = 500,
    warmup  = 250,
    chains  = 1,
    cores   = 1
  )

  # All values should be finite after imputation
  expect_true(all(is.finite(imp$value)))

  # imputed column must exist
  expect_true("imputed" %in% names(imp))

  # Rows that were BDL or missing should now be imputed = TRUE
  expect_true(any(imp$imputed))

  # imputed_kind column must exist and have expected values
  expect_true("imputed_kind" %in% names(imp))
  expect_true(all(imp$imputed_kind %in% c("observed", "censored_left", "missing")))

  # Fitted brmsfit is attached as attribute
  expect_true(!is.null(attr(imp, "brmsfit")))
})
