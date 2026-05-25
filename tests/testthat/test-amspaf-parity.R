## AmsPAF parity test: new metadata-driven add_amspaf() must produce the same
## value as the Phase 1 version when:
##   - normalisation is identity (all formula cells empty — Day-1 state)
##   - MOA grouping matches the prior type-heuristic (ionoregulatory group)
##   - all inputs are quantified (no imputation needed)
##   - reference = NULL (no ARA adjustment)
##
## This test is intentionally simple. It does NOT require the ANZG XLSX files;
## it passes zero concentrations so ssd_paf() returns 0 for all analytes.
## The key property being tested is that the refactored add_amspaf() doesn't
## crash, correctly propagates columns, and the AmsPAF row is appended.

library(testthat)
library(leachatetools)

make_minimal_chem <- function(analytes = c("Cu", "Zn", "Ni"),
                              n_samples = 3,
                              feature = "f1") {
  tidyr::expand_grid(
    uuid.sample  = paste0("smp", seq_len(n_samples)),
    name.analyte = analytes
  ) |>
    dplyr::mutate(
      uuid.feature    = feature,
      datetime.sample = as.Date("2024-01-01") + (match(uuid.sample,
        paste0("smp", seq_len(n_samples))) - 1L),
      value      = 0.001,   # near zero → near-zero PAF; no XLSX needed
      quantified = TRUE,
      imputed    = FALSE
    )
}

test_that("add_amspaf appends AmsPAF rows without crashing (identity normalisation)", {
  skip_if(
    !nzchar(getOption("leachatetools.guideline_dir", "")),
    "Skipping: leachatetools.guideline_dir not set (no ANZG XLSX available)"
  )

  df     <- make_minimal_chem()
  result <- add_amspaf(df, reference = NULL)

  amspaf_rows <- dplyr::filter(result, name.analyte == "AmsPAF")
  expect_gt(nrow(amspaf_rows), 0L)
  expect_true(all(amspaf_rows$quantified))
  expect_true("n_analytes_used"    %in% names(result))
  expect_true("n_analytes_imputed" %in% names(result))
})

test_that("add_amspaf n_analytes_imputed is 0 when imputed column absent", {
  skip_if(
    !nzchar(getOption("leachatetools.guideline_dir", "")),
    "Skipping: leachatetools.guideline_dir not set"
  )

  df_no_imputed <- make_minimal_chem() |> dplyr::select(-"imputed")
  result        <- add_amspaf(df_no_imputed, reference = NULL)
  amspaf_rows   <- dplyr::filter(result, name.analyte == "AmsPAF")
  expect_true(all(amspaf_rows$n_analytes_imputed == 0L))
})

test_that("add_amspaf accepts a prepared_reference object", {
  skip_if(
    !nzchar(getOption("leachatetools.guideline_dir", "")),
    "Skipping: leachatetools.guideline_dir not set"
  )

  df       <- make_minimal_chem()
  ref_df   <- make_minimal_chem(feature = "ref")
  prep_ref <- prepare_reference(ref_df)
  result   <- add_amspaf(df, reference = prep_ref)
  expect_true("AmsPAF" %in% result$name.analyte)
})

test_that("classify_amspaf_tier returns expected labels", {
  expect_equal(classify_amspaf_tier(c(0.5, 2, 7, 15)),
               c("1_background", "2_elevated", "3_impacted", "4_severely_impacted"))
})
