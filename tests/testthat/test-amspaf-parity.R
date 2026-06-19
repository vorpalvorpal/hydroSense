## AmsPAF unit tests: metadata-driven add_amspaf() with bundled SSD data.
##
## These tests do NOT require the ANZG XLSX files — the bundled
## anzg_xlsx_observations.csv ships with the package and is used automatically
## when guideline_dir is not set.
##
## Co-analytes (pH, DOC, Ca, Mg, hardness) are included in make_minimal_chem()
## because the normalisation formulas for Cu (DOC), Zn (pH/hardness/DOC), and
## Ni (pH/DOC/Ca/Mg) require them.

library(testthat)
library(hydroSense)

make_minimal_chem <- function(analytes = c("Cu", "Zn", "Ni"),
                              n_samples = 3,
                              feature = "f1") {
  ## Include co-analytes required by normalisation formulas for Cu/Zn/Ni.
  all_analytes <- union(analytes, c("pH", "DOC", "Ca", "Mg", "hardness"))

  tidyr::expand_grid(
    sample_id = paste0("smp", seq_len(n_samples)),
    analyte   = all_analytes
  ) |>
    dplyr::mutate(
      site_id  = feature,
      datetime = as.Date("2024-01-01") + (match(sample_id,
        paste0("smp", seq_len(n_samples))) - 1L),
      value = dplyr::case_when(
        analyte == "pH"       ~ 7.5,
        analyte == "DOC"      ~ 2.0,
        analyte == "Ca"       ~ 6.0,
        analyte == "Mg"       ~ 4.0,
        analyte == "hardness" ~ 30.0,
        TRUE                  ~ 0.001   # near zero → near-zero PAF for metals
      ),
      detected = TRUE,
      imputed  = FALSE
    )
}

test_that("add_amspaf appends AmsPAF rows without crashing", {
  df     <- make_minimal_chem()
  result <- add_amspaf(df, reference = NULL, conc_units = "ug/L")

  amspaf_rows <- dplyr::filter(result, analyte == "AmsPAF")
  expect_gt(nrow(amspaf_rows), 0L)
  expect_true(all(amspaf_rows$detected))
  expect_true("n_analytes_used"    %in% names(result))
  expect_true("n_analytes_imputed" %in% names(result))
})

test_that("add_amspaf n_analytes_imputed is 0 when imputed column absent", {
  df_no_imputed <- make_minimal_chem() |> dplyr::select(-"imputed")
  result        <- add_amspaf(df_no_imputed, reference = NULL, conc_units = "ug/L")
  amspaf_rows   <- dplyr::filter(result, analyte == "AmsPAF")
  expect_true(all(amspaf_rows$n_analytes_imputed == 0L))
})

test_that("add_amspaf accepts a prepared_reference object", {
  df       <- make_minimal_chem()
  ref_df   <- make_minimal_chem(feature = "ref")
  prep_ref <- prepare_reference(ref_df, conc_units = "ug/L")
  result   <- add_amspaf(df, reference = prep_ref, conc_units = "ug/L")
  expect_true("AmsPAF" %in% result$analyte)
})

