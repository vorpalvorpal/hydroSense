## Tests for prescreen_analytes()

library(testthat)
library(leachatetools)

make_chem <- function() {
  # expand_grid varies name.analyte slowest (outer) so rows are grouped by sample.
  # Use case_when on name.analyte to set per-analyte detection rates:
  #   Cu:     19/20 quantified
  #   Zn:     19/20 quantified
  #   Rarely:  1/20 quantified (5 %)
  tidyr::expand_grid(
    uuid.sample  = paste0("s", 1:20),
    name.analyte = c("Cu", "Zn", "Rarely")
  ) |>
    dplyr::mutate(
      uuid.feature = "f1",
      value        = runif(dplyr::n(), 0.1, 10),
      quantified   = dplyr::case_when(
        name.analyte == "Cu"     ~ uuid.sample != "s20",
        name.analyte == "Zn"     ~ uuid.sample != "s20",
        name.analyte == "Rarely" ~ uuid.sample == "s1"    # only 1/20
      )
    )
}

test_that("prescreen_analytes returns vector of passing analytes", {
  df  <- make_chem()
  inc <- prescreen_analytes(df, k = 0.10)   # Rarely fails at 5 %
  expect_type(inc, "character")
  expect_true("Cu"  %in% inc)
  expect_true("Zn"  %in% inc)
  expect_false("Rarely" %in% inc)
})

test_that("excluded attribute captures dropped analytes", {
  df  <- make_chem()
  inc <- prescreen_analytes(df, k = 0.10)
  expect_equal(attr(inc, "excluded"), "Rarely")
})

test_that("return = 'table' gives a tibble with all expected columns", {
  df  <- make_chem()
  tbl <- prescreen_analytes(df, k = 0.10, return = "table")
  expect_s3_class(tbl, "tbl_df")
  expect_true(all(c("name.analyte","n_samples","n_quantified","detect_freq","included") %in% names(tbl)))
  expect_equal(tbl$included[tbl$name.analyte == "Rarely"], FALSE)
})

test_that("k = 0 includes all analytes", {
  df  <- make_chem()
  inc <- prescreen_analytes(df, k = 0)
  expect_setequal(inc, c("Cu", "Zn", "Rarely"))
})

test_that("k = 1 excludes analytes that are not 100 % detected", {
  df  <- make_chem()
  inc <- prescreen_analytes(df, k = 1)
  expect_length(inc, 0L)
})

test_that("group_by_feature requires uuid.feature column", {
  df_no_feat <- tibble::tibble(
    name.analyte = "Cu", value = 1, quantified = TRUE, uuid.sample = "s1"
  )
  expect_error(
    prescreen_analytes(df_no_feat, group_by_feature = TRUE),
    regexp = "uuid.feature"
  )
})
