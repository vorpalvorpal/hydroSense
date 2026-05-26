## Tests for prescreen_analytes()

library(testthat)
library(leachatetools)

make_chem <- function() {
  # expand_grid varies analyte slowest (outer) so rows are grouped by sample.
  # Use case_when on analyte to set per-analyte detection rates:
  #   Cu:     19/20 detected
  #   Zn:     19/20 detected
  #   Rarely:  1/20 detected (5 %)
  tidyr::expand_grid(
    sample_id = paste0("s", 1:20),
    analyte   = c("Cu", "Zn", "Rarely")
  ) |>
    dplyr::mutate(
      site_id  = "f1",
      value    = runif(dplyr::n(), 0.1, 10),
      detected = dplyr::case_when(
        analyte == "Cu"     ~ sample_id != "s20",
        analyte == "Zn"     ~ sample_id != "s20",
        analyte == "Rarely" ~ sample_id == "s1"    # only 1/20
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
  expect_true(all(c("analyte","n_samples","n_detected","detect_freq","included") %in% names(tbl)))
  expect_equal(tbl$included[tbl$analyte == "Rarely"], FALSE)
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

test_that("group_by_feature requires site_id column", {
  df_no_feat <- tibble::tibble(
    analyte = "Cu", value = 1, detected = TRUE, sample_id = "s1"
  )
  expect_error(
    prescreen_analytes(df_no_feat, group_by_feature = TRUE),
    regexp = "site_id"
  )
})

test_that("return = 'table' includes protected column", {
  df  <- make_chem()
  tbl <- prescreen_analytes(df, k = 0.10, return = "table")
  expect_true("protected" %in% names(tbl))
})

test_that("co-analyte drivers are protected from exclusion at k = 1", {
  # pH and DOC are listed in coanalytes_required for NH3-N / Cu / etc.
  # They should be kept even when never detected.
  df <- tidyr::expand_grid(
    sample_id = paste0("s", 1:10),
    analyte   = c("Cu", "pH", "DOC")
  ) |>
    dplyr::mutate(
      site_id  = "f1",
      value    = 1,
      detected = analyte == "Cu"  # pH and DOC never detected
    )
  inc <- prescreen_analytes(df, k = 1)   # k=1 → only 100 % detected pass
  # Cu has 100 % detection → passes on merit
  # pH and DOC have 0 % detection → would fail but are protected
  expect_true("pH"  %in% inc)
  expect_true("DOC" %in% inc)
})
