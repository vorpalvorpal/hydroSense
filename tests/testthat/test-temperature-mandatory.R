## Temperature is mandatory for any sample assessed for NH3-N: the un-ionised
## fraction normalisation is undefined without water temperature. add_amspaf()
## must fail loudly (not silently drop ammonia) when it is missing.

library(testthat)
library(leachatetools)

test_that(".assert_temperature_present passes when every ammonia sample has temperature", {
  df <- tibble::tibble(
    sample_id = c("S1", "S1", "S2", "S2"),
    site_id   = "A",
    analyte   = c("NH3-N", "temperature", "NH3-N", "temperature"),
    value     = c(500, 18, 300, 21)
  )
  expect_true(leachatetools:::.assert_temperature_present(df))
})

test_that(".assert_temperature_present errors when an ammonia sample lacks temperature", {
  df <- tibble::tibble(
    sample_id = c("S1", "S1", "S2"),
    site_id   = "A",
    analyte   = c("NH3-N", "temperature", "NH3-N"),  # S2 has no temperature
    value     = c(500, 18, 300)
  )
  expect_error(
    leachatetools:::.assert_temperature_present(df),
    "temperature"
  )
})

test_that("a present-but-NA temperature does not satisfy the requirement", {
  df <- tibble::tibble(
    sample_id = c("S1", "S1"),
    site_id   = "A",
    analyte   = c("NH3-N", "temperature"),
    value     = c(500, NA_real_)
  )
  expect_error(leachatetools:::.assert_temperature_present(df), "S1")
})

test_that("datasets without ammonia are unaffected", {
  df <- tibble::tibble(
    sample_id = c("S1", "S1"),
    site_id   = "A",
    analyte   = c("Cu", "Zn"),
    value     = c(5, 20)
  )
  expect_true(leachatetools:::.assert_temperature_present(df))
})

test_that("add_amspaf rejects an ammonia sample with no temperature by default", {
  df <- tibble::tibble(
    sample_id = "S1",
    site_id   = "A",
    analyte   = "NH3-N",
    value     = 500
  )
  expect_error(add_amspaf(df), "temperature")
})

test_that("add_amspaf(require_temperature = FALSE) bypasses the check", {
  df <- tibble::tibble(
    sample_id = "S1",
    site_id   = "A",
    analyte   = "NH3-N",
    value     = 500
  )
  # Should not raise the temperature error; it may warn/return for other
  # reasons (too few analytes), but must get past the temperature gate.
  expect_no_error(
    suppressWarnings(add_amspaf(df, require_temperature = FALSE,
                                conc_units    = "ug/L",
                                guideline_dir = getOption("leachatetools.guideline_dir")))
  )
})
