## Tests for .parse_normalisation_formula() and .apply_normalisation()

library(testthat)
library(leachatetools)

test_that("NULL/empty formula is identity", {
  parsed <- leachatetools:::.parse_normalisation_formula("")
  expect_null(parsed)
  expect_equal(leachatetools:::.apply_normalisation(NULL, 5.0), 5.0)
})

test_that("NA formula is identity", {
  parsed <- leachatetools:::.parse_normalisation_formula(NA_character_)
  expect_null(parsed)
})

test_that("simple expression is evaluated correctly", {
  parsed <- leachatetools:::.parse_normalisation_formula("C * 2")
  expect_equal(leachatetools:::.apply_normalisation(parsed, 3.0), 6.0)
})

test_that("co-analyte variables are available in expression", {
  parsed <- leachatetools:::.parse_normalisation_formula("C / DOC")
  result <- leachatetools:::.apply_normalisation(parsed, 10.0, c(DOC = 2.0))
  expect_equal(result, 5.0)
})

test_that("missing co-analyte returns NA", {
  parsed <- leachatetools:::.parse_normalisation_formula("C / DOC")
  result <- leachatetools:::.apply_normalisation(parsed, 10.0, c(pH = 7))
  expect_true(is.na(result))
})

test_that("malformed formula raises informative error", {
  expect_error(
    leachatetools:::.parse_normalisation_formula("C +** 2"),
    class = "rlang_error"
  )
})

test_that("parsed formulas are cached (same object on second parse)", {
  f1 <- leachatetools:::.parse_normalisation_formula("C + 1")
  f2 <- leachatetools:::.parse_normalisation_formula("C + 1")
  expect_true(identical(f1, f2))
})
