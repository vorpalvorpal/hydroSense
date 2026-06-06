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

## ── Ammonia pH/temperature correction direction ──────────────────────────────
## Toxicity tracks the un-ionised NH3 fraction, which rises with pH/temperature.
## The correction maps total ammonia-N TO the pH 7.0 / 20 °C reference, so a
## high-pH sample must normalise UP and a low-pH sample DOWN. An inverted
## formula (C * f_ref / f_sample) understates ammonia risk at high pH — these
## tests fail loudly if that regression is ever reintroduced.

test_that("correct_ammonia_ph_temp is identity at the reference condition", {
  expect_equal(correct_ammonia_ph_temp(900, conc_units = "ug/L",
                                       pH = 7.0, temperature_C = 20), 900)
})

test_that("high pH normalises ammonia upward, low pH downward", {
  hi <- correct_ammonia_ph_temp(900, conc_units = "ug/L", pH = 8.5, temperature_C = 20)
  lo <- correct_ammonia_ph_temp(900, conc_units = "ug/L", pH = 6.5, temperature_C = 20)
  expect_gt(hi, 900)   # more un-ionised NH3 than reference → more toxic
  expect_lt(lo, 900)   # less un-ionised NH3 than reference → less toxic
  expect_gt(hi, lo)
})

test_that("higher temperature normalises ammonia upward", {
  warm <- correct_ammonia_ph_temp(900, conc_units = "ug/L", pH = 7.0, temperature_C = 30)
  cool <- correct_ammonia_ph_temp(900, conc_units = "ug/L", pH = 7.0, temperature_C = 10)
  expect_gt(warm, 900)
  expect_lt(cool, 900)
})

test_that("the metadata NH3-N normalisation_formula matches the exported helper", {
  meta <- leachatetools:::.load_analyte_metadata()
  row  <- meta[meta$analyte == "NH3-N", , drop = FALSE]
  skip_if(nrow(row) != 1L, "NH3-N row not found in metadata")

  parsed <- leachatetools:::.parse_normalisation_formula(row$normalisation_formula)
  for (cond in list(c(pH = 6.5, t = 12), c(pH = 7.0, t = 20), c(pH = 8.7, t = 28))) {
    via_meta <- leachatetools:::.apply_normalisation(
      parsed, 900, c(pH = cond[["pH"]], temperature = cond[["t"]])
    )
    via_help <- correct_ammonia_ph_temp(900, conc_units = "ug/L",
                                        pH = cond[["pH"]], temperature_C = cond[["t"]])
    expect_equal(via_meta, via_help, tolerance = 1e-9)
  }
})
