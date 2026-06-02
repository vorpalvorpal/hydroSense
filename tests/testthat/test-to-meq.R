## to_meq(): mass-concentration -> milliequivalent-per-litre conversion.
## meq/L = (mmol of substance per litre) * |valence|. Converted rows are
## appended with the analyte name suffixed "_"; originals are preserved.

library(testthat)
library(leachatetools)

base_df <- function() {
  tibble::tibble(
    analyte             = c("Ca", "Mg", "Na", "Cl", "SO4"),
    value               = c(40.08, 24.305, 22.99, 35.45, 96.06),  # mg/L
    units.analyte       = "mg/L",
    valence.analyte     = c(2, 2, 1, -1, -2),
    atomic_mass.analyte = c(40.08, 24.305, 22.99, 35.45, 96.06)   # g/mol
  )
}

test_that("conversion equals mmol/L * |valence|", {
  out <- to_meq(base_df())
  conv <- out[grepl("_$", out$analyte), ]
  # Each input is 1 molar-mass worth -> 1 mmol/L, so meq/L == |valence|.
  expect_equal(conv$value[conv$analyte == "Ca_"],  2)
  expect_equal(conv$value[conv$analyte == "Mg_"],  2)
  expect_equal(conv$value[conv$analyte == "Na_"],  1)
  expect_equal(conv$value[conv$analyte == "Cl_"],  1)
  expect_equal(conv$value[conv$analyte == "SO4_"], 2)
})

test_that("original rows are preserved unchanged and converted rows appended", {
  df  <- base_df()
  out <- to_meq(df)
  expect_equal(nrow(out), 2L * nrow(df))
  # Originals retained verbatim.
  orig <- out[!grepl("_$", out$analyte), ]
  expect_equal(orig$value, df$value)
})

test_that("input units other than mg/L convert correctly (ug/L)", {
  df <- tibble::tibble(analyte = "Ca", value = 40080, units.analyte = "ug/L",
                       valence.analyte = 2, atomic_mass.analyte = 40.08)
  out <- to_meq(df)
  expect_equal(out$value[out$analyte == "Ca_"], 2)
})

test_that("rows lacking valence or atomic mass are passed through unconverted", {
  df <- tibble::tibble(
    analyte             = c("Ca", "pH"),
    value               = c(40.08, 7.5),
    units.analyte       = c("mg/L", "pH"),
    valence.analyte     = c(2, NA),
    atomic_mass.analyte = c(40.08, NA)
  )
  out <- to_meq(df)
  # Only Ca gets a converted "Ca_" row; pH does not.
  expect_true("Ca_" %in% out$analyte)
  expect_false("pH_" %in% out$analyte)
  expect_equal(nrow(out), 3L)  # 2 originals + 1 converted
})

test_that("a frame with no convertible rows is returned unchanged", {
  df <- tibble::tibble(analyte = "pH", value = 7, units.analyte = "pH",
                       valence.analyte = NA_real_, atomic_mass.analyte = NA_real_)
  out <- to_meq(df)
  expect_equal(nrow(out), 1L)
  expect_identical(out$analyte, "pH")
})

test_that("missing required columns raise an informative error", {
  bad <- tibble::tibble(analyte = "Ca", value = 40)
  expect_error(to_meq(bad))
})
