## Tests for unit-resolution helpers and the units.analyte column path.
##
## Issue #6: adopt units package for dimensional correctness throughout.
## Covers:
##   - .resolve_to()             scalar/vector conversion and error paths
##   - .convert_df_tox_to_ugL()  units.analyte column and conc_units fallback
##   - add_amspaf()              units.analyte → same PAF as conc_units µg/L × 1000
##   - add_amspaf() / prepare_reference() error when no units information at all

library(testthat)
library(hydroSense)


# ── .resolve_to() ─────────────────────────────────────────────────────────────

test_that(".resolve_to() converts bare numeric + units_str to the target unit", {
  # 1 mg/L == 1000 ug/L
  val <- hydroSense:::.resolve_to(1, "ug/L", units_str = "mg/L")
  expect_equal(val, 1000)
})

test_that(".resolve_to() accepts a units object and auto-converts", {
  x   <- units::set_units(1, "mg/L")
  val <- hydroSense:::.resolve_to(x, "ug/L")
  expect_equal(val, 1000)
})

test_that(".resolve_to() is a no-op when source and target units match", {
  val <- hydroSense:::.resolve_to(500, "ug/L", units_str = "ug/L")
  expect_equal(val, 500)
})

test_that(".resolve_to() errors when bare numeric has no units_str", {
  expect_error(
    hydroSense:::.resolve_to(1, "ug/L"),
    regexp = "bare numeric"
  )
})


# ── .convert_df_tox_to_ugL() ─────────────────────────────────────────────────

make_tox_df <- function(value, unit, analyte = "Cu") {
  tibble::tibble(
    analyte       = analyte,
    value         = value,
    units.analyte = unit,
    detected      = TRUE
  )
}

test_that(".convert_df_tox_to_ugL() converts mg/L rows to µg/L via units.analyte", {
  df  <- make_tox_df(1, "mg/L")
  out <- hydroSense:::.convert_df_tox_to_ugL(df, ssd_analytes = "Cu")
  expect_equal(out$value, 1000)
})

test_that(".convert_df_tox_to_ugL() is a no-op when value is already in µg/L", {
  df  <- make_tox_df(500, "ug/L")
  out <- hydroSense:::.convert_df_tox_to_ugL(df, ssd_analytes = "Cu")
  expect_equal(out$value, 500)
})

test_that(".convert_df_tox_to_ugL() uses conc_units fallback when no units.analyte column", {
  df  <- tibble::tibble(analyte = "Cu", value = 1, detected = TRUE)
  out <- hydroSense:::.convert_df_tox_to_ugL(df, ssd_analytes = "Cu", conc_units = "mg/L")
  expect_equal(out$value, 1000)
})

test_that(".convert_df_tox_to_ugL() errors when units.analyte is NA for a tox row", {
  df <- make_tox_df(1, NA_character_)
  expect_error(
    hydroSense:::.convert_df_tox_to_ugL(df, ssd_analytes = "Cu"),
    regexp = "Missing.*units.analyte"
  )
})

test_that(".convert_df_tox_to_ugL() errors when no units.analyte column and no conc_units", {
  df <- tibble::tibble(analyte = "Cu", value = 1, detected = TRUE)
  expect_error(
    hydroSense:::.convert_df_tox_to_ugL(df, ssd_analytes = "Cu"),
    regexp = "units.analyte"
  )
})

test_that(".convert_df_tox_to_ugL() does not convert non-SSD-eligible rows", {
  df <- tibble::tibble(
    analyte       = c("Cu", "pH"),
    value         = c(1, 7.5),
    units.analyte = c("mg/L", "pH"),
    detected      = TRUE
  )
  out <- hydroSense:::.convert_df_tox_to_ugL(df, ssd_analytes = "Cu")
  expect_equal(out$value[out$analyte == "Cu"], 1000)    # converted
  expect_equal(out$value[out$analyte == "pH"],  7.5)    # untouched
})

test_that(".convert_df_tox_to_ugL() returns empty df unchanged", {
  df  <- tibble::tibble(analyte = character(), value = numeric(), detected = logical())
  out <- hydroSense:::.convert_df_tox_to_ugL(df, ssd_analytes = "Cu")
  expect_equal(nrow(out), 0L)
})


# ── add_amspaf() via units.analyte column ─────────────────────────────────────
#
# Build a one-sample chemistry frame twice:
#   Route A — toxicant values in mg/L; units.analyte column carries "mg/L"
#   Route B — same values converted to µg/L (×1000); conc_units = "ug/L"
#
# Both should produce the same AmsPAF because .convert_df_tox_to_ugL() converts
# route A to µg/L before the SSD lookup.

make_amspaf_chem <- function(cu_val) {
  analytes <- c("Cu", "Zn", "pH", "DOC", "Ca", "Mg", "hardness")
  tidyr::expand_grid(sample_id = "s1", analyte = analytes) |>
    dplyr::mutate(
      site_id  = "f1",
      datetime = as.Date("2024-01-01"),
      value    = dplyr::case_when(
        analyte == "pH"       ~ 7.5,
        analyte == "DOC"      ~ 2.0,
        analyte == "Ca"       ~ 6.0,
        analyte == "Mg"       ~ 4.0,
        analyte == "hardness" ~ 30.0,
        TRUE                  ~ cu_val
      ),
      detected = TRUE
    )
}

test_that("add_amspaf() via units.analyte mg/L gives same AmsPAF as conc_units ug/L × 1000", {
  # Route A: 0.001 mg/L toxicants, declared via units.analyte column
  df_a <- make_amspaf_chem(0.001) |>
    dplyr::mutate(units.analyte = dplyr::case_when(
      analyte %in% c("Cu", "Zn") ~ "mg/L",
      TRUE                        ~ "ug/L"   # co-analytes not SSD-eligible; ignored by converter
    ))
  paf_a <- add_amspaf(df_a, reference = NULL)   # no conc_units needed

  # Route B: 1 µg/L toxicants (= 0.001 mg/L × 1000), via conc_units
  df_b <- make_amspaf_chem(1)
  paf_b <- add_amspaf(df_b, reference = NULL, conc_units = "ug/L")

  amspaf_a <- dplyr::filter(paf_a, analyte == "AmsPAF")$value
  amspaf_b <- dplyr::filter(paf_b, analyte == "AmsPAF")$value

  expect_equal(amspaf_a, amspaf_b, tolerance = 1e-9)
})

test_that("add_amspaf() errors when no units.analyte column and no conc_units", {
  df <- make_amspaf_chem(1)
  expect_error(
    add_amspaf(df, reference = NULL),
    regexp = "units.analyte|conc_units"
  )
})


# ── prepare_reference() errors without units ──────────────────────────────────

test_that("prepare_reference() errors when no units.analyte column and no conc_units", {
  df <- tibble::tibble(
    sample_id = "r1", site_id = "ref",
    datetime  = as.Date("2024-01-01"),
    analyte   = "Cu",
    value     = 1,
    detected  = TRUE
  )
  expect_error(
    prepare_reference(df),
    regexp = "units.analyte|conc_units"
  )
})
