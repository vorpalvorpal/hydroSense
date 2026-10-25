## ============================================================================
## mspaf_pipeline() end to end on plain data frames (no DB, no Stan)
## ============================================================================
##
## Contract: mspaf_pipeline() with its defaults (impute = TRUE, the "marginal"
## mgcv imputation method, chronic = TRUE) runs on plain long-format data
## frames and returns the chronic msPAF frame — one row per (focal_date x
## site), analyte "msPAF", a numeric `value` column — with the fitted models
## as attributes.
##
## The first spec locks today's behaviour for package-named input carrying a
## units.analyte column. The second (TDD) requires the same result for
## sampleTidy-shaped input: sampleTidy analyte names (Temperature,
## Hardness-total-CaCO3) and a character `units` column instead of
## units.analyte (toxicants in "µg/L").
##
## Kept small (2 target sites x 6 samples; 12 reference samples, the minimum
## the imputation fit accepts is 10) so each run takes a few seconds.

library(testthat)
library(hydroSense)

.plain_inputs <- function() {
  list(
    target    = .pi_mspaf_chem(c("ds1", "ds2"), n = 6L, seed = 2L, mult = 5),
    reference = .pi_mspaf_chem("ref", n = 12L, seed = 1L),
    hydro     = .pi_hydro()
  )
}

## sampleTidy shape: stored names + `units` (no units.analyte).
.to_sampletidy_mspaf <- function(df) {
  df <- .pi_mspaf_sampletidy_names(df)
  natural <- c(pH = "pH units", EC = "µS/cm", DOC = "mg/L", Ca = "mg/L",
               Mg = "mg/L", "Hardness-total-CaCO3" = "mg/L",
               Temperature = "°C")
  units <- ifelse(df$analyte %in% c("Cu", "Zn", "Ni", "NH3-N"), "µg/L",
                  unname(natural[df$analyte]))
  tibble::tibble(
    sample_id = df$sample_id, site_id = df$site_id, analyte = df$analyte,
    value = df$value, detected = df$detected, datetime = df$datetime,
    units = units
  )
}

.run_pipeline <- function(inp) {
  suppressMessages(suppressWarnings(
    mspaf_pipeline(inp$target, inp$reference, hydro = inp$hydro)
  ))
}

.expect_chronic_frame <- function(out, sites) {
  expect_s3_class(out, "data.frame")
  expect_true(all(c("focal_date", "site_id", "analyte", "value") %in% names(out)))
  expect_s3_class(out$focal_date, "Date")
  expect_type(out$value, "double")
  expect_true(all(out$analyte == "msPAF"))
  expect_setequal(unique(out$site_id), sites)

  ## One row per (focal_date x site).
  expect_false(anyDuplicated(out[, c("focal_date", "site_id")]) > 0)
  ## Daily focal grid per site (focal_by = "day").
  for (s in sites) {
    fd <- sort(out$focal_date[out$site_id == s])
    expect_gt(length(fd), 100L)
    expect_true(all(as.numeric(diff(fd)) == 1))
  }

  ## msPAF is a percentage; most focal dates carry a value.
  v <- out$value[!is.na(out$value)]
  expect_gt(length(v), 0.9 * nrow(out))
  expect_true(all(is.finite(v) & v >= 0 & v <= 100))

  im <- attr(out, "imputation_model")
  expect_false(is.null(im))
  expect_identical(im$impute_method, "marginal")
  expect_false(is.null(attr(out, "reference_model")))
}

test_that("mspaf_pipeline runs end to end on plain data frames (package names)", {
  inp <- .plain_inputs()
  out <- .run_pipeline(inp)
  .expect_chronic_frame(out, c("ds1", "ds2"))
})

test_that("mspaf_pipeline gives the same chronic msPAF on sampleTidy-shaped input", {
  inp <- .plain_inputs()
  pkg <- .run_pipeline(inp)

  st_inp <- list(
    target    = .to_sampletidy_mspaf(inp$target),
    reference = .to_sampletidy_mspaf(inp$reference),
    hydro     = inp$hydro
  )
  expect_false("units.analyte" %in% names(st_inp$target))
  st <- .run_pipeline(st_inp)
  .expect_chronic_frame(st, c("ds1", "ds2"))

  key <- function(d) d[order(d$site_id, d$focal_date), c("focal_date", "site_id", "value")]
  expect_equal(as.data.frame(key(st)), as.data.frame(key(pkg)),
               tolerance = 1e-6, ignore_attr = TRUE)
})
