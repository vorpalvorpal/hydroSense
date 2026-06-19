## Regression tests for the single-substance SSD/PAF engine. These exercise the
## bundled ANZG observation data (no external "guideline data" folder needed),
## so they run offline and pin the numeric behaviour of ssd_hc50()/ssd_pct().

library(testthat)
library(hydroSense)

test_that("ssd_hc50() returns the ANZG index value for total ammonia-N", {
  # NH3-N SSD is expressed as total ammonia-N at the pH 7.0 / 20 degC index
  # condition; its HC50 is a fixed property of the bundled dataset.
  expect_equal(ssd_hc50("NH3-N"), 9321, tolerance = 1e-3)
})

test_that("ssd_hc50() returns a stable copper HC50", {
  expect_equal(ssd_hc50("Cu"), 4.2258, tolerance = 1e-3)
})

test_that("ssd_hc50() returns NA for analytes with no SSD", {
  expect_true(is.na(ssd_hc50("Notreal")))
})

test_that("PAF at the HC50 is 50%", {
  hc50 <- ssd_hc50("Cu")
  expect_equal(ssd_pct("Cu", hc50, conc_units = "ug/L"), 50, tolerance = 1e-2)
})

test_that("ssd_pct() is monotonic increasing in concentration", {
  lo <- ssd_pct("Zn", 5,  conc_units = "ug/L")
  hi <- ssd_pct("Zn", 30, conc_units = "ug/L")
  expect_gt(hi, lo)
  expect_gte(lo, 0)
  expect_lte(hi, 100)
})

test_that("ssd_paf() rejects non-positive concentrations", {
  expect_error(ssd_paf("Cu", 0,  conc_units = "ug/L"))
  expect_error(ssd_paf("Cu", -1, conc_units = "ug/L"))
})

test_that("ssd_paf() returns the documented list shape", {
  res <- ssd_paf("Cu", 4, conc_units = "ug/L")
  expect_type(res, "list")
  expect_named(res, c("analyte", "conc_ug_L", "method", "pct",
                      "lower", "upper", "note"))
  expect_true(is.numeric(res$pct) && res$pct >= 0 && res$pct <= 100)
})

## --- NO3-N probabilistic hardness weighting (§10a) -------------------------

test_that(".no3_weights are a valid distribution and track hardness", {
  expect_equal(sum(hydroSense:::.no3_weights(90)), 1)
  # Soft water -> soft class dominates; hard water -> hard class dominates.
  expect_gt(hydroSense:::.no3_weights(15)["soft"],  0.9)
  expect_gt(hydroSense:::.no3_weights(500)["hard"], 0.9)
  # On a boundary (150) mod and hard split evenly.
  w150 <- hydroSense:::.no3_weights(150)
  expect_equal(unname(w150["mod"]), 0.5, tolerance = 1e-6)
  expect_equal(unname(w150["hard"]), 0.5, tolerance = 1e-6)
  # Missing hardness -> all NA.
  expect_true(all(is.na(hydroSense:::.no3_weights(NA))))
})

test_that("ssd_paf NO3-N blend equals the weighted mean of class PAFs", {
  conc <- 50000
  blend <- ssd_paf("NO3-N", conc, conc_units = "ug/L",
                   hardness = 150, hardness_units = "mg/L")
  mod   <- ssd_paf("NO3-N_mod",  conc, conc_units = "ug/L")$pct
  hard  <- ssd_paf("NO3-N_hard", conc, conc_units = "ug/L")$pct
  expect_equal(blend$pct, 0.5 * mod + 0.5 * hard, tolerance = 1e-6)
  expect_match(blend$note, "hardness blend")
})

test_that("ssd_paf NO3-N is continuous across the hardness class boundary", {
  conc <- 50000
  lo <- ssd_paf("NO3-N", conc, conc_units = "ug/L",
                hardness = 149, hardness_units = "mg/L")$pct
  hi <- ssd_paf("NO3-N", conc, conc_units = "ug/L",
                hardness = 151, hardness_units = "mg/L")$pct
  # The smooth blend changes only slightly across 150 (the hard cutoff would
  # jump by tens of percent between the moderate and hard SSDs).
  expect_lt(abs(lo - hi), 10)
})

test_that("ssd_paf NO3-N without hardness falls back to the soft class", {
  expect_warning(res <- ssd_paf("NO3-N", 50000, conc_units = "ug/L"), "hardness")
  expect_equal(res$analyte, "NO3-N_soft")
})
