## Regression tests for the single-substance SSD/PAF engine. These exercise the
## bundled ANZG observation data (no external "guideline data" folder needed),
## so they run offline and pin the numeric behaviour of ssd_hc50()/ssd_pct().

library(testthat)
library(leachatetools)

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
  expect_equal(ssd_pct("Cu", hc50), 50, tolerance = 1e-2)
})

test_that("ssd_pct() is monotonic increasing in concentration", {
  lo <- ssd_pct("Zn", 5)
  hi <- ssd_pct("Zn", 30)
  expect_gt(hi, lo)
  expect_gte(lo, 0)
  expect_lte(hi, 100)
})

test_that("ssd_paf() rejects non-positive concentrations", {
  expect_error(ssd_paf("Cu", 0))
  expect_error(ssd_paf("Cu", -1))
})

test_that("ssd_paf() returns the documented list shape", {
  res <- ssd_paf("Cu", 4)
  expect_type(res, "list")
  expect_named(res, c("analyte", "conc_ug_L", "method", "pct",
                      "lower", "upper", "note"))
  expect_true(is.numeric(res$pct) && res$pct >= 0 && res$pct <= 100)
})
