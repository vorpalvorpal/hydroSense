## Tests for prepare_reference() — summary statistics and CI

library(testthat)
library(leachatetools)

## Mn and Hg have no normalisation formula in the bundled metadata, so the
## reference computation is a clean test of the summary statistics without
## interference from missing co-analytes.
make_ref <- function(n_per_analyte = 12, seed = 1) {
  set.seed(seed)
  analytes <- c("Mn", "Hg")
  tidyr::expand_grid(
    sample_id = paste0("r", seq_len(n_per_analyte)),
    analyte   = analytes
  ) |>
    dplyr::mutate(
      site_id  = "ref",
      datetime = as.Date("2024-01-01") + (match(sample_id,
        paste0("r", seq_len(n_per_analyte))) - 1L),
      value    = dplyr::case_when(
        analyte == "Mn" ~ exp(rnorm(dplyr::n(), log(1.5), 0.4)),
        analyte == "Hg" ~ exp(rnorm(dplyr::n(), log(0.5), 0.3))
      ),
      detected = TRUE
    )
}

test_that("default summary is geom_mean", {
  prep <- prepare_reference(make_ref(), conc_units = "ug/L")
  expect_equal(prep$summary, "geom_mean")
})

test_that("ref_table has the expected columns including n_obs", {
  prep <- prepare_reference(make_ref(), conc_units = "ug/L")
  expect_true(all(c("analyte", "ref_norm", "n_obs") %in% names(prep$ref_table)))
  expect_true(all(prep$ref_table$n_obs > 0L))
})

test_that("geom_mean ≤ arith_mean for positive concentrations (AM-GM)", {
  ref <- make_ref()
  gm  <- prepare_reference(ref, summary = "geom_mean", conc_units = "ug/L")
  am  <- prepare_reference(ref, summary = "arith_mean", conc_units = "ug/L")
  gm_vals <- gm$ref_table$ref_norm[order(gm$ref_table$analyte)]
  am_vals <- am$ref_table$ref_norm[order(am$ref_table$analyte)]
  expect_true(all(gm_vals <= am_vals + 1e-9))
})

test_that("median ≤ p80 ≤ p90 ≤ p95", {
  ref  <- make_ref(n_per_analyte = 50, seed = 99)
  med  <- prepare_reference(ref, summary = "median",    conc_units = "ug/L")
  p80  <- prepare_reference(ref, summary = "p80",       conc_units = "ug/L")
  p90  <- prepare_reference(ref, summary = "p90",       conc_units = "ug/L")
  p95  <- prepare_reference(ref, summary = "p95",       conc_units = "ug/L")
  by_an <- function(p) p$ref_table$ref_norm[order(p$ref_table$analyte)]
  expect_true(all(by_an(med) <= by_an(p80) + 1e-9))
  expect_true(all(by_an(p80) <= by_an(p90) + 1e-9))
  expect_true(all(by_an(p90) <= by_an(p95) + 1e-9))
})

test_that("bootstrap_ci = TRUE adds ref_lower / ref_upper", {
  prep <- prepare_reference(make_ref(), bootstrap_ci = TRUE, n_boot = 200L,
                            conc_units = "ug/L")
  expect_true(all(c("ref_lower", "ref_upper") %in% names(prep$ref_table)))
  # CI brackets the point estimate
  expect_true(all(prep$ref_table$ref_lower <= prep$ref_table$ref_norm + 1e-9))
  expect_true(all(prep$ref_table$ref_upper >= prep$ref_table$ref_norm - 1e-9))
})

test_that("low-n warning fires when analyte has < 5 observations", {
  ref <- make_ref(n_per_analyte = 3)
  expect_warning(
    prepare_reference(ref, conc_units = "ug/L"),
    regexp = "< 5 reference observations"
  )
})

test_that("BDL observations contribute 0 to the summary", {
  # All Mn rows BDL → ref_norm = 0 (geom_mean of zeros via eps shift)
  ref <- make_ref()
  ref$detected[ref$analyte == "Mn"] <- FALSE
  prep <- suppressMessages(prepare_reference(ref, conc_units = "ug/L"))
  mn_ref <- prep$ref_table$ref_norm[prep$ref_table$analyte == "Mn"]
  expect_lt(mn_ref, 1e-6)
})

test_that("empty input yields a usable empty prepared_reference", {
  ref <- make_ref()[0, ]
  expect_warning(prep <- prepare_reference(ref), regexp = "No reference")
  expect_s3_class(prep, "prepared_reference")
  expect_equal(nrow(prep$ref_table), 0L)
})

test_that("print method runs without error", {
  prep <- prepare_reference(make_ref(), conc_units = "ug/L")
  expect_output(print(prep), regexp = "prepared_reference")
})
