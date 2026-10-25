## ============================================================================
## Self-contained ion chemistry for to_meq() / add_lmf()
## ============================================================================
##
## Contract (TDD): the package supplies valence and molar mass for the LMF
## panel ions itself. to_meq() and add_lmf() accept sampleTidy-shaped input —
## columns sample_id, site_id, analyte, value, detected, datetime and a
## character `units` column (e.g. "mg/L") — with NO valence.analyte /
## atomic_mass.analyte columns. A `units.analyte` column remains accepted, and
## caller-supplied valence/atomic-mass columns keep working (test-to-meq.R).
##
## Expected meq/L use 1 molar-mass-worth of each species, so meq/L == |z|.
## Tolerance 2e-3 absorbs differences in the molar masses the package chooses
## (e.g. 22.99 vs 22.98977 for Na).

library(testthat)
library(hydroSense)

## One row per LMF panel ion: value = molar mass in mg/L, expected meq/L.
.ion_expect <- function() {
  tibble::tribble(
    ~analyte,       ~value,   ~meq,
    "Na",            22.990,   1,
    "K",             39.098,   1,
    "Ca",            40.078,   2,
    "Mg",            24.305,   2,
    "Cl",            35.453,   1,
    "F",             18.998,   1,
    "SO4²⁻", 96.06,  2,
    "NH3-N",         14.007,   1,   # as N
    "NO3-N",         14.007,   1,   # as N
    "NO2-N",         14.007,   1,   # as N
    "CO3-CaCO3",    100.087,   2,   # as CaCO3
    "HCO3-CaCO3",   100.087,   2    # as CaCO3
  )
}

.st_rows <- function(analyte, value, units = "mg/L") {
  tibble::tibble(
    sample_id = "s1", site_id = "site", analyte = analyte, value = value,
    detected = TRUE, datetime = .pi_utc("2022-06-01"), units = units
  )
}

test_that("to_meq converts every LMF panel ion with only a `units` column", {
  ex  <- .ion_expect()
  df  <- .st_rows(ex$analyte, ex$value)
  out <- to_meq(df)

  conv <- out[out$analyte %in% paste0(ex$analyte, "_"), ]
  expect_equal(nrow(conv), nrow(ex))
  got <- conv$value[match(paste0(ex$analyte, "_"), conv$analyte)]
  expect_equal(got, ex$meq, tolerance = 2e-3)

  ## Originals preserved unchanged.
  orig <- out[!out$analyte %in% paste0(ex$analyte, "_"), ]
  expect_equal(orig$analyte, df$analyte)
  expect_equal(orig$value, df$value)
})

test_that("to_meq honours the `units` string (ug/L and µg/L)", {
  out_u <- to_meq(.st_rows("Ca", 40078, units = "ug/L"))
  expect_equal(out_u$value[out_u$analyte == "Ca_"], 2, tolerance = 2e-3)

  out_mu <- to_meq(.st_rows("Ca", 40078, units = "µg/L"))
  expect_equal(out_mu$value[out_mu$analyte == "Ca_"], 2, tolerance = 2e-3)
})

test_that("to_meq leaves non-ionic analytes unconverted with only `units`", {
  df  <- .st_rows(c("Ca", "pH", "EC"), c(40.078, 7.2, 300),
                  units = c("mg/L", "pH units", "uS/cm"))
  out <- to_meq(df)
  expect_true("Ca_" %in% out$analyte)
  expect_false(any(c("pH_", "EC_") %in% out$analyte))
  expect_equal(nrow(out), 4L)   # 3 originals + 1 converted
})

test_that("to_meq still requires a units column of some kind", {
  bad <- tibble::tibble(analyte = "Ca", value = 40)
  expect_error(to_meq(bad))
})

test_that("add_lmf works on sampleTidy-shaped input with no valence/mass columns", {
  set.seed(101)
  ref <- .pi_samples(rep("2022-06-01", 12), 0, "ref", "refsite")
  lea <- .pi_samples(rep("2022-06-01", 12), 1, "lea", "leasite")
  fvec <- c(0, 0.25, 0.5, 0.75, 1)
  set.seed(7)
  q <- dplyr::bind_rows(lapply(seq_along(fvec), function(i) {
    .pi_samples("2022-06-01", fvec[i], paste0("q", i), "query")
  }))

  native <- .pi_lmf(add_lmf(q, leachate_data = lea, reference_data = ref))

  q_st   <- .pi_as_sampletidy(q)
  ref_st <- .pi_as_sampletidy(ref)
  lea_st <- .pi_as_sampletidy(lea)
  expect_false(any(c("valence.analyte", "atomic_mass.analyte", "units.analyte") %in%
                   c(names(q_st), names(ref_st), names(lea_st))))

  out_st <- add_lmf(q_st, leachate_data = lea_st, reference_data = ref_st)
  st     <- .pi_lmf(out_st)

  expect_equal(nrow(st), length(fvec))
  expect_true(all(is.na(st$lmf_reason)))
  ## LMF is invariant to the per-ion meq factor, so package-supplied molar
  ## masses must reproduce the caller-supplied-column result.
  expect_equal(st$value, native$value, tolerance = 1e-6)
  expect_equal(st$sigma_lmf, native$sigma_lmf, tolerance = 1e-6)
  expect_equal(st$value, 100 * fvec, tolerance = 3)

  ## The caller's rows come back unchanged (names, values, units).
  orig <- out_st[out_st$analyte != "LMF", names(q_st)]
  orig <- orig[order(orig$sample_id, orig$analyte), ]
  q_st <- q_st[order(q_st$sample_id, q_st$analyte), ]
  expect_equal(as.data.frame(orig), as.data.frame(q_st), ignore_attr = TRUE)
})

test_that("add_lmf works on sampleTidy-shaped input in trailing mode", {
  ref <- .pi_as_sampletidy(.pi_samples(.pi_monthly("2017-01-01", "2022-05-01"), 0, "ref", "refsite", seed = 1))
  lea <- .pi_as_sampletidy(.pi_samples(.pi_monthly("2017-01-01", "2022-05-01"), 1, "lea", "leasite", seed = 2))
  q   <- .pi_as_sampletidy(.pi_samples("2022-06-01", 0.5, "q", "query", seed = 3))
  out <- add_lmf(q, leachate_data = lea, reference_data = ref, calibration = "trailing")
  lmf <- .pi_lmf(out)
  expect_equal(nrow(lmf), 1L)
  expect_equal(lmf$value, 50, tolerance = 5)
})
