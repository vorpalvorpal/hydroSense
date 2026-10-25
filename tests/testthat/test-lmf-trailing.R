## ============================================================================
## add_lmf(calibration = "trailing"): per-sample trailing calibration window
## ============================================================================
##
## Contract (TDD, written before the implementation):
##
##   add_lmf() gains `calibration = c("window", "trailing")`. The default
##   "window" is the existing behaviour. With "trailing", EACH sample's LMF is
##   computed using only the leachate_data and reference_data rows whose
##   datetime lies in
##
##       [sample datetime - calibration_window_years, sample datetime]
##
##   inclusive at both ends. Everything that is calibrated from the end-members
##   (R, sigma_R, L, sigma_L, and the high-information ion classification) is
##   therefore a function of that sample's window alone, so new data dated
##   after a sample can never move that sample's LMF.
##
##   A sample whose window holds no usable reference or leachate data gets an
##   LMF row with value NA and an lmf_reason starting
##   "insufficient_calibration_data" (it does not abort the whole call).
##
##   "trailing" requires leachate_data and reference_data (data frames with a
##   datetime column).
##
## Window length: calibration_window_years years before the sample. The
## boundary tests below are written so that both calendar years
## (lubridate::years()) and 365.25-day years satisfy them.

library(testthat)
library(hydroSense)

## Mid-mixture query sample on a given day.
.tr_query <- function(day, f = 0.4, prefix = "q", seed = 21) {
  .pi_samples(day, f = f, prefix = prefix, site = "query", seed = seed)
}

## Generous precision gate so value is never NA'd by sigma in these specs.
.tr_lmf <- function(df, lea, ref, ...) {
  add_lmf(df, leachate_data = lea, reference_data = ref,
          calibration = "trailing", max_sigma_lsi = 1e6, ...)
}


## ---------------------------------------------------------------------------
## API
## ---------------------------------------------------------------------------

test_that("add_lmf has a calibration argument with choices window/trailing", {
  fm <- formals(add_lmf)
  expect_true("calibration" %in% names(fm))
  expect_identical(eval(fm$calibration), c("window", "trailing"))
})

test_that("calibration defaults to 'window' (existing behaviour)", {
  set.seed(101)
  ref <- .pi_samples(rep("2022-06-01", 12), 0, "ref", "refsite")
  lea <- .pi_samples(rep("2022-06-01", 12), 1, "lea", "leasite")
  q   <- .pi_samples(c("2022-06-01", "2022-06-01"), 0.5, "q", "query", seed = 7)

  out_default <- .pi_lmf(add_lmf(q, leachate_data = lea, reference_data = ref))
  out_window  <- .pi_lmf(add_lmf(q, leachate_data = lea, reference_data = ref,
                                 calibration = "window"))
  expect_identical(out_default$value, out_window$value)
  expect_identical(out_default$sigma_lmf, out_window$sigma_lmf)
})

test_that("an unknown calibration value is rejected", {
  set.seed(101)
  ref <- .pi_samples(rep("2022-06-01", 12), 0, "ref", "refsite")
  lea <- .pi_samples(rep("2022-06-01", 12), 1, "lea", "leasite")
  q   <- .tr_query("2022-06-01")
  expect_error(
    add_lmf(q, leachate_data = lea, reference_data = ref, calibration = "bogus"),
    "should be one of|[Mm]ust be element of"
  )
})


## ---------------------------------------------------------------------------
## (a) New data never move old values
## ---------------------------------------------------------------------------

test_that("appending later samples and end-member data leaves earlier LMF bit-identical", {
  ## Early end-members: 2015-01 .. 2019-12. Late end-members (2020-02 ..
  ## 2023-12) are deliberately DRIFTED, so if they leaked into q1's window its
  ## LMF would move.
  ref_early <- .pi_samples(.pi_monthly("2015-01-01", "2019-12-01"), 0, "refE", "refsite", seed = 11)
  lea_early <- .pi_samples(.pi_monthly("2015-01-01", "2019-12-01"), 1, "leaE", "leasite", seed = 12)
  ref_late  <- .pi_samples(.pi_monthly("2020-02-01", "2023-12-01"), 0, "refL", "refsite",
                           scale = 1.6, seed = 13)
  lea_late  <- .pi_samples(.pi_monthly("2020-02-01", "2023-12-01"), 1, "leaL", "leasite",
                           scale = 0.6, seed = 14)

  q1 <- .tr_query("2020-01-01", f = 0.4, prefix = "q1", seed = 21)
  q2 <- .tr_query("2022-06-01", f = 0.6, prefix = "q2", seed = 22)
  sid1 <- unique(q1$sample_id)
  sid2 <- unique(q2$sample_id)

  run_a <- .tr_lmf(q1, lea_early, ref_early)
  run_b <- .tr_lmf(dplyr::bind_rows(q1, q2),
                   dplyr::bind_rows(lea_early, lea_late),
                   dplyr::bind_rows(ref_early, ref_late))

  a1 <- .pi_lmf_for(run_a, sid1)
  b1 <- .pi_lmf_for(run_b, sid1)
  expect_false(is.na(a1$value))
  expect_identical(b1$value, a1$value)
  expect_identical(b1$sigma_lmf, a1$sigma_lmf)
  expect_identical(b1$lmf_naive, a1$lmf_naive)
  expect_identical(b1$n_ions_used, a1$n_ions_used)

  ## The appended later sample is itself quantified (its window has data).
  b2 <- .pi_lmf_for(run_b, sid2)
  expect_false(is.na(b2$value))
  expect_true(is.na(b2$lmf_reason))
})

test_that("the same late data DO move the earlier value under calibration = 'window'", {
  ## Discrimination check for the spec above: the drifted late data are
  ## material, so the bit-identity under "trailing" is not vacuous.
  ref_early <- .pi_samples(.pi_monthly("2015-01-01", "2019-12-01"), 0, "refE", "refsite", seed = 11)
  lea_early <- .pi_samples(.pi_monthly("2015-01-01", "2019-12-01"), 1, "leaE", "leasite", seed = 12)
  ref_late  <- .pi_samples(.pi_monthly("2020-02-01", "2023-12-01"), 0, "refL", "refsite",
                           scale = 1.6, seed = 13)
  lea_late  <- .pi_samples(.pi_monthly("2020-02-01", "2023-12-01"), 1, "leaL", "leasite",
                           scale = 0.6, seed = 14)
  q1   <- .tr_query("2020-01-01", f = 0.4, prefix = "q1", seed = 21)
  sid1 <- unique(q1$sample_id)

  win_a <- add_lmf(q1, leachate_data = lea_early, reference_data = ref_early,
                   calibration = "window", max_sigma_lsi = 1e6)
  win_b <- add_lmf(q1, leachate_data = dplyr::bind_rows(lea_early, lea_late),
                   reference_data = dplyr::bind_rows(ref_early, ref_late),
                   calibration = "window", max_sigma_lsi = 1e6)
  expect_false(isTRUE(all.equal(.pi_lmf_for(win_b, sid1)$value,
                                .pi_lmf_for(win_a, sid1)$value)))
})

test_that("a sample's trailing LMF does not depend on which other samples are in df", {
  ref <- .pi_samples(.pi_monthly("2015-01-01", "2023-12-01"), 0, "ref", "refsite", seed = 31)
  lea <- .pi_samples(.pi_monthly("2015-01-01", "2023-12-01"), 1, "lea", "leasite", seed = 32)
  q1  <- .tr_query("2021-03-01", f = 0.3, prefix = "q1", seed = 41)
  q2  <- .tr_query("2023-09-01", f = 0.7, prefix = "q2", seed = 42)
  sid1 <- unique(q1$sample_id)

  alone    <- .pi_lmf_for(.tr_lmf(q1, lea, ref), sid1)
  together <- .pi_lmf_for(.tr_lmf(dplyr::bind_rows(q2, q1), lea, ref), sid1)
  expect_identical(together$value, alone$value)
  expect_identical(together$sigma_lmf, alone$sigma_lmf)
})


## ---------------------------------------------------------------------------
## Window bounds: [t - calibration_window_years, t], inclusive
## ---------------------------------------------------------------------------

test_that("end-member rows outside [t - 5 y, t] are ignored", {
  ## Three segments; only the middle one lies in the 2022-06-01 sample's
  ## window. The out-of-window segments have very different compositions.
  seg_old <- function(f, p, s) .pi_samples(.pi_monthly("2010-01-01", "2017-05-01"), f, paste0(p, "O"), s, scale = 2.0, seed = 51)
  seg_in  <- function(f, p, s) .pi_samples(.pi_monthly("2017-06-01", "2022-06-01"), f, paste0(p, "I"), s, scale = 1.0, seed = 52)
  seg_new <- function(f, p, s) .pi_samples(.pi_monthly("2022-07-01", "2025-12-01"), f, paste0(p, "N"), s, scale = 0.5, seed = 53)

  ref_all <- dplyr::bind_rows(seg_old(0, "ref", "refsite"), seg_in(0, "ref", "refsite"), seg_new(0, "ref", "refsite"))
  lea_all <- dplyr::bind_rows(seg_old(1, "lea", "leasite"), seg_in(1, "lea", "leasite"), seg_new(1, "lea", "leasite"))
  ref_in  <- seg_in(0, "ref", "refsite")
  lea_in  <- seg_in(1, "lea", "leasite")

  q   <- .tr_query("2022-06-01", f = 0.5, seed = 61)
  sid <- unique(q$sample_id)

  full     <- .pi_lmf_for(.tr_lmf(q, lea_all, ref_all), sid)
  filtered <- .pi_lmf_for(.tr_lmf(q, lea_in, ref_in), sid)
  expect_false(is.na(full$value))
  expect_identical(full$value, filtered$value)
  expect_identical(full$sigma_lmf, filtered$sigma_lmf)
})

test_that("the window is inclusive at both ends and excludes the day outside each end", {
  q   <- .tr_query("2022-06-01", f = 0.5, seed = 61)
  sid <- unique(q$sample_id)

  ref_base <- .pi_samples(.pi_monthly("2017-07-01", "2022-05-01"), 0, "ref", "refsite", seed = 71)
  lea      <- .pi_samples(.pi_monthly("2017-07-01", "2022-05-01"), 1, "lea", "leasite", seed = 72)

  ## A single strongly atypical reference sample on a given day.
  odd_ref <- function(day) {
    .pi_samples(day, 0, paste0("odd", format(as.Date(day), "%Y%m%d")),
                "refsite", scale = 3, seed = 73)
  }
  lmf_with <- function(extra) {
    .pi_lmf_for(.tr_lmf(q, lea, dplyr::bind_rows(ref_base, extra)), sid)$value
  }
  base <- .pi_lmf_for(.tr_lmf(q, lea, ref_base), sid)$value
  expect_false(is.na(base))

  ## Lower bound: exactly t - 5 years is IN; the day before is OUT.
  expect_false(isTRUE(all.equal(lmf_with(odd_ref("2017-06-01")), base)))
  expect_identical(lmf_with(odd_ref("2017-05-31")), base)

  ## Upper bound: the sample's own day is IN; the day after is OUT.
  expect_false(isTRUE(all.equal(lmf_with(odd_ref("2022-06-01")), base)))
  expect_identical(lmf_with(odd_ref("2022-06-02")), base)
})

test_that("calibration_window_years sets the trailing window length", {
  q   <- .tr_query("2022-06-01", f = 0.5, seed = 61)
  sid <- unique(q$sample_id)
  ref <- .pi_samples(.pi_monthly("2019-07-01", "2022-05-01"), 0, "ref", "refsite", seed = 81)
  lea <- .pi_samples(.pi_monthly("2019-07-01", "2022-05-01"), 1, "lea", "leasite", seed = 82)
  ## Atypical reference sample 4 years before the sample: inside a 5-year
  ## window, outside a 3-year one.
  odd <- .pi_samples("2018-06-01", 0, "odd", "refsite", scale = 3, seed = 83)

  v3_base <- .pi_lmf_for(.tr_lmf(q, lea, ref, calibration_window_years = 3), sid)$value
  v3_odd  <- .pi_lmf_for(.tr_lmf(q, lea, dplyr::bind_rows(ref, odd),
                                 calibration_window_years = 3), sid)$value
  v5_base <- .pi_lmf_for(.tr_lmf(q, lea, ref, calibration_window_years = 5), sid)$value
  v5_odd  <- .pi_lmf_for(.tr_lmf(q, lea, dplyr::bind_rows(ref, odd),
                                 calibration_window_years = 5), sid)$value

  expect_identical(v3_odd, v3_base)                          # excluded at 3 y
  expect_false(isTRUE(all.equal(v5_odd, v5_base)))           # included at 5 y
})


## ---------------------------------------------------------------------------
## (b) No usable calibration data in the window -> NA with a reason
## ---------------------------------------------------------------------------

test_that("a sample dated before all end-member data gets NA with insufficient_calibration_data", {
  ref <- .pi_samples(.pi_monthly("2015-01-01", "2020-12-01"), 0, "ref", "refsite", seed = 91)
  lea <- .pi_samples(.pi_monthly("2015-01-01", "2020-12-01"), 1, "lea", "leasite", seed = 92)
  q_early <- .tr_query("2014-06-01", prefix = "qe", seed = 93)
  q_ok    <- .tr_query("2020-01-01", prefix = "qo", seed = 94)

  out <- expect_no_error(.tr_lmf(dplyr::bind_rows(q_early, q_ok), lea, ref))
  bad  <- .pi_lmf_for(out, unique(q_early$sample_id))
  good <- .pi_lmf_for(out, unique(q_ok$sample_id))

  expect_true(is.na(bad$value))
  expect_true(is.na(bad$sigma_lmf))
  expect_match(bad$lmf_reason, "^insufficient_calibration_data")
  expect_false(is.na(good$value))
  expect_true(is.na(good$lmf_reason))

  ## Input rows are preserved and each sample gets exactly one LMF row.
  expect_equal(sum(out$analyte != "LMF"), nrow(q_early) + nrow(q_ok))
  expect_equal(sum(out$analyte == "LMF"), 2L)
})

test_that("reference data in the window but no leachate data gives insufficient_calibration_data", {
  ref <- .pi_samples(.pi_monthly("2015-01-01", "2023-12-01"), 0, "ref", "refsite", seed = 101)
  ## Leachate sampling only starts in 2021.
  lea <- .pi_samples(.pi_monthly("2021-01-01", "2023-12-01"), 1, "lea", "leasite", seed = 102)
  q_before <- .tr_query("2020-06-01", prefix = "qb", seed = 103)
  q_after  <- .tr_query("2023-06-01", prefix = "qa", seed = 104)

  out <- .tr_lmf(dplyr::bind_rows(q_before, q_after), lea, ref)
  bad  <- .pi_lmf_for(out, unique(q_before$sample_id))
  good <- .pi_lmf_for(out, unique(q_after$sample_id))
  expect_true(is.na(bad$value))
  expect_match(bad$lmf_reason, "^insufficient_calibration_data")
  expect_false(is.na(good$value))
})

test_that("leachate data in the window but no reference data gives insufficient_calibration_data", {
  ## Reference sampling only starts in 2021.
  ref <- .pi_samples(.pi_monthly("2021-01-01", "2023-12-01"), 0, "ref", "refsite", seed = 111)
  lea <- .pi_samples(.pi_monthly("2015-01-01", "2023-12-01"), 1, "lea", "leasite", seed = 112)
  q_before <- .tr_query("2020-06-01", prefix = "qb", seed = 113)

  out <- .tr_lmf(q_before, lea, ref)
  bad <- .pi_lmf_for(out, unique(q_before$sample_id))
  expect_true(is.na(bad$value))
  expect_match(bad$lmf_reason, "^insufficient_calibration_data")
})

test_that("too few reference samples in the window to estimate spread gives insufficient_calibration_data", {
  ## Only two reference samples fall in the 2020-06-01 window: no ion has the
  ## >= 3 observations a reference SD needs.
  ref <- dplyr::bind_rows(
    .pi_samples(c("2019-01-01", "2019-06-01"), 0, "refA", "refsite", seed = 121),
    .pi_samples(.pi_monthly("2021-01-01", "2023-12-01"), 0, "refB", "refsite", seed = 122)
  )
  lea <- .pi_samples(.pi_monthly("2015-01-01", "2023-12-01"), 1, "lea", "leasite", seed = 123)
  q   <- .tr_query("2020-06-01", seed = 124)

  out <- expect_no_error(.tr_lmf(q, lea, ref))
  bad <- .pi_lmf_for(out, unique(q$sample_id))
  expect_true(is.na(bad$value))
  expect_match(bad$lmf_reason, "^insufficient_calibration_data")
})

test_that("every sample lacking calibration data still returns one NA row each, no error", {
  ref <- .pi_samples(.pi_monthly("2021-01-01", "2023-12-01"), 0, "ref", "refsite", seed = 131)
  lea <- .pi_samples(.pi_monthly("2021-01-01", "2023-12-01"), 1, "lea", "leasite", seed = 132)
  q   <- dplyr::bind_rows(.tr_query("2010-01-01", prefix = "qa", seed = 133),
                          .tr_query("2012-01-01", prefix = "qb", seed = 134))
  out <- expect_no_error(.tr_lmf(q, lea, ref))
  lmf <- .pi_lmf(out)
  expect_equal(nrow(lmf), 2L)
  expect_true(all(is.na(lmf$value)))
  expect_true(all(grepl("^insufficient_calibration_data", lmf$lmf_reason)))
})


## ---------------------------------------------------------------------------
## (c) Trailing calibration needs explicit, dated end-member data
## ---------------------------------------------------------------------------

test_that("calibration = 'trailing' errors without leachate_data or reference_data", {
  ref <- .pi_samples(.pi_monthly("2018-01-01", "2022-12-01"), 0, "ref", "refsite", seed = 141)
  lea <- .pi_samples(.pi_monthly("2018-01-01", "2022-12-01"), 1, "lea", "leasite", seed = 142)
  q   <- .tr_query("2022-06-01", seed = 143)
  .pi_stub_dashboard_globals()

  expect_error(
    add_lmf(q, reference_data = ref, calibration = "trailing"),
    "leachate_data"
  )
  expect_error(
    add_lmf(q, leachate_data = lea, calibration = "trailing"),
    "reference_data"
  )
})

test_that("calibration = 'trailing' errors when end-member data have no datetime", {
  ref <- .pi_samples(.pi_monthly("2018-01-01", "2022-12-01"), 0, "ref", "refsite", seed = 151)
  lea <- .pi_samples(.pi_monthly("2018-01-01", "2022-12-01"), 1, "lea", "leasite", seed = 152)
  q   <- .tr_query("2022-06-01", seed = 153)

  expect_error(
    add_lmf(q, leachate_data = lea, reference_data = ref[, setdiff(names(ref), "datetime")],
            calibration = "trailing"),
    "datetime"
  )
  expect_error(
    add_lmf(q, leachate_data = lea[, setdiff(names(lea), "datetime")], reference_data = ref,
            calibration = "trailing"),
    "datetime"
  )
})
