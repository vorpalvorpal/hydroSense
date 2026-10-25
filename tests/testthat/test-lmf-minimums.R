## ============================================================================
## add_lmf(): end-member minimums are enforced in both calibration modes
## ============================================================================
##
## Ruling (coordinator, after the dashboard-global path was removed):
## min_ref_samples and min_leachate_total_n are enforced in BOTH modes.
##   * "trailing": per sample window. A window below either minimum gives that
##     sample an NA LMF with lmf_reason starting "insufficient_calibration_data".
##   * "window": for the whole call. Every sample gets that NA + reason; no
##     error.
## Defaults unchanged: min_ref_samples = 10, min_leachate_total_n = NULL
## (-> 20 mg/L), min_leachate_samples = 10.
##
## Counting, matched to the removed builders (origin/main R/lmf.R):
##   * Reference (build_reference_endmember): n = distinct sample_id carrying
##     at least one LMF panel ion (nrow of the collapsed wide table); usable iff
##     n >= min_ref_samples.
##   * Leachate (build_leachate_endmember): a sample is valid iff its total N,
##     sum(NH3-N, NO3-N, NO2-N) in mg/L, is >= min_leachate_total_n; usable iff
##     n_distinct(valid sample_id) >= min_leachate_samples. Invalid samples do
##     not contribute to the end-member at all.
## The old code summed raw values in implicit mg/L; here total N must respect
## the row's units (the bundled demo data carry N species in ug/L).
##
## All end-member data are dated 2021, before the 2022-06-01 query samples, so
## the same fixtures sit inside the trailing 5-year window and the specs run
## identically in both modes.

library(testthat)
library(hydroSense)

## n monthly dates from 2021-01-01.
.mn_dates <- function(n, from = "2021-01-01") {
  seq(as.Date(from), by = "month", length.out = n)
}

.mn_ref <- function(n, seed = 201, prefix = "ref", from = "2021-01-01") {
  .pi_samples(.mn_dates(n, from), 0, prefix, "refsite", seed = seed)
}

.mn_lea <- function(n, seed = 202, prefix = "lea", from = "2021-01-01") {
  .pi_samples(.mn_dates(n, from), 1, prefix, "leasite", seed = seed)
}

## One atypical leachate sample with chosen N species (exact in floating
## point: e.g. 15 + 4 + 1 == 20) and Na/Ca tripled, so its inclusion visibly
## moves the end-member.
.mn_odd_lea <- function(nh3, no3 = 4, no2 = 1, units = "mg/L",
                        day = "2021-12-15", prefix = "odd") {
  s <- .pi_samples(day, 1, prefix, "leasite", seed = 203)
  s$value[s$analyte == "NH3-N"] <- nh3
  s$value[s$analyte == "NO3-N"] <- no3
  s$value[s$analyte == "NO2-N"] <- no2
  s$units.analyte[s$analyte %in% c("NH3-N", "NO3-N", "NO2-N")] <- units
  s$value[s$analyte %in% c("Na", "Ca")] <- 3 * s$value[s$analyte %in% c("Na", "Ca")]
  s
}

## Two query samples on 2022-06-01.
.mn_query <- function() {
  dplyr::bind_rows(
    .pi_samples("2022-06-01", 0.3, "qa", "query", seed = 204),
    .pi_samples("2022-06-01", 0.7, "qb", "query", seed = 205)
  )
}

.mn_run <- function(mode, lea, ref, q = .mn_query(), ...) {
  add_lmf(q, leachate_data = lea, reference_data = ref, calibration = mode,
          max_sigma_lsi = 1e6, ...)
}

.mn_expect_all_insufficient <- function(out, q = .mn_query()) {
  lmf <- .pi_lmf(out)
  expect_equal(nrow(lmf), dplyr::n_distinct(q$sample_id))
  expect_true(all(is.na(lmf$value)))
  expect_true(all(grepl("^insufficient_calibration_data", lmf$lmf_reason)))
  expect_equal(sum(out$analyte != "LMF"), nrow(q))
}

.mn_expect_all_quantified <- function(out, q = .mn_query()) {
  lmf <- .pi_lmf(out)
  expect_equal(nrow(lmf), dplyr::n_distinct(q$sample_id))
  expect_false(anyNA(lmf$value))
  expect_true(all(is.na(lmf$lmf_reason)))
}

.mn_values <- function(out) .pi_lmf(out)$value


for (mode in c("window", "trailing")) {

  ## -------------------------------------------------------------------------
  ## min_ref_samples
  ## -------------------------------------------------------------------------

  test_that(paste0("[", mode, "] exactly min_ref_samples (default 10) reference samples is usable"), {
    out <- .mn_run(mode, .mn_lea(12), .mn_ref(10))
    .mn_expect_all_quantified(out)
  })

  test_that(paste0("[", mode, "] one reference sample below the default minimum gives NA, no error"), {
    out <- expect_no_error(.mn_run(mode, .mn_lea(12), .mn_ref(9)))
    .mn_expect_all_insufficient(out)
  })

  test_that(paste0("[", mode, "] a non-default min_ref_samples is honoured at its boundary"), {
    .mn_expect_all_quantified(.mn_run(mode, .mn_lea(12), .mn_ref(5), min_ref_samples = 5))
    out <- expect_no_error(.mn_run(mode, .mn_lea(12), .mn_ref(4), min_ref_samples = 5))
    .mn_expect_all_insufficient(out)
  })

  test_that(paste0("[", mode, "] reference count is distinct samples carrying panel ions, not rows"), {
    ## 9 panel samples (108 rows, far more than 10) plus one sample with only
    ## non-panel analytes: 9 countable samples, below the minimum of 10.
    ref9 <- .mn_ref(9)
    extra <- tibble::tibble(
      sample_id = "ref_nonpanel", site_id = "refsite",
      analyte = c("pH", "EC"), value = c(7.2, 250), detected = TRUE,
      datetime = .pi_utc("2021-11-15"), units.analyte = c("pH", "uS/cm"),
      valence.analyte = NA_real_, atomic_mass.analyte = NA_real_
    )
    out <- expect_no_error(.mn_run(mode, .mn_lea(12), dplyr::bind_rows(ref9, extra)))
    .mn_expect_all_insufficient(out)
  })

  ## -------------------------------------------------------------------------
  ## min_leachate_total_n (per-sample validity filter)
  ## -------------------------------------------------------------------------

  test_that(paste0("[", mode, "] default 20 mg/L total-N filter: exactly 20 is valid, 19.99 is not"), {
    ref  <- .mn_ref(12)
    lea  <- .mn_lea(12)
    base <- .mn_values(.mn_run(mode, lea, ref))
    expect_false(anyNA(base))

    at    <- .mn_values(.mn_run(mode, dplyr::bind_rows(lea, .mn_odd_lea(nh3 = 15)), ref))
    below <- .mn_values(.mn_run(mode, dplyr::bind_rows(lea, .mn_odd_lea(nh3 = 14.99)), ref))
    expect_false(isTRUE(all.equal(at, base)))   # total N 20.00: used
    expect_identical(below, base)               # total N 19.99: ignored
  })

  test_that(paste0("[", mode, "] an explicit min_leachate_total_n is honoured at its boundary"), {
    ref  <- .mn_ref(12)
    lea  <- .mn_lea(12)
    run  <- function(l) {
      .mn_values(.mn_run(mode, l, ref, min_leachate_total_n = 50,
                         min_leachate_total_n_units = "mg/L"))
    }
    base  <- run(lea)
    at    <- run(dplyr::bind_rows(lea, .mn_odd_lea(nh3 = 45)))      # 45 + 4 + 1 = 50
    below <- run(dplyr::bind_rows(lea, .mn_odd_lea(nh3 = 44.99)))   # 49.99
    expect_false(isTRUE(all.equal(at, base)))
    expect_identical(below, base)
  })

  test_that(paste0("[", mode, "] leachate total N respects the N species' units"), {
    ref  <- .mn_ref(12)
    lea  <- .mn_lea(12)
    base <- .mn_values(.mn_run(mode, lea, ref))
    ## 15000 + 4000 + 1000 ug/L = 20 mg/L: valid. 14990 ug/L NH3-N: 19.99 mg/L.
    at    <- .mn_values(.mn_run(mode, dplyr::bind_rows(
      lea, .mn_odd_lea(nh3 = 15000, no3 = 4000, no2 = 1000, units = "ug/L")), ref))
    below <- .mn_values(.mn_run(mode, dplyr::bind_rows(
      lea, .mn_odd_lea(nh3 = 14990, no3 = 4000, no2 = 1000, units = "ug/L")), ref))
    expect_false(isTRUE(all.equal(at, base)))
    expect_identical(below, base)
  })

  test_that(paste0("[", mode, "] no leachate sample passing the total-N filter gives NA, no error"), {
    ## Every leachate sample has total N of 10 mg/L (< 20).
    lea <- dplyr::bind_rows(lapply(seq_len(12), function(i) {
      .mn_odd_lea(nh3 = 5, day = as.character(.mn_dates(12)[i]), prefix = paste0("lowN", i))
    }))
    out <- expect_no_error(.mn_run(mode, lea, .mn_ref(12)))
    .mn_expect_all_insufficient(out)
  })

  ## -------------------------------------------------------------------------
  ## min_leachate_samples counts VALID (total-N passing) samples
  ## (the removed build_leachate_endmember() gate; see header)
  ## -------------------------------------------------------------------------

  test_that(paste0("[", mode, "] exactly min_leachate_samples valid leachate samples is usable; one fewer is not"), {
    ref <- .mn_ref(12)
    .mn_expect_all_quantified(.mn_run(mode, .mn_lea(10), ref))

    ## 12 leachate samples, 3 of them below the total-N filter: 9 valid.
    low <- dplyr::bind_rows(lapply(1:3, function(i) {
      .mn_odd_lea(nh3 = 5, day = as.character(.mn_dates(3, "2021-02-10")[i]),
                  prefix = paste0("lowN", i))
    }))
    out <- expect_no_error(.mn_run(mode, dplyr::bind_rows(.mn_lea(9), low), ref))
    .mn_expect_all_insufficient(out)
  })
}


## ---------------------------------------------------------------------------
## Trailing: minimums apply to each sample's own window
## ---------------------------------------------------------------------------

test_that("[trailing] min_ref_samples counts only reference samples inside the window", {
  lea     <- .mn_lea(12)
  old_ref <- .mn_ref(30, prefix = "refold", from = "2010-01-01")   # outside the window
  q       <- .mn_query()

  out_ok <- .mn_run("trailing", lea, dplyr::bind_rows(old_ref, .mn_ref(10)), q)
  .mn_expect_all_quantified(out_ok, q)

  out_short <- expect_no_error(
    .mn_run("trailing", lea, dplyr::bind_rows(old_ref, .mn_ref(9)), q)
  )
  .mn_expect_all_insufficient(out_short, q)
})

test_that("[trailing] the minimum is judged per sample: an early sample fails, a later one passes", {
  ## Reference samples monthly through 2021 (12). A sample on 2021-06-15 sees
  ## only 6 of them; a sample on 2022-06-01 sees all 12.
  lea <- .mn_lea(12, from = "2020-01-01")
  ref <- .mn_ref(12)
  q   <- dplyr::bind_rows(
    .pi_samples("2021-06-15", 0.5, "qearly", "query", seed = 206),
    .pi_samples("2022-06-01", 0.5, "qlate",  "query", seed = 207)
  )
  out <- expect_no_error(.mn_run("trailing", lea, ref, q))
  early <- .pi_lmf_for(out, unique(q$sample_id[grepl("^qearly", q$sample_id)]))
  late  <- .pi_lmf_for(out, unique(q$sample_id[grepl("^qlate", q$sample_id)]))
  expect_true(is.na(early$value))
  expect_match(early$lmf_reason, "^insufficient_calibration_data")
  expect_false(is.na(late$value))
})

test_that("[trailing] total-N-invalid leachate samples in the window are not counted", {
  ## The 2022-06-01 window holds 9 valid + 3 low-N leachate samples (9 valid
  ## < 10); 20 further valid samples from 2010-2011 lie outside it and must
  ## not make up the shortfall.
  ref <- .mn_ref(12)
  low <- dplyr::bind_rows(lapply(1:3, function(i) {
    .mn_odd_lea(nh3 = 5, day = as.character(.mn_dates(3, "2021-02-10")[i]),
                prefix = paste0("lowN", i))
  }))
  old_valid <- .mn_lea(20, prefix = "leaold", from = "2010-01-01")   # outside the window
  out <- expect_no_error(
    .mn_run("trailing", dplyr::bind_rows(old_valid, .mn_lea(9), low), ref)
  )
  .mn_expect_all_insufficient(out)
})
