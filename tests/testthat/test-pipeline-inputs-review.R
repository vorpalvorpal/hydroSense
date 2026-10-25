## ============================================================================
## Specs from the review of a14c865 (coordinator rulings 1, 2 and 4)
## ============================================================================
##
## 1. Combined species. A combined analyte (NO2+NO3-N; Alkalinity-total-CaCO3)
##    is dropped only when ALL the split species it stands for are present in
##    that sample. If only some are present, the combined value is used and
##    the partial split species are dropped. Applies to the LMF ions and to the
##    leachate total-N filter.
## 2. add_mspaf() returns the caller's rows with their ORIGINAL values and
##    units; unit conversion to ug/L is internal only. msPAF values unchanged.
## 4. Trailing window bounds are compared at whole-calendar-day resolution:
##    an end-member on the same calendar day as the target (as.Date in the
##    datetime's own time zone, UTC when none) is inside the window whatever
##    its clock time; the lower bound is inclusive by day.

library(testthat)
library(hydroSense)


## ---------------------------------------------------------------------------
## 1. Combined species with a partial split
## ---------------------------------------------------------------------------

.cs_em <- function() {
  set.seed(101)
  list(
    ref = .pi_samples(rep("2022-06-01", 12), 0, "ref", "refsite"),
    lea = .pi_samples(rep("2022-06-01", 12), 1, "lea", "leasite")
  )
}

.cs_query <- function() {
  fvec <- c(0.2, 0.5, 0.8)
  set.seed(7)
  dplyr::bind_rows(lapply(seq_along(fvec), function(i) {
    .pi_samples("2022-06-01", fvec[i], paste0("q", i), "query")
  }))
}

## Per sample: replace `drop` rows by one `combined` row whose value is the
## sum of the `sum_of` species; rows in neither set are kept (so a partial
## split species can be left in place).
.cs_partial <- function(df, sum_of, drop, combined) {
  sums <- df[df$analyte %in% sum_of, ] |>
    dplyr::group_by(sample_id) |>
    dplyr::summarise(value = sum(value), .groups = "drop")
  tmpl <- df[df$analyte == sum_of[1], ]
  tmpl$value   <- sums$value[match(tmpl$sample_id, sums$sample_id)]
  tmpl$analyte <- combined
  dplyr::bind_rows(df[!df$analyte %in% drop, ], tmpl)
}

.cs_lmf <- function(q, em) {
  out <- .pi_lmf(add_lmf(q, leachate_data = em$lea, reference_data = em$ref,
                         max_sigma_lsi = 1e6))
  testthat::expect_false(anyNA(out$value))
  out
}

test_that("NO2+NO3-N with only NO2-N present: the combined value is used, NO2-N dropped", {
  em   <- .cs_em()
  q    <- .cs_query()
  full <- .cs_lmf(q, em)

  ## Samples carry NH3-N, NO2-N and NO2+NO3-N (= NO3-N + NO2-N); no NO3-N.
  ## Total N must be NH3-N + (NO2+NO3-N), i.e. the full-split total.
  part <- function(d) .cs_partial(d, c("NO3-N", "NO2-N"), "NO3-N", "NO2+NO3-N")
  q_p  <- part(q)
  expect_true(all(c("NO2-N", "NO2+NO3-N") %in% q_p$analyte))
  expect_false("NO3-N" %in% q_p$analyte)

  got <- .cs_lmf(q_p, list(ref = part(em$ref), lea = part(em$lea)))
  expect_equal(got$value, full$value, tolerance = 1e-6)
  expect_equal(got$sigma_lmf, full$sigma_lmf, tolerance = 1e-6)
})

test_that("NO2+NO3-N with only NO3-N present: the combined value is used, NO3-N dropped", {
  em   <- .cs_em()
  q    <- .cs_query()
  full <- .cs_lmf(q, em)
  part <- function(d) .cs_partial(d, c("NO3-N", "NO2-N"), "NO2-N", "NO2+NO3-N")
  got  <- .cs_lmf(part(q), list(ref = part(em$ref), lea = part(em$lea)))
  expect_equal(got$value, full$value, tolerance = 1e-6)
})

test_that("worked example: NH3-N 1, NO2-N 0.01, NO2+NO3-N 5 mg/L books total N as 1 + 5 mg/L", {
  ## Same LMF as a sample reporting NH3-N 1 and NO3-N 5 (no NO2-N): total N
  ## (1 + 5) / 14.0067 ~= 0.428 meq/L either way.
  em <- .cs_em()
  base <- .pi_samples("2022-06-01", 0.5, "qx", "query", seed = 9)
  base <- base[!base$analyte %in% c("NH3-N", "NO3-N", "NO2-N"), ]
  n_row <- function(analyte, value) {
    r <- base[base$analyte == "Na", ][1, ]
    r$analyte <- analyte; r$value <- value
    r$valence.analyte <- 1; r$atomic_mass.analyte <- 14.0067
    r
  }
  q_comb <- dplyr::bind_rows(base, n_row("NH3-N", 1), n_row("NO2-N", 0.01),
                             n_row("NO2+NO3-N", 5))
  q_ref  <- dplyr::bind_rows(base, n_row("NH3-N", 1), n_row("NO3-N", 5))
  expect_equal(.cs_lmf(q_comb, em)$value, .cs_lmf(q_ref, em)$value, tolerance = 1e-6)
})

test_that("Alkalinity-total-CaCO3 with only CO3-CaCO3 present: total alkalinity is used", {
  em   <- .cs_em()
  q    <- .cs_query()
  full <- .cs_lmf(q, em)
  part <- function(d) .cs_partial(d, c("CO3-CaCO3", "HCO3-CaCO3"), "HCO3-CaCO3",
                                  "Alkalinity-total-CaCO3")
  q_p  <- part(q)
  expect_true(all(c("CO3-CaCO3", "Alkalinity-total-CaCO3") %in% q_p$analyte))
  got <- .cs_lmf(q_p, list(ref = part(em$ref), lea = part(em$lea)))
  expect_equal(got$value, full$value, tolerance = 1e-6)
})

test_that("worked example: CO3-CaCO3 <1 (non-detect) with total alkalinity 200 mg/L gives alkalinity 200", {
  em   <- .cs_em()
  base <- .pi_samples("2022-06-01", 0.5, "qy", "query", seed = 10)
  base <- base[!base$analyte %in% c("CO3-CaCO3", "HCO3-CaCO3"), ]
  alk_row <- function(analyte, value, detected = TRUE) {
    r <- base[base$analyte == "Ca", ][1, ]
    r$analyte <- analyte; r$value <- value; r$detected <- detected
    r$valence.analyte <- 2; r$atomic_mass.analyte <- 100.0869
    r
  }
  q_comb <- dplyr::bind_rows(base, alk_row("CO3-CaCO3", 1, detected = FALSE),
                             alk_row("Alkalinity-total-CaCO3", 200))
  q_ref  <- dplyr::bind_rows(base, alk_row("HCO3-CaCO3", 200))
  a <- .cs_lmf(q_comb, em)
  b <- .cs_lmf(q_ref, em)
  expect_equal(a$value, b$value, tolerance = 1e-6)
  expect_equal(a$sigma_lmf, b$sigma_lmf, tolerance = 1e-6)
})

test_that("leachate total-N filter uses the combined value when the split is partial", {
  ## Default filter: total N >= 20 mg/L. An atypical leachate sample reports
  ## NH3-N 15, NO2-N 1 and NO2+NO3-N 5 (no NO3-N): total N = 15 + 5 = 20, so
  ## it is valid and moves the end-member. With NO2+NO3-N 4.99 it is 19.99 and
  ## ignored.
  ref <- .pi_samples(seq(as.Date("2021-01-01"), by = "month", length.out = 12),
                     0, "ref", "refsite", seed = 201)
  lea <- .pi_samples(seq(as.Date("2021-01-01"), by = "month", length.out = 12),
                     1, "lea", "leasite", seed = 202)
  q   <- .pi_samples("2022-06-01", 0.5, "q", "query", seed = 204)

  odd <- function(nox) {
    s <- .pi_samples("2021-12-15", 1, "odd", "leasite", seed = 203)
    s <- s[s$analyte != "NO3-N", ]
    s$value[s$analyte == "NH3-N"] <- 15
    s$value[s$analyte == "NO2-N"] <- 1
    comb <- s[s$analyte == "NO2-N", ]
    comb$analyte <- "NO2+NO3-N"; comb$value <- nox
    s$value[s$analyte %in% c("Na", "Ca")] <- 3 * s$value[s$analyte %in% c("Na", "Ca")]
    dplyr::bind_rows(s, comb)
  }
  for (mode in c("window", "trailing")) {
    run <- function(l) {
      .pi_lmf(add_lmf(q, leachate_data = l, reference_data = ref,
                      calibration = mode, max_sigma_lsi = 1e6))$value
    }
    base  <- run(lea)
    at    <- run(dplyr::bind_rows(lea, odd(5)))
    below <- run(dplyr::bind_rows(lea, odd(4.99)))
    expect_false(isTRUE(all.equal(at, base)), info = mode)   # 20.00: valid
    expect_identical(below, base, info = mode)                # 19.99: ignored
  }
})


## ---------------------------------------------------------------------------
## 2. add_mspaf returns the caller's original values and units
## ---------------------------------------------------------------------------

## Toxicants converted from ug/L to mg/L in the caller's own unit column.
.to_mg <- function(df, unit_col) {
  tox <- df$analyte %in% c("Cu", "Zn", "Ni", "NH3-N")
  df$value[tox] <- df$value[tox] / 1000
  df[[unit_col]][tox] <- "mg/L"
  df
}

.input_rows_match <- function(out, inp) {
  got <- as.data.frame(out[out$analyte != "msPAF", names(inp)])
  exp <- as.data.frame(inp)
  got <- got[order(got$sample_id, got$analyte), ]
  exp <- exp[order(exp$sample_id, exp$analyte), ]
  expect_equal(got, exp, ignore_attr = TRUE)
}

test_that("add_mspaf returns input rows with the caller's values and units.analyte", {
  tg_ug <- .pi_mspaf_chem("ds1", n = 3L, seed = 2L, mult = 5)
  tg_mg <- .to_mg(tg_ug, "units.analyte")

  out_mg <- suppressMessages(suppressWarnings(add_mspaf(tg_mg)))
  out_ug <- suppressMessages(suppressWarnings(add_mspaf(tg_ug)))

  .input_rows_match(out_mg, tg_mg)
  cu_in  <- tg_mg$value[tg_mg$analyte == "Cu"]
  cu_out <- out_mg$value[out_mg$analyte == "Cu"]
  expect_equal(sort(cu_out), sort(cu_in))
  expect_true(all(out_mg$units.analyte[out_mg$analyte == "Cu"] == "mg/L"))

  ## msPAF itself is unit-independent.
  m_mg <- out_mg[out_mg$analyte == "msPAF", ]
  m_ug <- out_ug[out_ug$analyte == "msPAF", ]
  expect_equal(m_mg$value[order(m_mg$sample_id)], m_ug$value[order(m_ug$sample_id)],
               tolerance = 1e-8)
})

test_that("add_mspaf returns sampleTidy-shaped input rows unchanged (values and `units`)", {
  tg <- .pi_mspaf_sampletidy_names(.pi_mspaf_chem("ds1", n = 3L, seed = 2L, mult = 5))
  tg <- tibble::tibble(
    sample_id = tg$sample_id, site_id = tg$site_id, analyte = tg$analyte,
    value = tg$value, detected = tg$detected, datetime = tg$datetime,
    units = ifelse(is.na(tg$units.analyte), "mg/L", tg$units.analyte)
  )
  tg_mg <- .to_mg(tg, "units")

  out <- suppressMessages(suppressWarnings(add_mspaf(tg_mg)))
  .input_rows_match(out, tg_mg)
  expect_true(any(out$analyte == "msPAF"))
})

test_that("add_mspaf with a reference also leaves the caller's rows untouched", {
  tg  <- .to_mg(.pi_mspaf_chem("ds1", n = 3L, seed = 2L, mult = 5), "units.analyte")
  ref <- .to_mg(.pi_mspaf_chem("ref", n = 6L, seed = 1L), "units.analyte")
  out <- suppressMessages(suppressWarnings(add_mspaf(tg, reference = ref)))
  .input_rows_match(out, tg)
})


## ---------------------------------------------------------------------------
## 4. Trailing window at whole-calendar-day resolution
## ---------------------------------------------------------------------------

## Base end-members strictly inside the window of a 2022-06-01 target, never
## on a boundary day.
.day_em <- function() {
  list(
    ref = .pi_samples(.pi_monthly("2017-07-01", "2022-05-01"), 0, "ref", "refsite", seed = 71),
    lea = .pi_samples(.pi_monthly("2017-07-01", "2022-05-01"), 1, "lea", "leasite", seed = 72)
  )
}

## One strongly atypical reference sample stamped at `when` (POSIXct).
.odd_ref_at <- function(when) {
  s <- .pi_samples("2020-01-01", 0, paste0("odd", as.numeric(when)), "refsite",
                   scale = 3, seed = 73)
  s$datetime <- when
  s
}

.day_lmf <- function(q, em, extra = NULL) {
  ref <- if (is.null(extra)) em$ref else dplyr::bind_rows(em$ref, extra)
  out <- add_lmf(q, leachate_data = em$lea, reference_data = ref,
                 calibration = "trailing", max_sigma_lsi = 1e6)
  .pi_lmf(out)$value
}

.utc <- function(x) as.POSIXct(x, tz = "UTC")

test_that("Date target: a same-day 04:00 UTC end-member is inside the window", {
  em <- .day_em()
  q  <- .pi_samples("2022-06-01", 0.5, "q", "query", seed = 61)
  q$datetime <- as.Date("2022-06-01")
  base <- .day_lmf(q, em)
  expect_false(anyNA(base))

  same_day <- .day_lmf(q, em, .odd_ref_at(.utc("2022-06-01 04:00:00")))
  next_day <- .day_lmf(q, em, .odd_ref_at(.utc("2022-06-02 00:00:00")))
  expect_false(isTRUE(all.equal(same_day, base)))   # included
  expect_identical(next_day, base)                  # excluded
})

test_that("Date target: the lower bound is inclusive by calendar day", {
  em <- .day_em()
  q  <- .pi_samples("2022-06-01", 0.5, "q", "query", seed = 61)
  q$datetime <- as.Date("2022-06-01")
  base <- .day_lmf(q, em)

  on_lower  <- .day_lmf(q, em, .odd_ref_at(.utc("2017-06-01 23:00:00")))
  day_out   <- .day_lmf(q, em, .odd_ref_at(.utc("2017-05-31 23:59:00")))
  expect_false(isTRUE(all.equal(on_lower, base)))   # 2017-06-01: included
  expect_identical(day_out, base)                   # 2017-05-31: excluded
})

test_that("POSIXct target with a clock time: bounds are whole days at both ends", {
  em <- .day_em()
  q  <- .pi_samples("2022-06-01", 0.5, "q", "query", seed = 61)
  q$datetime <- .utc("2022-06-01 09:00:00")
  base <- .day_lmf(q, em)
  expect_false(anyNA(base))

  ## Later the same day than the target: still inside.
  expect_false(isTRUE(all.equal(.day_lmf(q, em, .odd_ref_at(.utc("2022-06-01 20:00:00"))), base)))
  ## Earlier in the day than the target's clock time, on the lower-bound day.
  expect_false(isTRUE(all.equal(.day_lmf(q, em, .odd_ref_at(.utc("2017-06-01 03:00:00"))), base)))
  ## Day after / day before: outside.
  expect_identical(.day_lmf(q, em, .odd_ref_at(.utc("2022-06-02 00:30:00"))), base)
  expect_identical(.day_lmf(q, em, .odd_ref_at(.utc("2017-05-31 20:00:00"))), base)
})

test_that("the calendar day of an end-member is taken in its own time zone", {
  ## End-members in Australia/Sydney (UTC+10 in June). 2022-06-02 08:00 AEST is
  ## 2022-06-01 22:00 UTC: its own calendar day is 2022-06-02, after a
  ## 2022-06-01 target, so it is OUTSIDE. 2022-06-01 08:00 AEST
  ## (2022-05-31 22:00 UTC) is on the target's day, so it is INSIDE.
  syd <- function(df, hhmm = "10:00") {
    df$datetime <- as.POSIXct(paste(format(df$datetime, "%Y-%m-%d"), hhmm),
                              tz = "Australia/Sydney")
    df
  }
  em <- .day_em()
  em <- list(ref = syd(em$ref), lea = syd(em$lea))
  q  <- .pi_samples("2022-06-01", 0.5, "q", "query", seed = 61)
  q$datetime <- as.Date("2022-06-01")
  base <- .day_lmf(q, em)

  odd_syd <- function(x) .odd_ref_at(as.POSIXct(x, tz = "Australia/Sydney"))
  inside  <- .day_lmf(q, em, odd_syd("2022-06-01 08:00:00"))
  outside <- .day_lmf(q, em, odd_syd("2022-06-02 08:00:00"))
  expect_false(isTRUE(all.equal(inside, base)))
  expect_identical(outside, base)
})
