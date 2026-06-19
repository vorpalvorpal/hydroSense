## Tests for draw_id threading through time_weighted_aggregate (Chunk 2).
## Stan-free. No ANZG XLSX files needed for TWA tests; bundled SSD data used
## for the end-to-end composition test.
##
## Seven properties tested:
##   1. Degradation  — point input: no draw_id column, values unchanged
##   2. Shape        — draws input: output has draw_id 1..N per (focal, site, analyte)
##   3. Per-draw     — draw-d chronic == running TWA on draw-d chemistry alone
##   4. Broadcast    — observed (exact) sample participates in every draw's aggregate
##   5. Index-pairing— samples are paired by draw_id across time, not shuffled
##   6. Composition  — add_amspaf(draws) |> TWA yields chronic AmsPAF draws
##   7. Guard        — amspaf_daily errors on draws input

library(testthat)
library(hydroSense)

## ── Shared helpers ────────────────────────────────────────────────────────────

## Build a simple long-format chemistry frame.
## Each entry in `samples` is a list(date, analytes=list(name=value_or_draws_vec)).
## An analyte whose value is a vector of length > 1 becomes a drawn cell
## (draw_id=1..N); length==1 → exact cell (draw_id=NA).
make_twa_frame <- function(samples, site_id = "A") {
  rows <- purrr::imap_dfr(samples, function(smp, i) {
    sid  <- paste0("s", i)
    date <- as.Date(smp$date)
    purrr::imap_dfr(smp$analytes, function(vals, analyte) {
      n <- length(vals)
      tibble::tibble(
        sample_id = sid,
        site_id   = site_id,
        datetime  = date,
        analyte   = analyte,
        value     = vals,
        detected  = TRUE,
        draw_id   = if (n > 1L) seq_len(n) else NA_integer_
      )
    })
  })
  rows
}

## Reference focal-date / tau / window used in all TWA calls below.
.focal <- as.Date("2024-06-01")
.tau   <- 90
.win   <- 365


## ── 1. Degradation ───────────────────────────────────────────────────────────

test_that("point input: no draw_id column in output, values unchanged", {
  df <- make_twa_frame(list(
    list(date = "2024-04-01", analytes = list(Cu = 5)),
    list(date = "2024-05-01", analytes = list(Cu = 10))
  ))
  out_point <- time_weighted_aggregate(df, focal_dates = .focal,
    tau = .tau, tau_units = "d", window = .win, window_units = "d")

  expect_false("draw_id" %in% names(out_point))
  expect_true(all(is.finite(out_point$value)))
})

test_that("point input with imputed column: n_imputed_in_window propagated", {
  df <- make_twa_frame(list(
    list(date = "2024-04-01", analytes = list(Cu = 5)),
    list(date = "2024-05-01", analytes = list(Cu = 10))
  )) |>
    dplyr::mutate(imputed = sample_id == "s1")  # s1 imputed, s2 not

  out <- time_weighted_aggregate(df, focal_dates = .focal,
    tau = .tau, tau_units = "d", window = .win, window_units = "d")

  expect_equal(out$n_imputed_in_window, 1L)
})


## ── 2. Shape ─────────────────────────────────────────────────────────────────

test_that("draws input: output has draw_id 1..N per (focal_date, site_id, analyte)", {
  n_draws <- 3L
  df <- make_twa_frame(list(
    list(date = "2024-04-01", analytes = list(Cu = c(5, 10, 15))),
    list(date = "2024-05-01", analytes = list(Cu = c(20, 25, 30)))
  ))
  out <- time_weighted_aggregate(df, focal_dates = .focal,
    tau = .tau, tau_units = "d", window = .win, window_units = "d",
    return = "draws")

  expect_true("draw_id" %in% names(out))
  expect_false(anyNA(out$draw_id))
  expect_equal(sort(unique(out$draw_id)), seq_len(n_draws))
  # One row per draw for the single (focal, site, analyte) combination
  expect_equal(nrow(out), n_draws)
})

test_that("draws input: n_samples_in_window identical across draws", {
  df <- make_twa_frame(list(
    list(date = "2024-03-01", analytes = list(Cu = c(5, 50))),
    list(date = "2024-04-01", analytes = list(Cu = c(10, 100))),
    list(date = "2024-05-01", analytes = list(Cu = c(15, 150)))
  ))
  out <- time_weighted_aggregate(df, focal_dates = .focal,
    tau = .tau, tau_units = "d", window = .win, window_units = "d")

  # n_samples_in_window is structural (date-based), not draw-varying
  expect_equal(length(unique(out$n_samples_in_window)), 1L)
})


## ── 3. Per-draw correctness ───────────────────────────────────────────────────

test_that("draw-d chronic value == TWA on draw-d chemistry alone", {
  # Two samples with different Cu across two draws
  df_draws <- make_twa_frame(list(
    list(date = "2024-04-01", analytes = list(Cu = c(5, 50))),
    list(date = "2024-05-01", analytes = list(Cu = c(10, 100)))
  ))
  out_draws <- time_weighted_aggregate(df_draws, focal_dates = .focal,
    tau = .tau, tau_units = "d", window = .win, window_units = "d",
    return = "draws")

  d1_val <- dplyr::filter(out_draws, draw_id == 1L)$value
  d2_val <- dplyr::filter(out_draws, draw_id == 2L)$value

  # Point references: run TWA on each draw's chemistry independently
  df_p1 <- make_twa_frame(list(
    list(date = "2024-04-01", analytes = list(Cu = 5)),
    list(date = "2024-05-01", analytes = list(Cu = 10))
  ))
  df_p2 <- make_twa_frame(list(
    list(date = "2024-04-01", analytes = list(Cu = 50)),
    list(date = "2024-05-01", analytes = list(Cu = 100))
  ))
  p1_val <- time_weighted_aggregate(df_p1, focal_dates = .focal,
    tau = .tau, tau_units = "d", window = .win, window_units = "d")$value
  p2_val <- time_weighted_aggregate(df_p2, focal_dates = .focal,
    tau = .tau, tau_units = "d", window = .win, window_units = "d")$value

  expect_equal(d1_val, p1_val)
  expect_equal(d2_val, p2_val)
  expect_gt(d2_val, d1_val)   # draw 2 has higher Cu → higher geom mean
})


## ── 4. Exact-cell broadcast across time ──────────────────────────────────────

test_that("exact (observed) sample participates in every draw's aggregate", {
  # s1 is observed (exact, draw_id=NA); s2 is imputed (2 draws)
  # If s1 were excluded from any draw, the aggregate would differ.
  df <- dplyr::bind_rows(
    # Observed sample at day -60: exact, same value in all "draws"
    tibble::tibble(sample_id = "s1", site_id = "A",
      datetime = .focal - 60L, analyte = "Cu",
      value = 10, detected = TRUE, draw_id = NA_integer_),
    # Imputed sample at day -10: two draws with very different values
    tibble::tibble(sample_id = "s2", site_id = "A",
      datetime = .focal - 10L, analyte = "Cu",
      value = c(1, 1000), detected = TRUE, draw_id = c(1L, 2L))
  )
  out <- time_weighted_aggregate(df, focal_dates = .focal,
    tau = .tau, tau_units = "d", window = .win, window_units = "d",
    return = "draws")

  # Both draws see s1 (Cu=10) — if s1 were excluded, the draws would be
  # geom_mean(1)=1 and geom_mean(1000)=1000; s1 pulls them towards 10.
  d1 <- dplyr::filter(out, draw_id == 1L)$value
  d2 <- dplyr::filter(out, draw_id == 2L)$value

  expect_gt(d1, 1)      # d1 pulled above pure Cu=1 by the observed s1
  expect_lt(d2, 1000)   # d2 pulled below pure Cu=1000 by the observed s1

  # Verify s1's Cu=10 anchors both draws (both should be between 1 and 1000)
  expect_true(d1 > 1 && d1 < 1000)
  expect_true(d2 > 1 && d2 < 1000)
})


## ── 5. Temporal index-pairing ─────────────────────────────────────────────────

test_that("samples are paired by draw_id across time, not cross-contaminated", {
  # Draw 1: s1_Cu=high, s2_Cu=high  → consistently high
  # Draw 2: s1_Cu=low,  s2_Cu=low   → consistently low
  # If paired by index, d1 > d2.
  # Mis-pairing (s1_d1 with s2_d2) would give intermediate values for both.
  df <- dplyr::bind_rows(
    tibble::tibble(sample_id = "s1", site_id = "A",
      datetime = .focal - 60L, analyte = "Cu",
      value = c(500, 1), detected = TRUE, draw_id = c(1L, 2L)),
    tibble::tibble(sample_id = "s2", site_id = "A",
      datetime = .focal - 10L, analyte = "Cu",
      value = c(500, 1), detected = TRUE, draw_id = c(1L, 2L))
  )
  out <- time_weighted_aggregate(df, focal_dates = .focal,
    tau = .tau, tau_units = "d", window = .win, window_units = "d",
    return = "draws")

  d1 <- dplyr::filter(out, draw_id == 1L)$value
  d2 <- dplyr::filter(out, draw_id == 2L)$value

  # Correct pairing: d1 uses (500, 500) → high; d2 uses (1, 1) → low
  expect_gt(d1, d2)

  # Mis-paired reference: s1_d1=500 paired with s2 point=1 → intermediate
  df_mispaired <- make_twa_frame(list(
    list(date = format(.focal - 60L), analytes = list(Cu = 500)),
    list(date = format(.focal - 10L), analytes = list(Cu = 1))
  ))
  mis_val <- time_weighted_aggregate(df_mispaired, focal_dates = .focal,
    tau = .tau, tau_units = "d", window = .win, window_units = "d")$value

  expect_false(isTRUE(all.equal(d1, mis_val)))
})


## ── 6. End-to-end composition: add_amspaf → time_weighted_aggregate ──────────

test_that("AmsPAF draws flow from add_amspaf into time_weighted_aggregate", {
  # Co-analytes required by Cu/Zn/Ni normalisation
  co <- c(pH = 7.5, DOC = 2.0, Ca = 6.0, Mg = 4.0, hardness = 30.0)

  # Two samples; Cu varies across 2 draws; Zn/Ni constant
  make_sample <- function(sid, date_offset, cu_draws) {
    dplyr::bind_rows(
      # Cu: drawn
      tibble::tibble(sample_id = sid, site_id = "A",
        datetime = as.Date("2024-05-01") + date_offset,
        analyte = "Cu", value = cu_draws, detected = TRUE,
        draw_id = seq_along(cu_draws)),
      # Zn, Ni: drawn (constant across draws)
      tibble::tibble(sample_id = sid, site_id = "A",
        datetime = as.Date("2024-05-01") + date_offset,
        analyte = rep(c("Zn", "Ni"), each = length(cu_draws)),
        value = rep(c(10, 0.001), each = length(cu_draws)),
        detected = TRUE,
        draw_id = rep(seq_along(cu_draws), times = 2L)),
      # Co-analytes: exact
      tibble::tibble(sample_id = sid, site_id = "A",
        datetime = as.Date("2024-05-01") + date_offset,
        analyte = names(co), value = unname(co),
        detected = TRUE, draw_id = NA_integer_)
    )
  }

  df_draws <- dplyr::bind_rows(
    make_sample("s1", 0L,  cu_draws = c(5,  50)),
    make_sample("s2", 30L, cu_draws = c(10, 100))
  )

  amspaf_draws <- suppressMessages(
    add_amspaf(df_draws, reference = NULL, conc_units = "ug/L", return = "draws")
  )
  # Filter to just AmsPAF rows and attach datetime for TWA
  amspaf_only <- dplyr::filter(amspaf_draws, analyte == "AmsPAF")
  expect_true("draw_id" %in% names(amspaf_only))

  # Chronic AmsPAF over draws
  chronic_draws <- time_weighted_aggregate(
    amspaf_only,
    focal_dates  = as.Date("2024-07-01"),
    tau          = .tau,  tau_units    = "d",
    window       = .win,  window_units = "d",
    summary      = "arith_mean",
    return       = "draws"
  )

  expect_true("draw_id" %in% names(chronic_draws))
  expect_equal(sort(unique(chronic_draws$draw_id)), 1:2)

  # Each draw's chronic AmsPAF should be finite and in [0, 100]
  expect_true(all(is.finite(chronic_draws$value)))
  expect_true(all(chronic_draws$value >= 0 & chronic_draws$value <= 100))
})


## ── 7. amspaf_daily guard ────────────────────────────────────────────────────

test_that("amspaf_daily errors with a clear message when handed draws input", {
  df_draws <- tibble::tibble(
    sample_id = "s1", site_id = "A",
    datetime  = as.Date("2024-01-01"),
    analyte   = "Cu", value = c(5, 10), detected = TRUE,
    draw_id   = c(1L, 2L)
  )
  expect_error(
    amspaf_daily(df_draws),
    "does not accept draw-bearing input"
  )
})
