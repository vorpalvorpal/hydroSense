## Tests for summarise_draws() and the return= param at the boundary
## functions (add_amspaf, time_weighted_aggregate).  Stan-free.
##
## Eleven properties tested:
##   1.  summarise_draws point identity     — no CI cols added to point frames
##   2.  Quantile correctness               — known draws, exact median + 5th/95th
##   3.  Exact cell degenerate interval     — draw_id=NA → lower=upper=value, n=1
##   4.  interval configurable              — 0.5 strictly narrower than 0.9
##   5.  central="mean" vs "median"         — differ on skewed cell
##   6.  add_amspaf(return="summary")       — AmsPAF rows have value_lower/upper, no draw_id
##   7.  add_amspaf(return="draws")         — raw per-draw rows (Chunk 1 parity)
##   8.  add_amspaf point + default         — byte-identical output (no CI cols)
##   9.  TWA(return="summary")              — collapses chronic draws → median+CI
##   10. End-to-end composition             — draws|>add_amspaf(draws)|>TWA() → CI cols, no draw_id
##   11. Median within interval             — value ∈ [value_lower, value_upper]

library(testthat)
library(hydroSense)

## ── Shared helpers ────────────────────────────────────────────────────────────

.co_vals <- c(pH = 7.5, DOC = 2.0, Ca = 6.0, Mg = 4.0, hardness = 30.0)

## Minimal draw-carrier frame: Cu drawn, co-analytes exact.
make_amspaf_draws_input <- function(cu_draws = c(5, 20, 50),
                                    sample_id = "s1",
                                    site_id   = "f1") {
  n <- length(cu_draws)
  dplyr::bind_rows(
    tibble::tibble(sample_id = sample_id, site_id = site_id,
      datetime = as.Date("2024-01-01"), analyte = "Cu",
      value = cu_draws, detected = TRUE, draw_id = seq_len(n)),
    tibble::tibble(sample_id = sample_id, site_id = site_id,
      datetime = as.Date("2024-01-01"), analyte = "Zn",
      value = rep(10, n), detected = TRUE, draw_id = seq_len(n)),
    tibble::tibble(sample_id = sample_id, site_id = site_id,
      datetime = as.Date("2024-01-01"), analyte = "Ni",
      value = rep(0.001, n), detected = TRUE, draw_id = seq_len(n)),
    tibble::tibble(sample_id = sample_id, site_id = site_id,
      datetime = as.Date("2024-01-01"), analyte = names(.co_vals),
      value = unname(.co_vals), detected = TRUE, draw_id = NA_integer_)
  )
}

## Minimal point chemistry frame (no draw_id).
make_point_input <- function(cu = 5, sample_id = "s1", site_id = "f1") {
  tibble::tibble(
    sample_id = sample_id, site_id = site_id,
    datetime  = as.Date("2024-01-01"),
    analyte   = c("Cu", "Zn", "Ni", names(.co_vals)),
    value     = c(cu, 10, 0.001, unname(.co_vals)),
    detected  = TRUE
  )
}


## ── 1. summarise_draws point identity ────────────────────────────────────────

test_that("summarise_draws is identity on a point frame", {
  df  <- make_point_input()
  out <- summarise_draws(df)

  expect_identical(out, df)
  expect_false("value_lower" %in% names(out))
  expect_false("n_draws"     %in% names(out))
})

test_that("summarise_draws is identity on all-NA draw_id frame", {
  df  <- make_point_input() |> dplyr::mutate(draw_id = NA_integer_)
  out <- summarise_draws(df)
  expect_identical(out, df)
})


## ── 2. Quantile correctness ───────────────────────────────────────────────────

test_that("summarise_draws returns correct median and 5th/95th percentiles", {
  # Known 100-draw cell: draws = 1:100
  df <- tibble::tibble(
    sample_id = "s1", site_id = "A", analyte = "Cu",
    value   = 1:100,
    draw_id = 1:100,
    detected = TRUE
  )
  out <- summarise_draws(df, interval = 0.90)

  expect_equal(out$value,       stats::median(1:100))
  expect_equal(out$value_lower, stats::quantile(1:100, 0.05, names = FALSE))
  expect_equal(out$value_upper, stats::quantile(1:100, 0.95, names = FALSE))
  expect_equal(out$n_draws, 100L)
})


## ── 3. Exact cell degenerate interval ────────────────────────────────────────

test_that("exact cells (draw_id=NA) produce degenerate interval", {
  # Mix: Cu drawn (3 draws), pH exact (draw_id=NA).
  # summarise_draws should degenerate pH to lower=upper=value, n_draws=1.
  df <- dplyr::bind_rows(
    tibble::tibble(sample_id="s1", site_id="A", analyte="Cu",
                   value=c(5,10,15), draw_id=1:3, detected=TRUE),
    tibble::tibble(sample_id="s1", site_id="A", analyte="pH",
                   value=7.5, draw_id=NA_integer_, detected=TRUE)
  )
  out <- summarise_draws(df)

  ph <- dplyr::filter(out, analyte == "pH")
  expect_equal(ph$value,       7.5)
  expect_equal(ph$value_lower, 7.5)
  expect_equal(ph$value_upper, 7.5)
  expect_equal(ph$n_draws, 1L)
})


## ── 4. interval configurable ─────────────────────────────────────────────────

test_that("interval=0.5 gives strictly narrower CI than interval=0.9", {
  df <- tibble::tibble(
    sample_id = "s1", site_id = "A", analyte = "Cu",
    value = 1:100, draw_id = 1:100, detected = TRUE
  )
  out90 <- summarise_draws(df, interval = 0.90)
  out50 <- summarise_draws(df, interval = 0.50)

  width90 <- out90$value_upper - out90$value_lower
  width50 <- out50$value_upper - out50$value_lower
  expect_gt(width90, width50)
})


## ── 5. central="mean" vs "median" ────────────────────────────────────────────

test_that("central=mean and median differ on a right-skewed cell", {
  # Skewed draws: 1..9 plus one large outlier
  df <- tibble::tibble(
    sample_id = "s1", site_id = "A", analyte = "Cu",
    value = c(1:9, 1000), draw_id = 1:10, detected = TRUE
  )
  med_val  <- summarise_draws(df, central = "median")$value
  mean_val <- summarise_draws(df, central = "mean")$value

  expect_false(isTRUE(all.equal(med_val, mean_val)))
  expect_gt(mean_val, med_val)  # mean pulled up by outlier
})


## ── 6. add_amspaf(return="summary") default ──────────────────────────────────

test_that("add_amspaf default collapses AmsPAF draws to median+CI, one row per sample", {
  df   <- make_amspaf_draws_input(cu_draws = c(5, 20, 50))
  out  <- suppressMessages(
    add_amspaf(df, reference = NULL, conc_units = "ug/L")
  )

  amspaf <- dplyr::filter(out, analyte == "AmsPAF")

  # One row per sample (collapsed)
  expect_equal(nrow(amspaf), 1L)

  # CI columns present
  expect_true("value_lower" %in% names(amspaf))
  expect_true("value_upper" %in% names(amspaf))
  expect_true("n_draws"     %in% names(amspaf))

  # No raw-draws columns
  expect_false("draw_id"          %in% names(amspaf))
  expect_false("dominant_analyte" %in% names(amspaf))
  expect_false("analyte_pafs"     %in% names(amspaf))

  # n_draws = 3
  expect_equal(amspaf$n_draws, 3L)

  # value ∈ [lower, upper]
  expect_gte(amspaf$value, amspaf$value_lower)
  expect_lte(amspaf$value, amspaf$value_upper)
})


## ── 7. add_amspaf(return="draws") passes through raw draws ───────────────────

test_that("add_amspaf(return='draws') emits raw per-draw AmsPAF rows", {
  df  <- make_amspaf_draws_input(cu_draws = c(5, 20, 50))
  out <- suppressMessages(
    add_amspaf(df, reference = NULL, conc_units = "ug/L", return = "draws")
  )

  amspaf <- dplyr::filter(out, analyte == "AmsPAF")

  expect_equal(nrow(amspaf), 3L)
  expect_true("draw_id" %in% names(amspaf))
  expect_equal(sort(amspaf$draw_id), 1:3)
  expect_false("value_lower" %in% names(amspaf))
})


## ── 8. add_amspaf point input + default = byte-identical to pre-draws ────────

test_that("add_amspaf default on point input: no CI columns, schema unchanged", {
  df  <- make_point_input(cu = 5)
  out <- suppressMessages(
    add_amspaf(df, reference = NULL, conc_units = "ug/L")
  )

  expect_false("draw_id"    %in% names(out))
  expect_false("value_lower" %in% names(out))
  expect_false("n_draws"    %in% names(out))

  amspaf <- dplyr::filter(out, analyte == "AmsPAF")
  expect_equal(nrow(amspaf), 1L)
  expect_true(all(is.finite(amspaf$value)))
})


## ── 9. time_weighted_aggregate(return="summary") collapses draws ─────────────

test_that("TWA(return='summary') collapses chronic AmsPAF draws to median+CI", {
  # Build AmsPAF draws for two samples
  make_s <- function(sid, dt_offset, cu_draws) {
    n <- length(cu_draws)
    dplyr::bind_rows(
      tibble::tibble(sample_id=sid, site_id="A",
        datetime=as.Date("2024-01-01")+dt_offset,
        analyte="Cu", value=cu_draws, detected=TRUE, draw_id=seq_len(n)),
      tibble::tibble(sample_id=sid, site_id="A",
        datetime=as.Date("2024-01-01")+dt_offset,
        analyte=c("Zn","Ni"), value=c(10,0.001), detected=TRUE,
        draw_id=rep(NA_integer_,2L)),
      tibble::tibble(sample_id=sid, site_id="A",
        datetime=as.Date("2024-01-01")+dt_offset,
        analyte=names(.co_vals), value=unname(.co_vals),
        detected=TRUE, draw_id=NA_integer_)
    )
  }
  df <- dplyr::bind_rows(
    make_s("s1", 0L,  c(5, 10, 20)),
    make_s("s2", 30L, c(5, 10, 20))
  )
  amspaf_draws <- suppressMessages(
    add_amspaf(df, reference = NULL, conc_units = "ug/L", return = "draws")
  )
  amspaf_only <- dplyr::filter(amspaf_draws, analyte == "AmsPAF")

  chronic <- time_weighted_aggregate(
    amspaf_only,
    focal_dates  = as.Date("2024-03-01"),
    tau = 90, tau_units = "d", window = 365, window_units = "d",
    summary = "arith_mean",
    return  = "summary"
  )

  expect_true("value_lower" %in% names(chronic))
  expect_true("value_upper" %in% names(chronic))
  expect_true("n_draws"     %in% names(chronic))
  expect_false("draw_id"    %in% names(chronic))

  expect_equal(chronic$n_draws, 3L)
  expect_true(all(is.finite(chronic$value)))
})


## ── 10. End-to-end composition: draws → add_amspaf(draws) → TWA() ────────────

test_that("draws → add_amspaf(draws) → TWA default → CI cols, no draw_id", {
  df <- make_amspaf_draws_input(cu_draws = c(5, 20, 50))
  amspaf_raw <- suppressMessages(
    add_amspaf(df, reference = NULL, conc_units = "ug/L", return = "draws")
  )
  amspaf_only <- dplyr::filter(amspaf_raw, analyte == "AmsPAF")

  chronic <- time_weighted_aggregate(
    amspaf_only,
    focal_dates  = as.Date("2024-03-01"),
    tau = 90, tau_units = "d", window = 365, window_units = "d",
    summary = "arith_mean"
    # return defaults to "summary"
  )

  expect_true("value_lower"  %in% names(chronic))
  expect_true("value_upper"  %in% names(chronic))
  expect_false("draw_id"     %in% names(chronic))
  expect_true(all(is.finite(chronic$value)))
})


## ── 11. Median within interval ────────────────────────────────────────────────

test_that("value (median) lies within [value_lower, value_upper] for all cells", {
  df   <- make_amspaf_draws_input(cu_draws = c(5, 20, 50, 100, 200))
  out  <- suppressMessages(
    add_amspaf(df, reference = NULL, conc_units = "ug/L")
  )

  amspaf <- dplyr::filter(out, analyte == "AmsPAF")
  expect_true(all(amspaf$value >= amspaf$value_lower - .Machine$double.eps))
  expect_true(all(amspaf$value <= amspaf$value_upper + .Machine$double.eps))
})
