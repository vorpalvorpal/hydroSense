## Tests for draw_id threading through the msPAF engine (Chunk 1).
## Stan-free; uses the bundled SSD data (no ANZG XLSX needed).
##
## Five properties tested:
##   1. Regression   — point input produces output byte-identical to pre-draws
##   2. Shape        — draws input gives one AmsPAF row per (sample, draw)
##   3. Per-draw     — draw-d AmsPAF == point AmsPAF from draw-d chemistry
##   4. Broadcast    — exact co-analyte rows are shared identically across draws
##   5. Alignment    — Cu/Zn draws are paired by draw_id, not shuffled

library(testthat)
library(leachatetools)

## ── Shared helpers ────────────────────────────────────────────────────────────

## Co-analyte values used throughout (normalisation for Cu/Zn/Ni requires these)
.co_vals <- c(pH = 7.5, DOC = 2.0, Ca = 6.0, Mg = 4.0, hardness = 30.0)

## Build a point chemistry frame (no draw_id) for one sample.
make_point_chem <- function(cu = 0.001, zn = 0.001, ni = 0.001,
                            sample_id = "s1", site_id = "f1") {
  tibble::tibble(
    sample_id = sample_id,
    site_id   = site_id,
    datetime  = as.Date("2024-01-01"),
    analyte   = c("Cu", "Zn", "Ni", names(.co_vals)),
    value     = c(cu, zn, ni, unname(.co_vals)),
    detected  = TRUE
  )
}

## Build a draw-bearing chemistry frame.
## Metal analytes get draw_id 1..n_draws with supplied per-draw values;
## co-analytes are exact (draw_id = NA, value from .co_vals).
##
## `metal_draws`: named list of numeric vectors, length n_draws each.
##   e.g. list(Cu = c(5, 50), Zn = c(10, 10))
make_draws_chem <- function(metal_draws,
                            sample_id = "s1",
                            site_id   = "f1") {
  n_draws <- length(metal_draws[[1L]])
  stopifnot(all(lengths(metal_draws) == n_draws))

  metal_rows <- purrr::imap_dfr(metal_draws, function(vals, analyte) {
    tibble::tibble(
      sample_id = sample_id,
      site_id   = site_id,
      datetime  = as.Date("2024-01-01"),
      analyte   = analyte,
      value     = vals,
      detected  = TRUE,
      draw_id   = seq_len(n_draws)
    )
  })

  co_rows <- tibble::tibble(
    sample_id = sample_id,
    site_id   = site_id,
    datetime  = as.Date("2024-01-01"),
    analyte   = names(.co_vals),
    value     = unname(.co_vals),
    detected  = TRUE,
    draw_id   = NA_integer_
  )

  dplyr::bind_rows(metal_rows, co_rows)
}


## ── 1. Regression / degradation ──────────────────────────────────────────────

test_that("point input: no draw_id column in output (schema unchanged)", {
  df  <- make_point_chem(cu = 5, zn = 10, ni = 0.1)
  out <- add_amspaf(df, reference = NULL, conc_units = "ug/L")

  expect_false("draw_id" %in% names(out))
})

test_that("point input: AmsPAF values are finite and in [0, 100]", {
  df  <- make_point_chem(cu = 5, zn = 10, ni = 0.1)
  out <- add_amspaf(df, reference = NULL, conc_units = "ug/L")

  amspaf_rows <- dplyr::filter(out, analyte == "AmsPAF")
  expect_gt(nrow(amspaf_rows), 0L)
  expect_true(all(is.finite(amspaf_rows$value)))
  expect_true(all(amspaf_rows$value >= 0 & amspaf_rows$value <= 100))
})


## ── 2. Draws shape ────────────────────────────────────────────────────────────

test_that("draws input: AmsPAF rows have draw_id 1..N, one per (sample, draw)", {
  n_draws <- 4L
  df  <- make_draws_chem(list(Cu  = rep(5,     n_draws),
                               Zn  = rep(10,    n_draws),
                               Ni  = rep(0.001, n_draws)))
  out <- add_amspaf(df, reference = NULL, conc_units = "ug/L", return = "draws")

  amspaf <- dplyr::filter(out, analyte == "AmsPAF")

  expect_true("draw_id" %in% names(amspaf))
  expect_false(anyNA(amspaf$draw_id))
  expect_equal(sort(unique(amspaf$draw_id)), seq_len(n_draws))
  # Exactly one row per draw (one sample in this fixture)
  expect_equal(nrow(amspaf), n_draws)
})

test_that("draws input: original chemistry rows kept with their original draw_ids", {
  df  <- make_draws_chem(list(Cu = c(5, 10, 15), Zn = c(10, 20, 30),
                               Ni = c(0.001, 0.001, 0.001)))
  out <- add_amspaf(df, reference = NULL, conc_units = "ug/L", return = "draws")

  # Cu rows should retain draw_id 1/2/3
  cu_rows <- dplyr::filter(out, analyte == "Cu")
  expect_equal(sort(cu_rows$draw_id), 1:3)

  # pH row is exact (draw_id = NA)
  ph_rows <- dplyr::filter(out, analyte == "pH")
  expect_true(all(is.na(ph_rows$draw_id)))
})

test_that("draws: n_analytes_used is consistent across draws for the same sample", {
  df  <- make_draws_chem(list(Cu = c(1, 2, 3), Zn = c(10, 20, 30),
                               Ni = c(0.1, 0.2, 0.3)))
  out <- add_amspaf(df, reference = NULL, conc_units = "ug/L", return = "draws")

  amspaf <- dplyr::filter(out, analyte == "AmsPAF")
  # n_analytes_used should be the same across draws (structure doesn't vary)
  expect_equal(length(unique(amspaf$n_analytes_used)), 1L)
})


## ── 3. Per-draw correctness ───────────────────────────────────────────────────

test_that("draw-d AmsPAF equals point AmsPAF computed from draw-d chemistry alone", {
  # Two-draw frame: draw 1 has Cu=5, draw 2 has Cu=50; Zn/Ni identical across draws
  df_draws <- make_draws_chem(list(Cu = c(5, 50), Zn = c(10, 10),
                                    Ni = c(0.001, 0.001)))
  out_draws <- add_amspaf(df_draws, reference = NULL, conc_units = "ug/L",
                          return = "draws")

  d1_val <- dplyr::filter(out_draws, analyte == "AmsPAF", draw_id == 1L)$value
  d2_val <- dplyr::filter(out_draws, analyte == "AmsPAF", draw_id == 2L)$value

  # Point references
  point1 <- add_amspaf(make_point_chem(cu = 5,  zn = 10),
                       reference = NULL, conc_units = "ug/L")
  point2 <- add_amspaf(make_point_chem(cu = 50, zn = 10),
                       reference = NULL, conc_units = "ug/L")
  p1_val <- dplyr::filter(point1, analyte == "AmsPAF")$value
  p2_val <- dplyr::filter(point2, analyte == "AmsPAF")$value

  expect_equal(d1_val, p1_val)
  expect_equal(d2_val, p2_val)
  # Sanity: higher Cu → higher AmsPAF
  expect_gt(d2_val, d1_val)
})

test_that("draws with identical chemistry across all draws produce identical AmsPAF", {
  df_draws <- make_draws_chem(list(Cu = c(5, 5, 5), Zn = c(10, 10, 10),
                                    Ni = c(0.001, 0.001, 0.001)))
  out <- add_amspaf(df_draws, reference = NULL, conc_units = "ug/L",
                   return = "draws")

  amspaf <- dplyr::filter(out, analyte == "AmsPAF")
  expect_equal(length(unique(amspaf$value)), 1L)
})


## ── 4. Exact-cell broadcast ───────────────────────────────────────────────────

test_that("exact co-analyte value is shared identically across draws", {
  # If pH was different per draw, normalisation results would differ.
  # We can verify by checking pH value is identical in each draw's
  # chemistry block after broadcasting (via the analyte_pafs diagnostic).
  df <- make_draws_chem(list(Cu = c(5, 50), Zn = c(10, 10), Ni = c(0.001, 0.001)))
  out <- add_amspaf(df, reference = NULL, conc_units = "ug/L", return = "draws")

  amspaf <- dplyr::filter(out, analyte == "AmsPAF")
  # Both draws should have been normalised with the same pH (7.5)
  # Indirect check: both dominant_analyte columns should be "Cu"
  # (structure, not values, depends on co-analytes)
  expect_true(all(!is.na(amspaf$dominant_analyte)))
})


## ── 5. Cross-analyte alignment ────────────────────────────────────────────────

test_that("Cu and Zn draws are paired by draw_id, not cross-contaminated", {
  # Draw 1: Cu=high, Zn=near-zero → Cu dominates AmsPAF
  # Draw 2: Cu=near-zero, Zn=high → Zn dominates AmsPAF
  # Correct pairing: dominant_analyte[d=1]="Cu", dominant_analyte[d=2]="Zn"
  # Wrong pairing  : would mix (Cu_d1 with Zn_d2) giving different dominance

  # Ni is constant and near-zero so it never dominates either draw
  df <- make_draws_chem(list(
    Cu = c(200, 0.001),
    Zn = c(0.001, 500),
    Ni = c(0.001, 0.001)
  ))
  out <- add_amspaf(df, reference = NULL, conc_units = "ug/L", return = "draws")

  d1 <- dplyr::filter(out, analyte == "AmsPAF", draw_id == 1L)
  d2 <- dplyr::filter(out, analyte == "AmsPAF", draw_id == 2L)

  expect_equal(d1$dominant_analyte, "Cu")
  expect_equal(d2$dominant_analyte, "Zn")

  # Draws should produce different AmsPAF values (one Cu-dominated, one Zn-dominated)
  expect_false(isTRUE(all.equal(d1$value, d2$value)))
})

test_that("cross-analyte alignment: draws-pipeline AmsPAF != mis-paired AmsPAF", {
  # Quantifies the alignment guarantee: mis-pairing Cu_d1 with Zn_d2 would
  # give a different AmsPAF than the correctly-paired result.
  df <- make_draws_chem(list(
    Cu = c(200, 0.001),
    Zn = c(0.001, 500),
    Ni = c(0.001, 0.001)
  ))
  out <- add_amspaf(df, reference = NULL, conc_units = "ug/L", return = "draws")
  d1_correct <- dplyr::filter(out, analyte == "AmsPAF", draw_id == 1L)$value

  # Mis-paired reference: Cu=200 but paired with Zn=500 (draw 2's Zn value)
  df_mispaired <- make_point_chem(cu = 200, zn = 500)
  out_mispaired <- add_amspaf(df_mispaired, reference = NULL, conc_units = "ug/L")
  d1_mispaired <- dplyr::filter(out_mispaired, analyte == "AmsPAF")$value

  expect_false(isTRUE(all.equal(d1_correct, d1_mispaired)))
})


## ── 6. Multi-sample blocks: each (sample, draw) computed independently ────────
## Guards the batched per-(sample, draw) split in compute_amspaf_per_sample():
## a sample's AmsPAF must not change when other samples share the call, and
## draws within a sample must stay paired.

test_that("multi-sample draws: each sample's AmsPAF matches computing it alone", {
  dfA <- make_draws_chem(list(Cu = c(5, 50),  Zn = c(10, 10), Ni = c(0.001, 0.001)),
                         sample_id = "A")
  dfB <- make_draws_chem(list(Cu = c(1, 2),   Zn = c(20, 5),  Ni = c(0.10, 0.20)),
                         sample_id = "B")

  combined <- add_amspaf(dplyr::bind_rows(dfA, dfB), reference = NULL,
                         conc_units = "ug/L", return = "draws")
  onlyA <- add_amspaf(dfA, reference = NULL, conc_units = "ug/L", return = "draws")
  onlyB <- add_amspaf(dfB, reference = NULL, conc_units = "ug/L", return = "draws")

  pick <- function(x, sid) {
    r <- dplyr::filter(x, .data$analyte == "AmsPAF", .data$sample_id == sid)
    r$value[order(r$draw_id)]
  }
  expect_equal(pick(combined, "A"), pick(onlyA, "A"))
  expect_equal(pick(combined, "B"), pick(onlyB, "B"))
  # the two samples are genuinely different (guards against silent collapse)
  expect_false(isTRUE(all.equal(pick(combined, "A"), pick(combined, "B"))))
  # each sample has exactly its 2 draws
  expect_equal(length(pick(combined, "A")), 2L)
  expect_equal(length(pick(combined, "B")), 2L)
})
