## Stan-free tests for the draw-carrier primitives (R/draws.R).
## No brms/Stan — pure data-shape checks.

library(testthat)
library(hydroSense)

# ── Helpers ───────────────────────────────────────────────────────────────────

# Minimal long frame: draw-bearing cells for analytes, exact cell for pH.
make_mixed_frame <- function(n_samples = 2L,
                             n_draws   = 3L,
                             analytes  = c("Cu", "Zn"),
                             seed      = 1L) {
  set.seed(seed)
  draw_rows <- tidyr::expand_grid(
    sample_id = paste0("s", seq_len(n_samples)),
    analyte   = analytes,
    draw_id   = seq_len(n_draws)
  ) |>
    dplyr::mutate(
      site_id      = "A",
      datetime     = as.Date("2024-01-01"),
      value        = stats::rnorm(dplyr::n(), 10, 2),
      detected     = TRUE,
      imputed      = TRUE,
      imputed_kind = "censored_left",
      draw_id      = as.integer(.data$draw_id)
    )

  exact_rows <- tidyr::expand_grid(
    sample_id = paste0("s", seq_len(n_samples)),
    analyte   = "pH"
  ) |>
    dplyr::mutate(
      site_id      = "A",
      datetime     = as.Date("2024-01-01"),
      value        = stats::runif(dplyr::n(), 6, 8),
      detected     = TRUE,
      imputed      = FALSE,
      imputed_kind = "observed",
      draw_id      = NA_integer_
    )

  dplyr::bind_rows(draw_rows, exact_rows)
}

# All-exact frame (no draw_id column at all).
make_exact_frame <- function(n_samples = 3L) {
  tidyr::expand_grid(
    sample_id = paste0("s", seq_len(n_samples)),
    analyte   = c("Cu", "pH")
  ) |>
    dplyr::mutate(
      site_id  = "A",
      datetime = as.Date("2024-01-01"),
      value    = stats::runif(dplyr::n(), 1, 10),
      detected = TRUE
    )
}


# ── .draw_domain ─────────────────────────────────────────────────────────────

test_that(".draw_domain returns integer(0) when draw_id column absent", {
  df <- make_exact_frame()
  expect_identical(hydroSense:::.draw_domain(df), integer(0))
})

test_that(".draw_domain returns integer(0) when all draw_ids are NA", {
  df <- make_exact_frame() |>
    dplyr::mutate(draw_id = NA_integer_)
  expect_identical(hydroSense:::.draw_domain(df), integer(0))
})

test_that(".draw_domain returns 1:N for a valid draws frame", {
  df <- make_mixed_frame(n_draws = 4L)
  expect_identical(hydroSense:::.draw_domain(df), 1:4)
})

test_that(".draw_domain aborts on non-contiguous draw IDs", {
  df <- make_mixed_frame(n_draws = 3L) |>
    dplyr::filter(.data$draw_id != 2L)   # remove draw 2 → gap
  expect_error(hydroSense:::.draw_domain(df), "contiguous")
})

test_that(".draw_domain aborts on ragged N across cells", {
  df <- make_mixed_frame(n_draws = 3L, analytes = c("Cu", "Zn"))
  # Drop draw 3 for Cu only → Cu has 2 draws, Zn has 3
  bad <- dplyr::filter(df,
    !(.data$analyte == "Cu" & .data$draw_id == 3L)
  )
  expect_error(hydroSense:::.draw_domain(bad), "[Rr]agged")
})

test_that(".draw_domain ragged check is skipped when sample_id absent", {
  # A frame with only draw_id and value (no sample_id / analyte)
  df <- tibble::tibble(draw_id = c(1L, 1L, 2L), value = 1:3)
  expect_identical(hydroSense:::.draw_domain(df), 1:2)
})


# ── .broadcast_draws ─────────────────────────────────────────────────────────

test_that(".broadcast_draws replicates exact cells N times with correct draw_ids", {
  # make_mixed_frame produces Cu draw cells + pH exact cell for 1 sample
  df <- make_mixed_frame(n_samples = 1L, n_draws = 3L, analytes = "Cu")

  out <- hydroSense:::.broadcast_draws(df)

  # All draw_ids should now be non-NA
  expect_false(anyNA(out$draw_id))

  # pH (1 exact row) should appear 3 times (one per draw)
  ph_rows <- dplyr::filter(out, .data$analyte == "pH")
  expect_equal(nrow(ph_rows), 3L)
  expect_equal(sort(ph_rows$draw_id), 1:3)

  # pH value identical across draws (same observation replicated)
  expect_equal(length(unique(ph_rows$value)), 1L)

  # Cu (draw cells) still has 3 rows, unchanged
  cu_rows <- dplyr::filter(out, .data$analyte == "Cu")
  expect_equal(nrow(cu_rows), 3L)
})

test_that(".broadcast_draws leaves draw cells unchanged (order + values)", {
  df <- make_mixed_frame(n_samples = 2L, n_draws = 3L, analytes = c("Cu", "Zn"))
  draws_before <- dplyr::filter(df, !is.na(.data$draw_id)) |>
    dplyr::arrange(.data$sample_id, .data$analyte, .data$draw_id)

  out <- hydroSense:::.broadcast_draws(df)

  draws_after <- dplyr::filter(out, .data$analyte %in% c("Cu", "Zn")) |>
    dplyr::arrange(.data$sample_id, .data$analyte, .data$draw_id)

  expect_equal(draws_after$value,   draws_before$value)
  expect_equal(draws_after$draw_id, draws_before$draw_id)
})

test_that(".broadcast_draws handles a deterministic (no draws) frame", {
  df  <- make_exact_frame()
  out <- hydroSense:::.broadcast_draws(df)

  expect_false(anyNA(out$draw_id))
  expect_true(all(out$draw_id == 1L))
  expect_equal(nrow(out), nrow(df))    # no row expansion
})

test_that(".broadcast_draws is idempotent", {
  df   <- make_mixed_frame(n_samples = 1L, n_draws = 2L)
  out1 <- hydroSense:::.broadcast_draws(df)
  out2 <- hydroSense:::.broadcast_draws(out1)

  expect_equal(nrow(out1), nrow(out2))
  expect_equal(sort(out1$draw_id), sort(out2$draw_id))
})

test_that(".broadcast_draws round-trip: draw cell (analyte, draw_id, value) map is bit-identical", {
  df <- make_mixed_frame(n_samples = 2L, n_draws = 4L, analytes = c("Cu", "Zn"))

  key_before <- df |>
    dplyr::filter(!is.na(.data$draw_id)) |>
    dplyr::arrange(.data$sample_id, .data$analyte, .data$draw_id) |>
    dplyr::select("sample_id", "analyte", "draw_id", "value")

  key_after <- hydroSense:::.broadcast_draws(df) |>
    dplyr::filter(.data$analyte %in% c("Cu", "Zn")) |>
    dplyr::arrange(.data$sample_id, .data$analyte, .data$draw_id) |>
    dplyr::select("sample_id", "analyte", "draw_id", "value")

  expect_identical(key_before, key_after)
})

test_that(".broadcast_draws degradation: all-observed frame produces deterministic output", {
  # A fully-observed, no-error monitoring frame (no draw_id column) should
  # produce N=1 pass identical in content to the point-estimate pipeline.
  exact <- make_exact_frame(n_samples = 4L)

  out <- hydroSense:::.broadcast_draws(exact)

  expect_equal(nrow(out), nrow(exact))
  expect_true(all(out$draw_id == 1L))
  expect_equal(out$value, exact$value)
})

test_that(".broadcast_draws works on a frame with only exact cells and draw_id NA col", {
  df <- make_exact_frame() |>
    dplyr::mutate(draw_id = NA_integer_)

  out <- hydroSense:::.broadcast_draws(df)
  expect_true(all(out$draw_id == 1L))
  expect_equal(nrow(out), nrow(df))
})
