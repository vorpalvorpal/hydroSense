## Tests for impute_coanalytes() draw-carrier extension (Chunk 4).
## Stan-free: all GAMs use mgcv on synthetic data; nipals PCA via
## .prepare_chem_pca() (internal helper, no brms required).
##
## Eight properties tested:
##   1. Point parity    — return="point" default is byte-identical
##   2. Draws shape     — draw_id 1..N on imputed cells; NA on observed
##   3. Observed exact  — observed co-analyte rows stay draw_id=NA
##   4. N from domain   — when df already has draws, N is inferred
##   5. N required      — point-input + return="draws" without ndraws → error
##   6. ndraws mismatch — conflicting ndraws vs existing domain → error
##   7. Seed repro      — same seed → identical draws; different seeds → different
##   8. draw_domain()   — output passes ragged-N validation

library(testthat)
library(leachatetools)

## ── Shared helpers ────────────────────────────────────────────────────────────

## Build a synthetic imputation_model using the real internal .prepare_chem_pca()
## so the nipals pca_obj is guaranteed compatible with .compute_pca_scores().
##
## n_samples total; first n_observed_doc have DOC measured.
## Optionally adds metals (Cu) with n_metals_draws draw rows per sample,
## simulating output from impute_chemistry(return="draws").
make_test_setup <- function(n_samples       = 40L,
                             n_observed_doc  = 30L,
                             n_metals_draws  = 0L,
                             seed            = 42L) {
  set.seed(seed)

  all_ids <- paste0("s", seq_len(n_samples))
  dates   <- seq(as.Date("2022-01-01"), by = "week", length.out = n_samples)

  ph_vals  <- rnorm(n_samples, 7.5, 0.4)
  ec_vals  <- exp(rnorm(n_samples, log(300), 0.3))
  doc_vals <- exp(0.4 * log(ec_vals) + rnorm(n_samples, 0, 0.3))

  # Base long-format frame — pH + EC for all samples; DOC for first n_observed_doc
  base_df <- dplyr::bind_rows(
    tibble::tibble(sample_id = all_ids, site_id = "A", datetime = dates,
                   analyte = "pH",  value = ph_vals,  detected = TRUE),
    tibble::tibble(sample_id = all_ids, site_id = "A", datetime = dates,
                   analyte = "EC",  value = ec_vals,  detected = TRUE),
    tibble::tibble(sample_id = all_ids[seq_len(n_observed_doc)],
                   site_id   = "A",
                   datetime  = dates[seq_len(n_observed_doc)],
                   analyte   = "DOC",
                   value     = doc_vals[seq_len(n_observed_doc)],
                   detected  = TRUE)
  )

  # Optionally add Cu draws to simulate impute_chemistry(return="draws") output
  df <- if (n_metals_draws > 0L) {
    cu_vals <- exp(rnorm(n_samples * n_metals_draws, log(5), 0.5))
    cu_rows <- tibble::tibble(
      sample_id = rep(all_ids, each = n_metals_draws),
      site_id   = "A",
      datetime  = rep(dates,   each = n_metals_draws),
      analyte   = "Cu",
      value     = cu_vals,
      detected  = TRUE,
      draw_id   = rep(seq_len(n_metals_draws), times = n_samples)
    )
    dplyr::bind_rows(base_df, cu_rows)
  } else {
    base_df
  }

  # Build real nipals pca_obj using the internal helper — guarantees compatibility
  # with .compute_pca_scores().  pca_vars fed to the model includes "DOC" so that
  # impute_coanalytes() recognises it as an imputeable target.
  pca_obj <- leachatetools:::.prepare_chem_pca(base_df, wq_vars = c("pH", "EC"))

  model <- structure(
    list(groups   = list(),
         pca      = pca_obj,
         pca_vars = c("pH", "EC", "DOC"),
         targets  = list()),
    class = "imputation_model"
  )

  list(
    df              = df,
    model           = model,
    n_missing_doc   = n_samples - n_observed_doc,
    observed_doc_ids = all_ids[seq_len(n_observed_doc)],
    all_ids         = all_ids
  )
}


## ── 1. Point parity ──────────────────────────────────────────────────────────

test_that("return='point' default is byte-identical to explicit return='point'", {
  s <- make_test_setup()

  out_default <- suppressMessages(
    impute_coanalytes(s$df, s$model, targets = "DOC")
  )
  out_point <- suppressMessages(
    impute_coanalytes(s$df, s$model, targets = "DOC", return = "point")
  )

  expect_identical(out_default, out_point)
})

test_that("point mode: no draw_id column in output", {
  s <- make_test_setup()
  out <- suppressMessages(
    impute_coanalytes(s$df, s$model, targets = "DOC")
  )
  expect_false("draw_id" %in% names(out))
})

test_that("point mode: imputed DOC rows are present and finite", {
  s <- make_test_setup()
  out <- suppressMessages(
    impute_coanalytes(s$df, s$model, targets = "DOC")
  )
  imp_doc <- dplyr::filter(out, analyte == "DOC", imputed == TRUE)
  expect_equal(nrow(imp_doc), s$n_missing_doc)
  expect_true(all(is.finite(imp_doc$value)))
  expect_true(all(imp_doc$value > 0))
})


## ── 2. Draws shape ───────────────────────────────────────────────────────────

test_that("draws mode: imputed DOC cells carry draw_id 1..N", {
  N <- 5L
  s <- make_test_setup()
  out <- suppressMessages(
    impute_coanalytes(s$df, s$model, targets = "DOC",
                      return = "draws", ndraws = N, seed = 1L)
  )

  imp_doc <- dplyr::filter(out, analyte == "DOC", imputed == TRUE)

  expect_true("draw_id" %in% names(imp_doc))
  expect_false(anyNA(imp_doc$draw_id))
  expect_equal(sort(unique(imp_doc$draw_id)), seq_len(N))
  # Each missing sample should have exactly N rows
  expect_equal(nrow(imp_doc), s$n_missing_doc * N)
})

test_that("draws mode: imputed DOC values are positive and finite", {
  s <- make_test_setup()
  out <- suppressMessages(
    impute_coanalytes(s$df, s$model, targets = "DOC",
                      return = "draws", ndraws = 10L, seed = 2L)
  )
  imp <- dplyr::filter(out, analyte == "DOC", imputed == TRUE)
  expect_true(all(is.finite(imp$value)))
  expect_true(all(imp$value > 0))
})

test_that("draws mode: values vary across draws (not all identical)", {
  s <- make_test_setup()
  out <- suppressMessages(
    impute_coanalytes(s$df, s$model, targets = "DOC",
                      return = "draws", ndraws = 20L, seed = 3L)
  )
  # First sample that has no observed DOC (guaranteed missing)
  first_missing <- setdiff(s$all_ids, s$observed_doc_ids)[1L]
  imp <- dplyr::filter(out, analyte == "DOC", imputed == TRUE,
                        sample_id == first_missing)
  # With 20 draws there must be some variation (extremely unlikely to be all equal)
  expect_gt(length(unique(imp$value)), 1L)
})


## ── 3. Observed cells stay exact ─────────────────────────────────────────────

test_that("draws mode: observed DOC rows keep draw_id=NA (not replicated)", {
  s <- make_test_setup()
  out <- suppressMessages(
    impute_coanalytes(s$df, s$model, targets = "DOC",
                      return = "draws", ndraws = 4L, seed = 4L)
  )

  obs_doc <- dplyr::filter(out, analyte == "DOC", !imputed)
  expect_true(all(is.na(obs_doc$draw_id)))
  # One row per observed sample — not replicated N times
  expect_equal(nrow(obs_doc), length(s$observed_doc_ids))
})


## ── 4. N inferred from existing draw domain ───────────────────────────────────

test_that("draws mode: N inferred from existing metals draw domain", {
  N <- 3L
  s <- make_test_setup(n_metals_draws = N)
  out <- suppressMessages(
    impute_coanalytes(s$df, s$model, targets = "DOC", return = "draws")
  )

  imp_doc <- dplyr::filter(out, analyte == "DOC", imputed == TRUE)
  expect_equal(sort(unique(imp_doc$draw_id)), seq_len(N))
})

test_that("draws mode: ndraws=NULL accepted when domain already set", {
  N <- 4L
  s <- make_test_setup(n_metals_draws = N)
  expect_no_error(suppressMessages(
    impute_coanalytes(s$df, s$model, targets = "DOC",
                      return = "draws", ndraws = NULL)
  ))
})

test_that("draws mode: ndraws matching existing domain is accepted", {
  N <- 4L
  s <- make_test_setup(n_metals_draws = N)
  expect_no_error(suppressMessages(
    impute_coanalytes(s$df, s$model, targets = "DOC",
                      return = "draws", ndraws = N)
  ))
})


## ── 5. ndraws required for point input ───────────────────────────────────────

test_that("draws mode on point input without ndraws errors clearly", {
  s <- make_test_setup()
  expect_error(
    suppressMessages(
      impute_coanalytes(s$df, s$model, targets = "DOC", return = "draws")
    ),
    "ndraws"
  )
})


## ── 6. ndraws mismatch errors ─────────────────────────────────────────────────

test_that("ndraws conflicting with existing domain errors clearly", {
  N <- 3L
  s <- make_test_setup(n_metals_draws = N)
  expect_error(
    suppressMessages(
      impute_coanalytes(s$df, s$model, targets = "DOC",
                        return = "draws", ndraws = N + 1L)
    ),
    "conflicts"
  )
})


## ── 7. Seed reproducibility ───────────────────────────────────────────────────

test_that("same seed produces identical draws", {
  s <- make_test_setup()
  out1 <- suppressMessages(
    impute_coanalytes(s$df, s$model, targets = "DOC",
                      return = "draws", ndraws = 10L, seed = 99L)
  )
  out2 <- suppressMessages(
    impute_coanalytes(s$df, s$model, targets = "DOC",
                      return = "draws", ndraws = 10L, seed = 99L)
  )
  imp1 <- dplyr::filter(out1, analyte == "DOC", imputed == TRUE)$value
  imp2 <- dplyr::filter(out2, analyte == "DOC", imputed == TRUE)$value
  expect_identical(imp1, imp2)
})

test_that("different seeds produce different draws", {
  s <- make_test_setup()
  out1 <- suppressMessages(
    impute_coanalytes(s$df, s$model, targets = "DOC",
                      return = "draws", ndraws = 10L, seed = 1L)
  )
  out2 <- suppressMessages(
    impute_coanalytes(s$df, s$model, targets = "DOC",
                      return = "draws", ndraws = 10L, seed = 2L)
  )
  imp1 <- dplyr::filter(out1, analyte == "DOC", imputed == TRUE)$value
  imp2 <- dplyr::filter(out2, analyte == "DOC", imputed == TRUE)$value
  expect_false(isTRUE(all.equal(imp1, imp2)))
})


## ── 8. draw_domain() validates the output ─────────────────────────────────────

test_that("draws output passes .draw_domain() ragged-N validation", {
  s <- make_test_setup()
  out <- suppressMessages(
    impute_coanalytes(s$df, s$model, targets = "DOC",
                      return = "draws", ndraws = 6L, seed = 7L)
  )
  # .draw_domain() aborts on ragged N — if this passes the domain is valid
  expect_no_error(leachatetools:::.draw_domain(out))
  expect_equal(leachatetools:::.draw_domain(out), 1:6)
})
