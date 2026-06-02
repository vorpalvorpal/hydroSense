## Regression tests for the chemistry PCA used by the imputation model.
##
## These cover the brms-independent core of fit_imputation_model(): the unified
## chemistry PCA (`.prepare_chem_pca`) and its scoring counterpart
## (`.compute_pca_scores`).  They exist because two bugs lived here undetected
## (the brms smoke tests skip when Stan is unavailable):
##   1. `cum_var <- pca_fit$R2cum` was NULL under nipals 1.0 (per-component
##      `$R2`, no `$R2cum`) → `.prepare_chem_pca` crashed / produced n_pcs = 0.
##   2. Training scores came from `nipals$scores` while prediction used the
##      regression projection in `.compute_pca_scores` — the two differ by a
##      per-component eigenvalue factor, so the brms model was trained and
##      predicted on different score scales.

library(testthat)
library(leachatetools)

make_pca_chem <- function(n = 30, seed = 1, with_missing = TRUE) {
  set.seed(seed)
  ids <- paste0("s", seq_len(n))
  mk <- function(an, vals) {
    tibble::tibble(sample_id = ids, site_id = "A",
                   datetime = as.Date("2023-01-01") + seq_len(n),
                   analyte = an, value = vals, detected = TRUE)
  }
  df <- dplyr::bind_rows(
    mk("pH",          runif(n, 6, 8)),
    mk("Cl",          rlnorm(n, 5, 1)),
    mk("SO4²⁻", rlnorm(n, 4, 1.5)),
    mk("Ca",          rlnorm(n, 3, 1)),
    mk("temperature", runif(n, 5, 20)),
    mk("DOC",         rlnorm(n, 1, 1))
  )
  if (with_missing) df <- df[-c(3, 17, 25), ]  # within-sample missing cells
  df
}

pca_vars_used <- c("pH", "Cl", "SO4²⁻", "Ca", "temperature", "DOC")


test_that(".prepare_chem_pca returns a valid variance curve (R2cum bug guard)", {
  df  <- make_pca_chem()
  pca <- leachatetools:::.prepare_chem_pca(
    df, wq_vars = pca_vars_used, min_var_explained = 0.75, max_pcs = 4L
  )
  # The R2cum-NULL bug produced n_pcs = 0 and an NA/empty var_explained.
  expect_gte(pca$n_pcs, 1L)
  expect_true(is.finite(pca$var_explained))
  expect_gte(pca$var_explained, 0)
  expect_lte(pca$var_explained, 1)
  expect_setequal(names(pca$pc_scores), c(paste0("PC", seq_len(pca$n_pcs)), "sample_id"))
})


test_that("training PC scores equal the prediction-time projection (scale bug guard)", {
  df  <- make_pca_chem()
  pca <- leachatetools:::.prepare_chem_pca(
    df, wq_vars = pca_vars_used, min_var_explained = 0.75, max_pcs = 4L
  )
  rescored <- leachatetools:::.compute_pca_scores(df, pca)

  a <- dplyr::arrange(pca$pc_scores, sample_id)
  b <- dplyr::arrange(rescored,      sample_id)
  pc_cols <- paste0("PC", seq_len(pca$n_pcs))

  # The whole point: the model is trained on `pca$pc_scores` and predicted on
  # `.compute_pca_scores()`.  They must be identical on the training data.
  expect_equal(as.matrix(a[, pc_cols]), as.matrix(b[, pc_cols]), tolerance = 1e-8)

  # And they must be on the regression-projection scale, NOT the (eigenvalue-
  # shrunk) raw nipals scores that caused the bug.
  nipals_scores <- pca$fit$scores[, seq_len(pca$n_pcs), drop = FALSE]
  expect_false(isTRUE(all.equal(
    unname(as.matrix(a[, pc_cols])), unname(nipals_scores), tolerance = 1e-6
  )))
})


test_that(".log_transform_pca logs concentrations but not pH/temperature/ORP/DO", {
  m <- matrix(
    c(7.0, 1000, 12.5, 8.0,
      6.5, 2000,  9.0, 4.0),
    nrow = 2, byrow = TRUE
  )
  colnames(m) <- c("pH", "Cl", "temperature", "DO")
  out <- leachatetools:::.log_transform_pca(m)

  expect_equal(out[, "pH"],          m[, "pH"])           # untouched
  expect_equal(out[, "temperature"], m[, "temperature"])  # untouched
  expect_equal(out[, "DO"],          m[, "DO"])           # untouched
  expect_equal(out[, "Cl"],          log10(m[, "Cl"]))    # log10'd

  # NAs are preserved (NIPALS relies on this).
  m2 <- m; m2[1, "Cl"] <- NA
  expect_true(is.na(leachatetools:::.log_transform_pca(m2)[1, "Cl"]))
})
