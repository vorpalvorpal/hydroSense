## ============================================================================
## Shared builders for the imputation code-review TDD specs
## (test-impute-review-findings.R)
## ============================================================================
##
## These build long-format chemistry frames with the columns the imputation
## engine expects (sample_id, site_id, datetime, analyte, value, detected) and
## expose the knobs each surviving review finding needs in isolation:
##   - below-detection (BDL) predictor cells at a known detection limit,
##   - duplicate (sample, analyte) rows for a log-normal quantity,
##   - entirely-absent (all-NA) analyte columns,
##   - collision-inducing analyte names for the safe-name map.
##
## All specs here are brms/Stan-free: they exercise the deterministic
## predictor-building, hurdle and co-analyte machinery that Route C keeps, so
## they run in the default suite.

## One analyte's worth of long rows for a set of sample ids.
.imp_rows <- function(ids, analyte, value, detected = TRUE,
                      site = "A", datetime = as.Date("2023-01-01")) {
  tibble::tibble(
    sample_id = ids,
    site_id   = site,
    datetime  = datetime,
    analyte   = analyte,
    value     = value,
    detected  = detected
  )
}

## A standard WQ + metals panel. `n` samples; the WQ block drives the PCA and
## the metals are imputation targets. Deterministic given `seed`.
.imp_chem <- function(n = 40, seed = 1) {
  set.seed(seed)
  ids <- paste0("s", seq_len(n))
  dplyr::bind_rows(
    .imp_rows(ids, "pH",          stats::runif(n, 6, 8)),
    .imp_rows(ids, "EC",          stats::rlnorm(n, 6, 0.5)),
    .imp_rows(ids, "Cl",          stats::rlnorm(n, 5, 1)),
    .imp_rows(ids, "SO4²⁻", stats::rlnorm(n, 4, 1.5)),
    .imp_rows(ids, "Ca",          stats::rlnorm(n, 3, 1)),
    .imp_rows(ids, "Mg",          stats::rlnorm(n, 2.5, 0.8)),
    .imp_rows(ids, "DOC",         stats::rlnorm(n, 1, 1)),
    .imp_rows(ids, "NH3-N",       stats::rlnorm(n, 2, 1)),
    .imp_rows(ids, "temperature", stats::runif(n, 5, 20)),
    .imp_rows(ids, "Zn",          stats::rlnorm(n, 1, 1)),
    .imp_rows(ids, "Cu",          stats::rlnorm(n, 0, 1))
  )
}

## Replace an analyte column with BDL rows reported *at* a detection limit
## (value == DL, detected == FALSE) — the pattern the DL/2 fix must halve.
.imp_make_bdl <- function(df, analyte, dl) {
  df$detected[df$analyte == analyte] <- FALSE
  df$value[df$analyte == analyte]    <- dl
  df
}

## The default PCA variable set used by fit_imputation_model().
.imp_pca_vars <- function() {
  c("pH", "EC", "NH3-N", hydroSense:::.WQ_BLOCK_CANDIDATES)
}
