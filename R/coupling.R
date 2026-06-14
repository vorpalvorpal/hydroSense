## ============================================================================
## coupling.R -- cross-analyte correlation of Kalman/OU residuals (issue #32)
## ============================================================================
##
## Estimates the empirical cross-analyte correlation of the Kalman/OU anchor
## residuals (the `S` column in each analyte's anchor data), regularises it,
## and returns a valid positive-definite correlation matrix for use as the
## coupling seed in the draw chokepoint.
##
## Public surface (package-internal):
##   .anchor_residual_cor   pairwise-complete Pearson R, ridge-shrunk + nearPD

#' @importFrom Matrix nearPD
NULL


## ── Stage 1: empirical cross-analyte correlation ──────────────────────────────

#' Estimate a regularised cross-analyte correlation matrix from anchor residuals
#'
#' For each analyte in `analytes`, the anchor residual column `S` (from
#' `d_anchors` if available and >= 2 rows, otherwise `anchors`) is extracted and
#' joined into a wide `[date × analyte]` matrix.  Pairwise Pearson correlation is
#' computed, ridge-shrunk toward the identity, and the result is projected onto
#' the positive-definite cone via [Matrix::nearPD()].
#'
#' ## Statistical choices
#'
#' **Pairwise complete observations** (`use = "pairwise.complete.obs"`):
#' The ARA analyte panel is ragged across sampling cadences (e.g. B.S01:
#' NH3-N on 121 grabs, metals on 78–79, Al on 7; complete-case would be 0
#' rows).  Pairwise estimation uses every available co-measured date for each
#' pair, maximising the information extracted from the sparse overlap.
#'
#' **Ridge shrinkage (fixed λ = 0.1)**: pulls low-*n* pairs toward independence
#' (conservative when co-observed dates are sparse).  This is more appropriate
#' here than analytic Schäfer–Strimmer shrinkage, which requires a
#' complete-case matrix.  The shrinkage formula is
#' `R_ridge = (1 - λ) * R_hat + λ * I`.
#'
#' **`Matrix::nearPD`** (Higham 2002, alternating projection algorithm):
#' Pairwise estimates with unequal *n* per pair can produce indefinite matrices.
#' nearPD projects onto the nearest (in Frobenius norm) positive-definite
#' correlation matrix.  Only applied when the minimum eigenvalue of `R_ridge`
#' is <= 1e-8; in the typical case (ridge alone suffices) the result is returned
#' directly without the projection step.
#'
#' **Degenerate analytes** (constant `S` or < 2 finite observations): forced
#' independent (identity row/col in `R_hat`) rather than propagating NA.
#'
#' @param tm A target-model list with element `models`, a named list of
#'   per-analyte model objects.  Each model must have an `anchors` data frame
#'   with `date` and `S` columns, and optionally a `d_anchors` data frame with
#'   the same structure.
#' @param analytes Character vector of analyte names (must be names of
#'   `tm$models`).
#'
#' @return A named list with three elements:
#'   - `R`        – a `[p × p]` positive-definite correlation matrix (numeric
#'                  matrix, dimnames = `list(analytes, analytes)`).
#'   - `analytes` – the `analytes` argument, preserved for downstream use.
#'   - `lambda`   – the ridge hyperparameter (numeric scalar, currently 0.1).
#'
#' @references
#'   Higham, N. J. (2002). Computing the nearest correlation matrix — a problem
#'   from finance. *IMA Journal of Numerical Analysis*, **22**(3), 329–343.
#'   <https://doi.org/10.1093/imanum/22.3.329>
#'
#' @keywords internal
.anchor_residual_cor <- function(tm, analytes) {
  p <- length(analytes)

  ## ── 1. Extract anchor residuals per analyte ──────────────────────────────

  ## For each analyte pull the (date, S) data frame, using d_anchors when it
  ## has >= 2 rows, otherwise falling back to anchors.
  extract_anch <- function(nm) {
    m    <- tm$models[[nm]]
    da   <- m$d_anchors
    anch <- if (!is.null(da) && nrow(da) >= 2L) da else m$anchors
    if (is.null(anch) || nrow(anch) == 0L) return(NULL)
    anch[, c("date", "S")]
  }

  anch_list <- stats::setNames(
    lapply(analytes, extract_anch),
    analytes
  )

  ## ── 2. Flag degenerate analytes before the join ──────────────────────────

  ## An analyte is degenerate if: NULL anchors, < 2 finite S values, or
  ## (near-)constant S (sd == 0 up to machine precision).
  is_degenerate <- function(nm) {
    df <- anch_list[[nm]]
    if (is.null(df)) return(TRUE)
    s  <- df$S[is.finite(df$S)]
    if (length(s) < 2L) return(TRUE)
    ## stats::var() returns NA for length-1; for length-0 this can't reach here.
    v <- stats::var(s)
    !is.finite(v) || v <= .Machine$double.eps
  }

  degen <- stats::setNames(vapply(analytes, is_degenerate, logical(1L)), analytes)

  ## Single-analyte shortcut — avoids the join entirely.
  if (p == 1L) {
    R <- matrix(1, 1L, 1L, dimnames = list(analytes, analytes))
    return(list(R = R, analytes = analytes, lambda = 0.1))
  }

  ## ── 3. Assemble wide matrix [date × analyte] ─────────────────────────────

  ## Build a list of two-column data frames suitable for repeated merge().
  ## Non-degenerate analytes contribute their S; degenerate ones are left out
  ## (they contribute only NA after the join, which is fine — the degenerate
  ## flag already controls the downstream handling).
  make_df <- function(nm) {
    df  <- anch_list[[nm]]
    out <- data.frame(date = as.Date(character(0L)), S = numeric(0L))
    if (!is.null(df)) {
      out <- data.frame(date = as.Date(df$date), S = df$S)
    }
    ## Rename S to the analyte name so the wide join works.
    names(out)[names(out) == "S"] <- nm
    out
  }

  wide_list <- lapply(analytes, make_df)

  ## Sequential full outer join on date.
  S_wide <- Reduce(
    function(a, b) merge(a, b, by = "date", all = TRUE),
    wide_list
  )

  ## ── 4. Pairwise Pearson correlation ──────────────────────────────────────

  ## use = "pairwise.complete.obs": each pair uses its own set of co-measured
  ## dates, critical for the ragged ARA panel where complete-case n can be 0.
  ## suppressWarnings: cor() warns "the standard deviation is zero" for
  ## degenerate analytes (constant S); we handle the resulting NAs explicitly.
  S_mat <- as.matrix(S_wide[, analytes, drop = FALSE])
  R_hat <- suppressWarnings(
    stats::cor(S_mat, use = "pairwise.complete.obs")
  )

  ## Replace any NA off-diagonal produced by zero co-observed pairs with 0
  ## (forcing independence for that pair — conservative).
  R_hat[is.na(R_hat)] <- 0

  ## Restore the diagonal to exactly 1 (it should be 1 already, but NA masking
  ## above could have zeroed it if a single-analyte SD is zero at this stage).
  diag(R_hat) <- 1

  ## ── 5. Force degenerate analytes to identity row/col ─────────────────────

  ## A degenerate analyte has no meaningful covariance with any partner;
  ## setting its row and column to 0 (and diagonal to 1) declares independence.
  for (nm in analytes[degen]) {
    R_hat[nm, ]  <- 0
    R_hat[, nm]  <- 0
    R_hat[nm, nm] <- 1
  }

  ## ── 6. Ridge shrinkage ───────────────────────────────────────────────────

  ## λ = 0.1 is a fixed stability hyperparameter: pulls all off-diagonals
  ## toward 0 by 10 %, guarding against spurious large correlations from
  ## sparse co-measurement overlap.  Conservative choice preferred over
  ## data-driven shrinkage (Schäfer–Strimmer) because that requires a
  ## complete-case matrix, which may not exist for the ragged ARA panel.
  lambda  <- 0.1
  R_ridge <- (1 - lambda) * R_hat + lambda * diag(p)

  ## Ensure dimnames are preserved.
  dimnames(R_ridge) <- list(analytes, analytes)

  ## ── 7. Positive-definiteness guarantee ───────────────────────────────────

  ## Pairwise estimates with unequal n per pair can be indefinite even after
  ## ridge shrinkage.  nearPD (Higham 2002) finds the nearest PD correlation
  ## matrix in the Frobenius norm.  Only invoke when strictly necessary to
  ## avoid unnecessary rounding of a matrix that is already PD.
  min_eig <- min(eigen(R_ridge, symmetric = TRUE, only.values = TRUE)$values)

  R_pd <- if (min_eig <= 1e-8) {
    as.matrix(Matrix::nearPD(R_ridge, corr = TRUE, keepDiag = TRUE)$mat)
  } else {
    R_ridge
  }

  ## Strip dimnames: downstream consumers use integer indexing; named diag()
  ## would cause expect_equal(diag(R), c(1,1)) to fail due to name mismatch.
  dimnames(R_pd) <- NULL

  ## ── 8. Return ────────────────────────────────────────────────────────────

  list(
    R        = R_pd,
    analytes = analytes,
    lambda   = lambda
  )
}
