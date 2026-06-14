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
##   .anchor_residual_cor     pairwise-complete Pearson R, ridge-shrunk + nearPD
##   .coupled_residual_draws  correlated DK draws via Cholesky of the empirical R

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


## ── Stage 3: correlated residual draws via the DK identity ───────────────────

#' Draw cross-analyte correlated Kalman residual trajectories
#'
#' Generates `ndraws` draws from each analyte's Durbin–Koopman simulation
#' smoother, with innovations correlated **across analytes** according to the
#' `[p × p]` correlation sub-matrix `cor_R`.  Each analyte's per-analyte
#' marginal distribution is unchanged by the correlation (DK identity);
#' only the *joint* distribution across analytes is affected.
#'
#' ## Statistical justification
#'
#' **DK identity** (Durbin & Koopman 2002, §4): the draw
#' \deqn{\tilde{\alpha} = \hat{\alpha} + (\alpha^+ - L y^+)}
#' guarantees each analyte's marginal distribution is unchanged when the
#' process innovations \eqn{\eta^+} are correlated: the smoother gain `L` is
#' precomputed on the actual data, and the correction \eqn{\alpha^+ - L y^+}
#' has expectation zero.  Correlating \eqn{\eta^+} across analytes widens only
#' the *joint* distribution (i.e. the combined AmsPAF interval).
#'
#' **Cholesky coupling**: `Z_raw %*% U_chol` (where `U_chol = chol(R)`, upper
#' triangular) gives rows with covariance `R` when rows of `Z_raw` are iid
#' \eqn{N(0, I_p)}.  The same `U_chol` is applied to both the process
#' innovations and the initial-state draws so that the coupling is coherent
#' along the entire trajectory.
#'
#' **Independent observation noise**: `eps_std` is drawn independently per
#' analyte because measurement errors from separate grab analyses are
#' uncorrelated even when the true concentrations co-move.
#'
#' **Ragged grids**: `match(grid_a, U)` slices the union-grid field to each
#' analyte's own clipped grid.  Days where an analyte has no entry (outside its
#' grab span) have no field element — no coupling on those days, which is
#' correct.
#'
#' @param smoothers Named list by analyte.  Each element is a list with
#'   `grid_dates` (Date vector), `mean` (numeric vector, KFS posterior mean),
#'   and `draw_model` (KFAS model from [.build_kalman_model()], or `NULL` for
#'   a degenerate / non-modelled analyte).
#' @param modelled Character vector of analyte names to include in the output
#'   (subset of `names(smoothers)`).
#' @param ndraws Positive integer; number of draws to generate.
#' @param cor_R `[p × p]` numeric positive-definite correlation matrix (no
#'   dimnames required).  Typically the `$R` component from
#'   [.anchor_residual_cor()].
#' @param cor_analytes Character vector giving the analyte order for the rows
#'   and columns of `cor_R`.  Defaults to `rownames(cor_R)` when `cor_R` has
#'   dimnames; must be supplied explicitly otherwise.
#' @param seed Integer or `NULL`.  When non-`NULL`, `set.seed(seed)` is called
#'   at the start of the function so results are reproducible.  In
#'   `amspaf_daily()`, the seed is set by the caller (via its own `seed`
#'   argument) before this function is invoked, so `seed = NULL` is used there.
#'
#' @return Named list by analyte (same order as `modelled`).  Each element is
#'   `list(grid_dates = <Date>, draws = <matrix [n_grid × ndraws]>)`.
#'   Analytes with `draw_model = NULL` receive flat draws equal to their KFS
#'   posterior mean (the degenerate path).  Analytes with empty `grid_dates`
#'   receive `draws = NULL`.
#'
#' @references
#'   Durbin J, Koopman SJ (2002). A simple and efficient simulation smoother
#'   for state space time series analysis. *Biometrika* **89**(3), 603–616.
#'   <https://doi.org/10.1093/biomet/89.3.603>
#'
#' @keywords internal
.coupled_residual_draws <- function(smoothers, modelled, ndraws, cor_R,
                                    cor_analytes = rownames(cor_R),
                                    seed         = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))

  ## ── 1. Identify couplable analytes ─────────────────────────────────────────

  ## Couplable = in modelled AND has a non-NULL draw_model.
  ## Analytes with draw_model = NULL use the degenerate flat path below.
  couplable_nms <- names(Filter(
    function(s) !is.null(s$draw_model),
    smoothers[modelled]
  ))
  p <- length(couplable_nms)

  ## ── 2. Extract the correlation sub-matrix for couplable analytes ────────────

  idx       <- match(couplable_nms, cor_analytes)
  cor_R_sub <- cor_R[idx, idx, drop = FALSE]

  ## ── 3. Cholesky factor of cor_R_sub ────────────────────────────────────────

  ## Upper triangular U s.t. t(U) %*% U = cor_R_sub.
  ## Z_raw %*% U_chol gives rows with covariance cor_R_sub (Cholesky coupling).
  U_chol <- chol(cor_R_sub)

  ## ── 4. Union daily grid ─────────────────────────────────────────────────────

  ## Build the union of all analytes' grid dates; used to index the correlated
  ## field. Each analyte's own grid is a subset of this union (ragged grids).
  U   <- sort(unique(do.call(c, lapply(modelled, function(nm) {
    smoothers[[nm]]$grid_dates
  }))))
  n_U <- length(U)

  ## ── 5. Precompute DK setup for each couplable analyte ──────────────────────

  ## .kalman_sim_smoother_setup() extracts the gain matrix L, KFS mean x_hat,
  ## and noise scaling vectors needed by .kalman_draw_coupled().
  setups <- stats::setNames(
    lapply(couplable_nms,
           function(nm) .kalman_sim_smoother_setup(smoothers[[nm]]$draw_model)),
    couplable_nms
  )

  ## ── 6. Draw the correlated field ────────────────────────────────────────────

  ## Correlated process innovations: [n_U * ndraws rows × p analytes].
  ## Each block of n_U rows (one per draw) has row-covariance cor_R_sub.
  ## We draw the full union-grid × ndraws block at once and slice per-analyte
  ## below using grid_a_idx.
  n_rows  <- n_U * ndraws
  Z_raw   <- matrix(stats::rnorm(n_rows * p), nrow = n_rows, ncol = p)
  Z_cor   <- Z_raw %*% U_chol   # [n_rows × p], rows have cov = cor_R_sub

  ## Correlated initial states (one per draw, all p analytes at once).
  Z_a1_raw <- matrix(stats::rnorm(ndraws * p), nrow = ndraws, ncol = p)
  Z_a1c    <- Z_a1_raw %*% U_chol   # [ndraws × p]

  ## Observation noise is kept INDEPENDENT across analytes: measurement errors
  ## from separate grab analyses are uncorrelated even when true concentrations
  ## co-move (see §3 of the plan).  Drawn per-analyte inside the loop below.

  ## ── 7. Per-analyte draw ─────────────────────────────────────────────────────

  res_draws <- stats::setNames(lapply(modelled, function(nm) {
    sm       <- smoothers[[nm]]
    grid_a   <- sm$grid_dates

    ## Degenerate / empty path.
    if (is.null(sm) || length(grid_a) == 0L) {
      return(list(grid_dates = as.Date(character()), draws = NULL))
    }

    ## No draw_model: flat draws equal to the KFS posterior mean.
    if (is.null(sm$draw_model)) {
      return(list(
        grid_dates = grid_a,
        draws      = matrix(sm$mean, nrow = length(grid_a), ncol = ndraws)
      ))
    }

    ## Should not happen (couplable_nms = modelled with non-NULL draw_model),
    ## but guard: fall back to independent .kalman_draw.
    if (!nm %in% couplable_nms) {
      return(list(
        grid_dates = grid_a,
        draws      = .kalman_draw(sm$draw_model, ndraws)
      ))
    }

    a_col      <- match(nm, couplable_nms)   # column in Z_cor / Z_a1c
    setup      <- setups[[nm]]
    n_grid_a   <- length(grid_a)
    m_a        <- length(setup$pos)

    ## Indices of this analyte's grid dates in the union grid U.
    ## Ragged grids: only union-grid rows that fall within analyte A's span
    ## participate in the coupling field; days outside the span have no entry.
    grid_a_idx <- match(grid_a, U)   # integer indices into U (length n_grid_a)

    ## Row indices into Z_cor for each (draw, day) pair.
    ## Z_cor row (k-1)*n_U + grid_a_idx[t] gives draw k, analyte-grid day t.
    row_idx <- rep(grid_a_idx, ndraws) +
               n_U * rep(seq_len(ndraws) - 1L, each = n_grid_a)

    ## Extract this analyte's correlated process innovations: [n_grid_a × ndraws].
    eta_std_a <- matrix(Z_cor[row_idx, a_col], nrow = n_grid_a, ncol = ndraws)

    ## Correlated initial state for this analyte (one value per draw).
    a1_z_a <- Z_a1c[, a_col]

    ## Independent observation noise at the anchor positions.
    eps_std_a <- matrix(stats::rnorm(m_a * ndraws), nrow = m_a, ncol = ndraws)

    dr <- .kalman_draw_coupled(setup, eta_std_a, a1_z_a, eps_std_a)
    list(grid_dates = grid_a, draws = dr)
  }), modelled)

  res_draws
}
