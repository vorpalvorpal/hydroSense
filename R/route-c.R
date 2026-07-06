## ============================================================================
## Route C - low-rank censored factor model (dev/plan-route-c.md)
## ============================================================================
##
## Stan-free conditional-prediction kernels for the censored factor imputation
## model: y = mu + Lambda %*% f + eps, f ~ N(0, I_k), eps ~ N(0, diag(psi)).
## Given a partially-observed sample (some analytes observed, some missing),
## these compute the Gaussian conditional of the missing block and draw from
## it, truncating any BDL target at its detection limit. See
## `.rc_ref_conditional()` in tests/testthat/helper-route-c.R for the
## brute-force full-Sigma reference these must match.

#' Gaussian conditional of missing analytes given observed ones under a
#' low-rank factor covariance
#'
#' Uses the k x k Woodbury form (posterior over the latent factor `f`) rather
#' than inverting the full J x J residual covariance `Sigma = Lambda %*% t(Lambda)
#' + diag(psi)` — the whole reason for the low-rank structure. Stable even
#' from a single observed analyte.
#'
#' @param y Named or positional length-J vector: observed residual (or raw
#'   value, on whatever scale `mu`/`Lambda`/`psi` share) for observed analytes,
#'   `NA` for missing/BDL ones.
#' @param mu Length-J mean vector.
#' @param Lambda J x k loading matrix.
#' @param psi Length-J idiosyncratic variances (not sd).
#' @return `list(mean, cov)` for the missing block, in the order `which(is.na(y))`.
#' @keywords internal
.factor_condition <- function(y, mu, Lambda, psi) {
  Lambda <- as.matrix(Lambda)
  k <- ncol(Lambda)
  miss <- which(is.na(y))
  obs  <- which(!is.na(y))

  Lambda_M <- Lambda[miss, , drop = FALSE]

  if (length(obs) == 0L) {
    cov_M <- Lambda_M %*% t(Lambda_M) + diag(psi[miss], length(miss))
    return(list(mean = mu[miss], cov = cov_M))
  }

  Lambda_O   <- Lambda[obs, , drop = FALSE]
  psi_O_inv  <- 1 / psi[obs]
  r_O        <- (y - mu)[obs]

  # V = (I_k + Lambda_O' Psi_O^-1 Lambda_O)^-1 ; m = V Lambda_O' Psi_O^-1 r_O
  M <- diag(k) + t(Lambda_O) %*% (Lambda_O * psi_O_inv)
  V <- solve(M)
  m <- V %*% (t(Lambda_O) %*% (r_O * psi_O_inv))

  mean_M <- mu[miss] + as.numeric(Lambda_M %*% m)
  cov_M  <- Lambda_M %*% V %*% t(Lambda_M) + diag(psi[miss], length(miss))

  list(mean = mean_M, cov = cov_M)
}

#' Draws from the conditional missing-analyte distribution, truncating BDL
#' targets at their detection limit
#'
#' @param upper Optional length-`length(which(is.na(y)))` vector aligned with
#'   the missing block (same order `.factor_condition()` returns): a finite
#'   value truncates that column's draws to `(-inf, upper]` (a BDL target's
#'   `log(DL)`); `NA` leaves the column unbounded.
#' @inheritParams .factor_condition
#' @param ndraws Number of draws.
#' @return `ndraws` x `length(missing)` matrix.
#' @keywords internal
.factor_condition_draw <- function(y, mu, Lambda, psi, ndraws, upper = NULL) {
  cond  <- .factor_condition(y, mu, Lambda, psi)
  cmean <- as.numeric(cond$mean)
  ccov  <- as.matrix(cond$cov)
  m     <- length(cmean)

  if (m == 0L) return(matrix(numeric(0), nrow = ndraws, ncol = 0))

  ccov <- (ccov + t(ccov)) / 2   # guard against asymmetric FP noise before chol()
  Lc   <- chol(ccov)            # Lc' Lc == ccov
  z    <- matrix(stats::rnorm(ndraws * m), nrow = ndraws, ncol = m)
  draws <- z %*% Lc + matrix(cmean, nrow = ndraws, ncol = m, byrow = TRUE)

  if (!is.null(upper) && any(!is.na(upper))) {
    if (!requireNamespace("truncnorm", quietly = TRUE)) {
      cli::cli_abort(c(
        "Truncating BDL target draws needs the {.pkg truncnorm} package.",
        "i" = "Install with {.code install.packages(\"truncnorm\")}."
      ))
    }
    sds <- sqrt(diag(ccov))
    for (j in seq_len(m)) {
      if (!is.na(upper[j])) {
        draws[, j] <- truncnorm::rtruncnorm(
          ndraws, a = -Inf, b = upper[j], mean = cmean[j], sd = sds[j]
        )
      }
    }
  }

  dimnames(draws) <- NULL
  draws
}
