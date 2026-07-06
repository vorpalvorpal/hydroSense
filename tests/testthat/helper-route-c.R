## ============================================================================
## Shared builders for the Route C BDD suite (test-route-c.R)
## ============================================================================
##
## Route C is the low-rank censored factor model that replaces the brms
## rescor_mi / cens / cens_factor imputation path (see dev/plan-route-c.md).
## These specs are written test-first: they describe the *target* behaviour and
## call functions that do not exist yet, so each is guarded with
## `.skip_route_c()` and skips until the implementation lands. Drive them
## red->green by deleting the guard as each piece is built.
##
## The prediction-core specs (the finding-3 conditioning) are deliberately
## Stan-free: they inject a KNOWN factor model (Lambda, Psi, mu) and exercise the
## pure linear-algebra kernels, so the borrowing-strength behaviour can be
## verified without a Stan toolchain. The fit and end-to-end specs need
## cmdstanr and additionally skip when it is absent.

## Mark a pending Route C spec. `what` is a short description of the missing
## piece; the message tells the implementer which seam to build.
.skip_route_c <- function(what) {
  testthat::skip(paste0(
    "Route C not implemented — ", what,
    " (see dev/plan-route-c.md; delete .skip_route_c() to activate)"
  ))
}

## A known low-rank factor model for the Stan-free prediction kernels.
##
## Returns the ground-truth pieces the conditional-prediction kernels consume:
##   - analytes : J original analyte names (order defines the row order of y/mu/Lambda)
##   - Lambda   : J x k loadings (known)
##   - psi      : length-J idiosyncratic VARIANCES (not sd)
##   - mu       : length-J per-analyte mean for a single sample (mu_j(X_i))
## Sigma = Lambda %*% t(Lambda) + diag(psi) is the induced residual covariance.
.rc_factor <- function(J = 4L, k = 2L, seed = 1L) {
  set.seed(seed)
  analytes <- c("Cu", "Zn", "Cd", "Pb", "Ni", "As")[seq_len(J)]
  Lambda <- matrix(stats::rnorm(J * k, 0, 0.8), nrow = J, ncol = k,
                   dimnames = list(analytes, paste0("F", seq_len(k))))
  psi <- stats::setNames(stats::runif(J, 0.05, 0.15), analytes)
  mu  <- stats::setNames(stats::rnorm(J, 0, 1), analytes)
  list(analytes = analytes, Lambda = Lambda, psi = psi, mu = mu, k = k)
}

## Brute-force reference for the Gaussian conditional of the missing block given
## the observed block, computed directly from the full Sigma = LL' + diag(psi).
## The implementation kernel must match this (it may use the cheaper k x k
## Woodbury form, but the ANSWER is this).
.rc_ref_conditional <- function(y, mu, Lambda, psi) {
  Sigma <- Lambda %*% t(Lambda) + diag(psi)
  miss  <- which(is.na(y))
  obs   <- which(!is.na(y))
  if (length(obs) == 0L) {
    return(list(mean = mu[miss], cov = Sigma[miss, miss, drop = FALSE]))
  }
  r_o   <- (y - mu)[obs]
  S_mo  <- Sigma[miss, obs, drop = FALSE]
  S_oo  <- Sigma[obs, obs, drop = FALSE]
  S_mm  <- Sigma[miss, miss, drop = FALSE]
  cmean <- mu[miss] + as.numeric(S_mo %*% solve(S_oo, r_o))
  ccov  <- S_mm - S_mo %*% solve(S_oo, t(S_mo))
  list(mean = cmean, cov = ccov)
}
