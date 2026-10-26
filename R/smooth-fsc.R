## ----------------------------------------------------------------------------
## "fsc": factor smooth with a canonical null-space basis
## ----------------------------------------------------------------------------

#' Factor smooth with a canonical null-space basis (`bs = "fsc"`)
#'
#' A drop-in variant of \pkg{mgcv}'s factor-smooth interaction
#' (`bs = "fs"`, see [mgcv::factor.smooth]) whose model matrix is stable under
#' tiny (ulp-level) changes in the covariate. Use it as
#' `s(x, fac, bs = "fsc", k = ...)`; arguments `xt`, `m` and `k` behave as for
#' `"fs"`.
#'
#' The `"fs"` constructor reparameterises the base smooth with
#' `mgcv:::nat.param()` and gives each of the null-space columns its own ridge
#' penalty. For the usual second-order base smooth the null space is
#' two-dimensional (constant and linear) and the two basis vectors that
#' `nat.param()` returns for it span a degenerate eigenspace, so their
#' orientation within that span is arbitrary and can rotate under a 1e-12
#' relative change in `x`. Because each column carries its own penalty, the
#' rotation changes the fitted model, which makes criteria such as AIC
#' discontinuous in the data.
#'
#' `"fsc"` replaces those two columns with
#' `r * [1, (x - mean(x)) / sd(x)]` -- the same span in a fixed orientation,
#' with `r` the root-mean-square size of the original null-space columns -- and
#' updates the reparameterisation matrix `P` to match, so the `"fs"`
#' prediction method ([mgcv::Predict.matrix()] for class `"fs.interaction"`)
#' reproduces the fit on new data. Column space, penalty structure and number
#' of smoothing parameters are those of `"fs"`. When the null space of the
#' base smooth is not two-dimensional, the smooth is not one-dimensional, or
#' the canonical columns are not in the null-space span, the result is plain
#' `"fs"`.
#'
#' @param object A smooth specification object of class `"fsc.smooth.spec"`,
#'   as produced by `s(..., bs = "fsc")`.
#' @param data A list or data frame containing the covariate and the factor.
#' @param knots Optional list of knots, passed to the base smooth.
#'
#' @return A smooth object of class `"fs.interaction"`, as returned by the
#'   `"fs"` constructor.
#'
#' @examples
#' d <- data.frame(x = runif(90, 40, 80), f = factor(rep(c("a", "b", "c"), 30)))
#' d$y <- sin(d$x / 10) + as.numeric(d$f) + rnorm(90, 0, 0.2)
#' fit <- mgcv::gam(y ~ s(x, f, bs = "fsc", k = 5), data = d, method = "REML")
#'
#' @exportS3Method mgcv::smooth.construct
#' @keywords internal
smooth.construct.fsc.smooth.spec <- function(object, data, knots) {
  class(object) <- "fs.smooth.spec"
  ob <- mgcv::smooth.construct(object, data, knots)

  ## Only the factor-by-metric case built by fs (not the re / plain fallbacks).
  if (!inherits(ob, "fs.interaction") || is.null(ob$base) ||
      length(ob$base$term) != 1L) {
    return(ob)
  }
  p  <- ncol(ob$P)
  rk <- ob$base$rank
  if (p - rk != 2L) return(ob)

  x <- data[[ob$base$term]]
  if (!is.numeric(x)) return(ob)
  sx <- stats::sd(x)
  if (!is.finite(sx) || sx <= 0) return(ob)

  null_ind <- (rk + 1L):p
  gamm <- is.null(ob$Xb) # the gamm path keeps the per-level basis in X
  Xb <- if (gamm) ob$X else ob$Xb
  N <- Xb[, null_ind, drop = FALSE]

  ## Canonical orientation of the same span.  The Frobenius norm of N is
  ## invariant to its arbitrary rotation, so r is stable too.
  r  <- sqrt(sum(N^2) / (2 * nrow(N)))
  Tm <- cbind(1, (x - mean(x)) / sx) * r
  qrN <- qr(N)
  if (qrN$rank < 2L) return(ob)
  Q <- qr.coef(qrN, Tm)
  ## Not the same span (e.g. an exotic base smooth): keep plain fs.
  if (max(abs(N %*% Q - Tm)) > 1e-8 * max(1, max(abs(Tm)))) return(ob)

  ## Keep X and P consistent (X = basis %*% P), so prediction on the
  ## training data reproduces the fit; N %*% Q equals Tm to rounding.
  ob$P[, null_ind] <- ob$P[, null_ind, drop = FALSE] %*% Q
  Xb[, null_ind] <- N %*% Q

  if (gamm) {
    ob$X <- Xb
  } else {
    ob$Xb <- Xb
    nf <- length(ob$flev)
    X <- matrix(0, nrow(Xb), p * nf)
    ind <- seq_len(p)
    for (i in seq_len(nf)) {
      X[, ind] <- Xb * as.numeric(ob$fac == ob$flev[i])
      ind <- ind + p
    }
    ob$X <- X
  }
  ## Discretised (bam) fits carry a copy of the base basis as a margin.
  if (!is.null(ob$margin) && length(ob$margin) >= 2L) {
    ob$margin[[2L]]$P <- ob$P
    ob$margin[[2L]]$X <- Xb
  }
  ob
}
