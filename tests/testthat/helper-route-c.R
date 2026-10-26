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

## Skip a spec that actually FITS a factor model unless a working CmdStan is
## available. cmdstanr the package can be installed while CmdStan (the backend
## binary) is not set up (as on CI), so `skip_if_not_installed("cmdstanr")` is
## not enough. Mirrors the package's own `.brms_backend()` availability probe.
.skip_if_no_cmdstan <- function() {
  testthat::skip_if_not_installed("cmdstanr")
  ok <- !is.null(tryCatch(suppressMessages(cmdstanr::cmdstan_version()),
                          error = function(e) NULL))
  if (!ok) testthat::skip("CmdStan not installed / path not set")
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

## Reference convergence metric for a Route C factor fit: the max Rhat over the
## unique entries of the IDENTIFIED implied covariance
## Sigma = Lambda Lambda' + diag(psi). `draws` is a posterior draws object
## holding `Lambda[j,k]` and `psi[j]`. In low-df groups Lambda Lambda' and psi
## trade off against each other (neither is identified on its own), so only
## the full Sigma is a meaningful convergence target.
.rc_full_sigma_rhat <- function(draws) {
  d  <- posterior::as_draws_array(draws)
  vn <- dimnames(d)[[3]]
  lam <- grep("^Lambda\\[", vn, value = TRUE)
  ij  <- do.call(rbind, lapply(strsplit(gsub("Lambda\\[|\\]", "", lam), ","),
                               as.integer))
  J <- max(ij[, 1]); K <- max(ij[, 2])
  r <- c()
  for (a in seq_len(J)) for (b in a:J) {
    acc <- 0
    for (k in seq_len(K)) {
      acc <- acc + d[, , sprintf("Lambda[%d,%d]", a, k)] *
                   d[, , sprintf("Lambda[%d,%d]", b, k)]
    }
    if (a == b) acc <- acc + d[, , sprintf("psi[%d]", a)]
    r <- c(r, posterior::rhat(matrix(acc, nrow = dim(d)[1])))
  }
  max(r)
}

## A stand-in for a fitted Route C Stan object, built from a draws array, so
## .route_c_convergence() can be unit-tested without Stan. It answers the
## CmdStanMCMC calls the convergence check makes -- `$draws(variables)` and
## `$summary(variables)` -- from the supplied draws.
.rc_fake_fit <- function(draws) {
  draws <- posterior::as_draws_array(draws)
  pick <- function(variables) {
    if (is.null(variables)) return(draws)
    posterior::subset_draws(draws, variable = variables)
  }
  list(
    draws = function(variables = NULL, inc_warmup = FALSE,
                     format = "draws_array", ...) {
      posterior::as_draws(pick(variables), format = format)
    },
    summary = function(variables = NULL, ...) {
      posterior::summarise_draws(pick(variables))
    }
  )
}

## Draws for a J = 2, k = 1 factor fit where Lambda Lambda' DISAGREES across
## chains but the identified Sigma = Lambda Lambda' + Psi AGREES: chains 1-2
## put the shared variance on analyte 1's loading, chains 3-4 on analyte 2's,
## with psi taking up the difference. The off-diagonal Lambda_1 * Lambda_2 and
## each diagonal Lambda_j^2 + psi_j have one common distribution in every
## chain; the diagonal Lambda_j^2 differs between chain pairs by ~16x.
.rc_tradeoff_draws <- function(n_iter = 250L, seed = 1L) {
  withr::local_seed(seed)
  n_chain <- 4L
  arr <- array(NA_real_, dim = c(n_iter, n_chain, 5L),
               dimnames = list(NULL, NULL,
                               c("Lambda[1,1]", "Lambda[2,1]",
                                 "psi[1]", "psi[2]", "lp__")))
  for (ch in seq_len(n_chain)) {
    big   <- 2 * (1 + stats::rnorm(n_iter, 0, 0.02))
    prod  <- 1 + stats::rnorm(n_iter, 0, 0.05)        # Sigma[1,2]
    s1    <- 5 + stats::rnorm(n_iter, 0, 0.1)         # Sigma[1,1]
    s2    <- 5 + stats::rnorm(n_iter, 0, 0.1)         # Sigma[2,2]
    if (ch <= 2L) { l1 <- big; l2 <- prod / big } else { l2 <- big; l1 <- prod / big }
    arr[, ch, "Lambda[1,1]"] <- l1
    arr[, ch, "Lambda[2,1]"] <- l2
    arr[, ch, "psi[1]"]      <- s1 - l1^2
    arr[, ch, "psi[2]"]      <- s2 - l2^2
    arr[, ch, "lp__"]        <- stats::rnorm(n_iter, -100, 1)
  }
  posterior::as_draws_array(arr)
}

## Draws for a J = 2, k = 1 factor fit whose identified Sigma DISAGREES between
## interleaved chains: chains 1 and 3 centre Sigma[1,1] near 5, chains 2 and 4
## near 6.3 (via Lambda[1,1]). Each chain mixes well on its own, so only a
## between-CHAIN comparison catches it -- splitting the chains-concatenated
## sequence into halves (chains 1-2 vs 3-4) puts both modes in each half and
## reports ~1.0.
.rc_interleaved_draws <- function(n_iter = 250L, seed = 2L) {
  withr::local_seed(seed)
  n_chain <- 4L
  arr <- array(NA_real_, dim = c(n_iter, n_chain, 5L),
               dimnames = list(NULL, NULL,
                               c("Lambda[1,1]", "Lambda[2,1]",
                                 "psi[1]", "psi[2]", "lp__")))
  for (ch in seq_len(n_chain)) {
    centre <- if (ch %% 2L == 1L) 2 else 2.3
    arr[, ch, "Lambda[1,1]"] <- stats::rnorm(n_iter, centre, 0.02)
    arr[, ch, "Lambda[2,1]"] <- stats::rnorm(n_iter, 0.5, 0.02)
    arr[, ch, "psi[1]"]      <- stats::rnorm(n_iter, 1, 0.02)
    arr[, ch, "psi[2]"]      <- stats::rnorm(n_iter, 1, 0.02)
    arr[, ch, "lp__"]        <- stats::rnorm(n_iter, -100, 1)
  }
  posterior::as_draws_array(arr)
}
