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
#' + diag(psi)` -- the whole reason for the low-rank structure. Stable even
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


## ----------------------------------------------------------------------------
## Stage-2 Stan model: compile-once cache
## ----------------------------------------------------------------------------

.route_c_model_env <- new.env(parent = emptyenv())

#' Stop with a friendly message if cmdstanr is unavailable
#'
#' The Route C `"factor"` imputation method uses \pkg{cmdstanr} + \pkg{mgcv}
#' (it does **not** use \pkg{brms}), so this is the method-appropriate engine
#' guard -- the analogue of [.require_brms()] for the three brms-based methods.
#' @keywords internal
.require_cmdstanr <- function() {
  if (requireNamespace("cmdstanr", quietly = TRUE)) return(invisible(TRUE))
  cli::cli_abort(c(
    "The Route C {.val factor} imputation method needs the {.pkg cmdstanr} \\
     package (plus a working CmdStan install).",
    "i" = "Install with {.code install.packages(\"cmdstanr\", repos = \\
           \"https://stan-dev.r-universe.dev\")} then {.code \\
           cmdstanr::install_cmdstan()}."
  ))
}

#' Compiled Stage-2 censored-factor Stan model (cached for the session)
#' @keywords internal
.route_c_stan_model <- function() {
  .require_cmdstanr()
  if (is.null(.route_c_model_env$mod)) {
    stan_file <- system.file("stan", "factor_censored.stan", package = "hydroSense")
    if (!nzchar(stan_file)) {
      cli::cli_abort("Could not locate {.file inst/stan/factor_censored.stan} \\
                       in the installed package.")
    }
    .route_c_model_env$mod <- cmdstanr::cmdstan_model(stan_file)
  }
  .route_c_model_env$mod
}


#' Rhat for a Route C Stage-2 Stan fit
#'
#' `.fit_group_model_factor()` tags its Stage-2 `CmdStanMCMC` fit with class
#' `route_c_stanfit` so generic diagnostics like `brms::rhat()` (itself just
#' `posterior::rhat()`, a plain S3 generic) work on it directly. The raw
#' cmdstanr object has no `rhat()` method of its own to dispatch to; this
#' delegates to its `$summary()` method, which already computes per-parameter
#' Rhat correctly (unlike calling `posterior::rhat()` on the whole multi-
#' parameter draws array, which silently blends every parameter into one
#' meaningless number).
#' @param x A `route_c_stanfit` object (a `CmdStanMCMC` with this class prepended).
#' @param ... Unused; present for S3 signature compatibility.
#' @return Named numeric vector of per-parameter Rhat.
#' @keywords internal
#' @exportS3Method posterior::rhat
rhat.route_c_stanfit <- function(x, ...) {
  s <- x$summary()
  stats::setNames(s$rhat, s$variable)
}

#' Rotation-invariant convergence diagnostics for a Route C factor fit
#'
#' A low-rank factor model is only identified up to rotation of `Lambda`
#' (`Lambda R` gives the same `Lambda Lambda'`). With the positive-lower-
#' triangular constraint this is mostly pinned, but at small `N` / weak signal
#' the sampler can still explore near-rotationally-equivalent modes, inflating
#' the per-element Rhat of `Lambda` **without any consequence for imputation** --
#' prediction uses only `Sigma = Lambda Lambda' + Psi`, which is rotation-
#' invariant (see [.factor_condition()]). So the meaningful convergence check
#' monitors the invariant functionals: the implied residual covariance `Sigma`,
#' the idiosyncratic variances `psi`, and the log-density `lp__` -- not raw
#' `Lambda`. Empirically `Sigma` converges cleanly (Rhat ~ 1.0) even when a
#' `Lambda` element trips a naive `max(rhat(fit))` gate.
#'
#' @param group A fitted `factor`-method group (from [.fit_group_model_factor()]).
#' @return A list with `sigma_rhat`, `psi_rhat`, `lp_rhat` (max Rhat over the
#'   `Sigma` entries, over `psi`, and for `lp__`); all `NA` for a degenerate
#'   (single-analyte, Stan-free) group.
#' @keywords internal
.route_c_convergence <- function(group) {
  na <- list(sigma_rhat = NA_real_, psi_rhat = NA_real_, lp_rhat = NA_real_)
  if (isTRUE(group$degenerate) || is.null(group$fit)) return(na)

  fit <- group$fit
  J   <- length(group$analytes)
  k   <- group$k

  s <- fit$summary(c("psi", "lp__"))
  psi_rhat <- max(s$rhat[grepl("^psi\\[", s$variable)], na.rm = TRUE)
  lp_rhat  <- max(s$rhat[s$variable == "lp__"], na.rm = TRUE)

  # Implied Sigma = Lambda Lambda' (rotation-invariant), Rhat per unique entry.
  ld   <- posterior::as_draws_array(fit$draws("Lambda"))  # iter x chain x param
  vn   <- dimnames(ld)[[3]]
  ij   <- do.call(rbind, lapply(vn, function(v)
    as.integer(strsplit(gsub("Lambda\\[|\\]", "", v), ",")[[1]])))
  pairs <- which(upper.tri(matrix(0, J, J), diag = TRUE), arr.ind = TRUE)
  sigma_rhats <- vapply(seq_len(nrow(pairs)), function(p) {
    a <- pairs[p, 1]; b <- pairs[p, 2]
    acc <- 0
    for (kk in seq_len(k)) {
      ia <- which(ij[, 1] == a & ij[, 2] == kk)
      ib <- which(ij[, 1] == b & ij[, 2] == kk)
      acc <- acc + ld[, , ia] * ld[, , ib]
    }
    posterior::rhat(acc)   # acc is iterations x chains
  }, numeric(1))

  list(sigma_rhat = max(sigma_rhats, na.rm = TRUE),
       psi_rhat = psi_rhat, lp_rhat = lp_rhat)
}


## ----------------------------------------------------------------------------
## Fit -- the two-stage censored factor model
## ----------------------------------------------------------------------------

#' Fit the Route C two-stage censored factor group model
#'
#' Stage 1: per-analyte `mgcv::gam(lv ~ s(PC1) + ...)` on detected
#' observations gives the spline mean `mu_j(X_i)`; detected residuals and BDL
#' censoring bounds (`log(DL) - mu_j`) are computed from it. Stage 2: a Stan
#' program (`inst/stan/factor_censored.stan`) fits a rank-`k` factor model on
#' those residuals, with BDL cells as upper-bounded latent parameters (proper
#' left-censoring, jointly with the factor and covariance).
#' @keywords internal
.fit_group_model_factor <- function(target_analytes, safe_analytes, base,
                                     pc_wide, pc_cols, log_floors,
                                     iter, warmup, chains, cores,
                                     group_name = "group", k = NULL,
                                     seed = NULL, ...) {
  J <- length(target_analytes)
  safe_vec <- unname(safe_analytes[target_analytes])

  # A single target analyte has no cross-analyte coupling to model at all: a
  # rank-K factor structure with K < J is structurally impossible (K must be
  # >= 1 and < J), so a Stan factor fit is neither necessary nor identifiable
  # here. Fall back to a Stage-1-only path: the GAM mean plus a plug-in
  # residual variance, with BDL/missing draws from that marginal Normal
  # (truncated at log(DL) for BDL cells). This must not error out -- a
  # single-analyte group is a realistic outcome of the default catch-all
  # group routing on real data.
  if (J == 1L) {
    if (!is.null(k)) {
      cli::cli_inform(c(
        "i" = "{.arg k} = {k} ignored for group {.val {group_name}}: a \\
               single-analyte group always uses the Stage-1-only marginal fit."
      ))
    }
    return(.fit_group_model_factor_degenerate(
      target_analytes = target_analytes, safe_analytes = safe_analytes,
      base = base, pc_wide = pc_wide, pc_cols = pc_cols,
      log_floors = log_floors, group_name = group_name,
      # Match a Stan group's posterior draw count so a mixed Stan+degenerate
      # model emits one consistent draw domain in return = "draws" mode.
      post_ndraws = chains * max(1L, iter - warmup)
    ))
  }

  if (is.null(k)) k <- min(2L, J - 1L)
  if (k < 1L) {
    cli::cli_abort(c(
      "Route C {.val factor} method needs at least 2 target analytes in \\
       group {.val {group_name}}; got {J}.",
      "i" = "Route it into a larger group or drop it from the model."
    ))
  }

  rhs <- paste(paste0("s(", pc_cols, ")"), collapse = " + ")

  # -- Stage 1: per-analyte GAM mean + residuals/censoring bounds ------------
  gams       <- stats::setNames(vector("list", J), target_analytes)
  resid_rows <- vector("list", J)

  for (i in seq_len(J)) {
    a   <- target_analytes[i]
    s   <- safe_vec[i]
    dat <- base[base$analyte == a, c("sample_id", "detected", "lv"), drop = FALSE]
    dat <- dplyr::left_join(dat, pc_wide, by = "sample_id")
    dat <- dat[stats::complete.cases(dat[, pc_cols, drop = FALSE]), , drop = FALSE]

    n_detected <- sum(dat$detected)
    if (n_detected < 4L) {
      cli::cli_abort(c(
        "Route C Stage-1 GAM for {.val {a}} needs at least 4 detected \\
         observations with PC scores; got {n_detected}."
      ))
    }

    gam_formula <- stats::as.formula(paste0("lv ~ ", rhs))
    gam_fit <- mgcv::gam(gam_formula, data = dat[dat$detected, , drop = FALSE],
                        family = stats::gaussian())
    gams[[a]] <- gam_fit

    dat$mu_pred <- as.numeric(stats::predict(gam_fit, newdata = dat, type = "response"))
    # Detected: genuine residual. BDL: `lv` already holds log(DL) (built from
    # the stored detection-limit `value`), so this is the censoring bound
    # b_ij = log(DL) - mu_j(X_i), not an observed residual.
    resid_rows[[i]] <- tibble::tibble(
      sample_id = dat$sample_id,
      safe      = s,
      detected  = dat$detected,
      resid     = dat$lv - dat$mu_pred
    )
  }

  resid_long <- dplyr::bind_rows(resid_rows)

  # -- Stage 2: censored factor model on the residuals -----------------------
  sample_levels <- sort(unique(resid_long$sample_id))
  N <- length(sample_levels)
  row_idx <- match(resid_long$sample_id, sample_levels)
  col_idx <- match(resid_long$safe, safe_vec)
  obs_mask <- resid_long$detected

  standata <- list(
    N = N, J = J, K = as.integer(k),
    N_obs = sum(obs_mask), N_cens = sum(!obs_mask),
    obs_row = row_idx[obs_mask], obs_col = col_idx[obs_mask],
    r_obs   = resid_long$resid[obs_mask],
    cen_row = row_idx[!obs_mask], cen_col = col_idx[!obs_mask],
    b_cens  = resid_long$resid[!obs_mask]
  )

  mod <- .route_c_stan_model()
  sample_args <- utils::modifyList(
    list(
      chains          = chains,
      parallel_chains = max(1L, min(chains, cores)),
      iter_warmup     = warmup,
      iter_sampling   = max(1L, iter - warmup),
      adapt_delta     = 0.95,
      seed            = seed,
      refresh         = 0
    ),
    list(...)
  )
  sample_args$data <- standata

  cli::cli_inform(c(
    "i" = "Stan {group_name} (factor, k = {k}): {J} analyte{?s} x {N} \\
           sample{?s}. This may take a few minutes."
  ))

  fit <- tryCatch(
    do.call(mod$sample, sample_args),
    error = function(e) cli::cli_abort(c(
      "Stan sampling failed for the Route C factor model.",
      "x" = "{conditionMessage(e)}"
    ))
  )
  class(fit) <- c("route_c_stanfit", class(fit))

  # -- Posterior-mean point summaries (diagnostics; draws live in `fit`) ------
  Lambda_point <- matrix(0, nrow = J, ncol = k,
                        dimnames = list(safe_vec, paste0("F", seq_len(k))))
  lam_summary <- fit$summary("Lambda")
  lam_ij <- regmatches(lam_summary$variable,
                      regexpr("(?<=\\[)[0-9]+,[0-9]+(?=\\])", lam_summary$variable, perl = TRUE))
  lam_ij <- do.call(rbind, strsplit(lam_ij, ","))
  Lambda_point[cbind(as.integer(lam_ij[, 1]), as.integer(lam_ij[, 2]))] <- lam_summary$mean

  psi_summary <- fit$summary("psi")
  psi_j <- as.integer(regmatches(psi_summary$variable,
                      regexpr("(?<=\\[)[0-9]+(?=\\])", psi_summary$variable, perl = TRUE)))
  Psi_point <- stats::setNames(numeric(J), safe_vec)
  Psi_point[psi_j] <- psi_summary$mean

  list(
    gams            = gams,
    Lambda          = Lambda_point,
    Psi             = Psi_point,
    k               = as.integer(k),
    fit             = fit,
    analytes        = target_analytes,
    safe_names      = safe_analytes,
    pc_cols         = pc_cols,
    wide_sample_ids = sample_levels,
    impute_method   = "factor",
    log_floors      = log_floors,
    post_ndraws     = chains * max(1L, iter - warmup)
  )
}

#' Stage-1-only fit for a single-analyte factor group (no cross-analyte
#' coupling is possible with J = 1)
#'
#' Fits the Stage-1 GAM as usual and estimates a plug-in residual variance
#' from the detected residuals. No Stan model is fit; `Lambda` is a J x 0
#' matrix (`k = 0`) so `.factor_condition()`'s existing `length(obs) == 0`
#' branch (always true here -- a single-analyte group's own target is either
#' fully observed, in which case there is nothing to predict, or missing, in
#' which case there is nothing else in the group to condition on) reduces
#' exactly to the marginal `N(mu_j, Psi_j)`, letting `.predict_factor_conditional()`
#' and `.route_c_draw_params()` handle this group with no special-casing.
#' @keywords internal
.fit_group_model_factor_degenerate <- function(target_analytes, safe_analytes,
                                                base, pc_wide, pc_cols,
                                                log_floors, group_name = "group",
                                                post_ndraws = 1000L) {
  a <- target_analytes[[1]]
  s <- unname(safe_analytes[a])

  rhs <- paste(paste0("s(", pc_cols, ")"), collapse = " + ")
  dat <- base[base$analyte == a, c("sample_id", "detected", "lv"), drop = FALSE]
  dat <- dplyr::left_join(dat, pc_wide, by = "sample_id")
  dat <- dat[stats::complete.cases(dat[, pc_cols, drop = FALSE]), , drop = FALSE]

  n_detected <- sum(dat$detected)
  if (n_detected < 4L) {
    cli::cli_abort(c(
      "Route C Stage-1 GAM for {.val {a}} needs at least 4 detected \\
       observations with PC scores; got {n_detected}."
    ))
  }

  cli::cli_inform(c(
    "i" = "{group_name} (factor, single-analyte {.val {a}}): fitting a \\
           Stage-1-only marginal model -- no cross-analyte coupling is \\
           possible with a single target analyte."
  ))

  gam_formula <- stats::as.formula(paste0("lv ~ ", rhs))
  gam_fit <- mgcv::gam(gam_formula, data = dat[dat$detected, , drop = FALSE],
                      family = stats::gaussian())

  detected_resid <- dat$lv[dat$detected] -
    as.numeric(stats::predict(gam_fit, newdata = dat[dat$detected, , drop = FALSE],
                              type = "response"))
  resid_var <- stats::var(detected_resid)

  list(
    gams            = stats::setNames(list(gam_fit), a),
    Lambda          = matrix(numeric(0), nrow = 1L, ncol = 0L, dimnames = list(s, NULL)),
    Psi             = stats::setNames(resid_var, s),
    k               = 0L,
    fit             = NULL,
    analytes        = target_analytes,
    safe_names      = safe_analytes,
    pc_cols         = pc_cols,
    wide_sample_ids = unique(dat$sample_id),
    impute_method   = "factor",
    log_floors      = log_floors,
    post_ndraws     = post_ndraws,
    degenerate      = TRUE
  )
}


## ----------------------------------------------------------------------------
## Marginal method -- per-analyte censored GAM, no cross-analyte factor
## ----------------------------------------------------------------------------

#' Fit the marginal (no-borrowing) group model
#'
#' One left-censored GAM per target analyte (`lv ~ s(PC1) + ...` on detected
#' observations), with NO shared factor. Uses \pkg{mgcv} only -- no brms, no
#' cmdstanr. The "borrowing" the factor model attempts is skipped entirely, so
#' this is the robust choice on panels where cross-analyte correlation is weak
#' or ragged (sparse analytes' spurious loadings can't mis-condition anything).
#' @keywords internal
.fit_group_model_marginal <- function(target_analytes, safe_analytes, base,
                                      pc_wide, pc_cols, log_floors,
                                      group_name = "group") {
  J   <- length(target_analytes)
  rhs <- paste(paste0("s(", pc_cols, ")"), collapse = " + ")
  gams <- stats::setNames(vector("list", J), target_analytes)

  for (a in target_analytes) {
    dat <- base[base$analyte == a, c("sample_id", "detected", "lv"), drop = FALSE]
    dat <- dplyr::left_join(dat, pc_wide, by = "sample_id")
    dat <- dat[stats::complete.cases(dat[, pc_cols, drop = FALSE]), , drop = FALSE]
    n_detected <- sum(dat$detected)
    if (n_detected < 4L) {
      cli::cli_abort(c(
        "Marginal Stage-1 GAM for {.val {a}} needs at least 4 detected \\
         observations with PC scores; got {n_detected}."
      ))
    }
    gams[[a]] <- mgcv::gam(stats::as.formula(paste0("lv ~ ", rhs)),
                          data = dat[dat$detected, , drop = FALSE],
                          family = stats::gaussian())
  }

  cli::cli_inform(c(
    "i" = "{group_name} (marginal): {J} per-analyte censored GAM{?s}, no \\
           shared factor (mgcv only, no Stan)."
  ))

  list(
    gams          = gams,
    analytes      = target_analytes,
    safe_names    = safe_analytes,
    pc_cols       = pc_cols,
    log_floors    = log_floors,
    impute_method = "marginal"
  )
}

#' Posterior-predictive prediction for the marginal method
#'
#' Per analyte, per cell that needs imputing (BDL or entirely absent): draw the
#' GAM mean with parameter uncertainty (`beta ~ N(coef, Vp)`, `Vp` unconditional
#' so smoothing-parameter uncertainty is included) plus residual noise
#' `N(0, sig2)`; BDL cells draw the residual left-truncated so the value never
#' exceeds `log(DL)`. This is a proper posterior predictive (mean + residual
#' uncertainty), not the plug-in residual-only draw. Returns the same `pm_long`
#' shape as the other predictors so `.predict_and_merge()` is unchanged.
#' @keywords internal
.predict_marginal <- function(group, pc_wide, df_eligible, return,
                              ndraws, batch_size) {
  target_analytes <- group$analytes
  log_floors <- group$log_floors
  N <- if (!is.null(ndraws)) as.integer(ndraws) else 1000L

  tv <- df_eligible |>
    dplyr::filter(.data$analyte %in% .env$target_analytes) |>
    dplyr::select("sample_id", "analyte", "value", "detected") |>
    dplyr::group_by(.data$sample_id, .data$analyte) |>
    dplyr::slice(1L) |>
    dplyr::ungroup()

  all_elig <- pc_wide$sample_id
  out_list <- list()

  for (a in target_analytes) {
    gam_a   <- group$gams[[a]]
    floor_a <- log_floors[[a]]
    rows_a  <- tv[tv$analyte == a, ]
    bdl     <- rows_a[!rows_a$detected, , drop = FALSE]
    missing_ids <- setdiff(all_elig, rows_a$sample_id)

    pred_ids <- c(bdl$sample_id, missing_ids)
    if (length(pred_ids) == 0L) next
    # Absolute log(DL) bound for BDL cells; NA (unbounded) for missing cells.
    upper <- c(log(pmax(bdl$value, floor_a)), rep(NA_real_, length(missing_ids)))

    pc_p <- pc_wide[match(pred_ids, pc_wide$sample_id), , drop = FALSE]
    Xp   <- stats::predict(gam_a, newdata = pc_p, type = "lpmatrix")
    sig  <- sqrt(gam_a$sig2)
    is_bdl <- !is.na(upper)

    if (return == "point") {
      mu  <- as.numeric(Xp %*% stats::coef(gam_a))
      est <- mu
      if (any(is_bdl)) {
        if (!requireNamespace("truncnorm", quietly = TRUE))
          cli::cli_abort(c("Point-mode BDL bounding needs the {.pkg truncnorm} package.",
                           "i" = "Install with {.code install.packages(\"truncnorm\")}."))
        est[is_bdl] <- truncnorm::etruncnorm(a = -Inf, b = upper[is_bdl],
                                             mean = mu[is_bdl], sd = sig)
      }
      out_list[[a]] <- tibble::tibble(sample_id = pred_ids, analyte = a,
                                      .post_mean = exp(est))
    } else {
      # Proper posterior predictive: draw sigma^2 from its scaled-inverse-chi^2
      # posterior (flat prior; df = residual df) so the predictive is Student-t,
      # not Normal -- this fattens the tails and calibrates coverage. The mean
      # draw (beta ~ N(coef, Vp)) is scaled by the same per-draw sigma so beta is
      # marginally multivariate-t and coherent with the residual. Falls back to
      # the plug-in sigma^2 when the residual df is unusable.
      cf       <- stats::coef(gam_a)
      Vp       <- stats::vcov(gam_a, unconditional = TRUE)
      sig2_hat <- gam_a$sig2
      df_res   <- gam_a$df.residual
      sig2_d <- if (is.null(df_res) || !is.finite(df_res) || df_res < 2)
        rep(sig2_hat, N) else df_res * sig2_hat / stats::rchisq(N, df_res)
      sd_d    <- sqrt(sig2_d)                                   # per-draw residual sd
      scale_d <- sd_d / sqrt(sig2_hat)                          # per-draw coef scale

      eta_mean <- as.numeric(Xp %*% cf)                         # n_cells
      B0       <- mgcv::rmvn(N, rep(0, length(cf)), Vp)         # centred coef draws
      eta_dev  <- sweep(Xp %*% t(B0), 2, scale_d, "*")          # n_cells x N
      eta      <- eta_dev + eta_mean                            # eta_mean per column
      eps      <- sweep(matrix(stats::rnorm(length(eta)), nrow = nrow(eta)),
                        2, sd_d, "*")                           # eps[,d] ~ N(0, sig2_d)
      if (any(is_bdl)) {
        if (!requireNamespace("truncnorm", quietly = TRUE))
          cli::cli_abort(c("BDL draw truncation needs the {.pkg truncnorm} package.",
                           "i" = "Install with {.code install.packages(\"truncnorm\")}."))
        ub  <- upper[is_bdl] - eta[is_bdl, , drop = FALSE]      # per-draw residual bound
        sdb <- rep(sd_d, each = sum(is_bdl))                    # per-draw sd, col-major
        eps[is_bdl, ] <- matrix(
          truncnorm::rtruncnorm(length(ub), a = -Inf, b = as.vector(ub),
                                mean = 0, sd = sdb),
          nrow = sum(is_bdl))
      }
      val <- exp(eta + eps)                                     # n_cells x N
      out_list[[a]] <- tibble::tibble(
        sample_id   = rep(pred_ids, times = N),
        analyte     = a,
        draw_id     = rep(seq_len(N), each = length(pred_ids)),
        .post_value = as.numeric(val))
    }
  }

  out <- dplyr::bind_rows(out_list)
  if (nrow(out) == 0L) {
    out <- if (return == "point")
      tibble::tibble(sample_id = character(), analyte = character(),
                     .post_mean = double())
    else
      tibble::tibble(sample_id = character(), analyte = character(),
                     draw_id = integer(), .post_value = double())
  }
  out
}


## ----------------------------------------------------------------------------
## Prediction -- closed-form conditional draws from the fitted factor model
## ----------------------------------------------------------------------------

#' Extract per-draw (Lambda, Psi) pairs from a fitted Route C group
#'
#' @param group A fitted `factor`-method group (from `.fit_group_model_factor()`).
#' @param ndraws Optional cap on the number of posterior draws to use (the
#'   first `ndraws`); `NULL` uses every draw.
#' @return `list(Lambda = list of J x k matrices, Psi = list of length-J
#'   vectors)`, one entry per used draw, plus `n` (the number used).
#' @keywords internal
.route_c_draw_params <- function(group, ndraws = NULL, return = "draws") {
  J <- length(group$analytes)
  k <- group$k
  safe_vec <- unname(group$safe_names[group$analytes])

  if (isTRUE(group$degenerate)) {
    # No Stan fit exists (single-analyte group): Lambda/Psi are a fixed plug-in
    # point (J x 0 loadings, so .factor_condition()'s marginal branch applies).
    # In draws mode, replicate across `post_ndraws` -- matched to the Stan
    # groups' posterior draw count (chains * (iter - warmup)) so a mixed
    # Stan+degenerate model emits ONE consistent draw domain -- letting
    # .factor_condition_draw()'s own per-draw Gaussian noise give a proper
    # predictive spread. In point mode the params are constant across draws, so
    # a single evaluation suffices (no redundant averaging over identical draws).
    post <- group$post_ndraws %||% 1000L
    n <- if (return == "point") 1L
         else if (is.null(ndraws)) post
         else min(as.integer(ndraws), post)
    return(list(
      Lambda = rep(list(group$Lambda), n),
      Psi    = rep(list(group$Psi), n),
      n      = n
    ))
  }

  .require_cmdstanr()
  draws_mat <- posterior::as_draws_matrix(group$fit$draws(c("Lambda", "psi")))
  vars <- colnames(draws_mat)
  lam_cols <- grep("^Lambda\\[", vars)
  psi_cols <- grep("^psi\\[", vars)

  lam_ij <- regmatches(vars[lam_cols],
                       regexpr("(?<=\\[)[0-9]+,[0-9]+(?=\\])", vars[lam_cols], perl = TRUE))
  lam_ij <- do.call(rbind, strsplit(lam_ij, ","))
  lam_i  <- as.integer(lam_ij[, 1]); lam_j <- as.integer(lam_ij[, 2])

  psi_idx <- as.integer(regmatches(vars[psi_cols],
                       regexpr("(?<=\\[)[0-9]+(?=\\])", vars[psi_cols], perl = TRUE)))

  D_total  <- nrow(draws_mat)
  draw_idx <- if (is.null(ndraws)) seq_len(D_total) else seq_len(min(as.integer(ndraws), D_total))

  Lambda_list <- vector("list", length(draw_idx))
  Psi_list    <- vector("list", length(draw_idx))
  for (dd in seq_along(draw_idx)) {
    d <- draw_idx[dd]
    Lam_d <- matrix(0, nrow = J, ncol = k, dimnames = list(safe_vec, NULL))
    Lam_d[cbind(lam_i, lam_j)] <- draws_mat[d, lam_cols]
    Lambda_list[[dd]] <- Lam_d
    psi_d <- stats::setNames(numeric(J), safe_vec)
    psi_d[psi_idx] <- draws_mat[d, psi_cols]
    Psi_list[[dd]] <- psi_d
  }

  list(Lambda = Lambda_list, Psi = Psi_list, n = length(draw_idx))
}

#' Closed-form conditional prediction for the Route C factor model
#'
#' For each eligible sample, builds the residual vector `y` (observed minus
#' the Stage-1 GAM mean for detected target cells; `NA` for BDL/missing
#' cells), then applies `.factor_condition()` / `.factor_condition_draw()`
#' per posterior draw of `(Lambda, Psi)` to predict the missing cells --
#' conditioned on this sample's *own* observed analytes (finding 3), with any
#' BDL target truncated at `log(DL)` (findings 1-2). Mirrors
#' `.predict_factor_long()`'s `pm_long` return shape so
#' `.predict_and_merge()`'s merge/fabricate logic is unchanged; only rows for
#' cells that actually need prediction (BDL or entirely absent) are emitted.
#' @param group A fitted `factor`-method group.
#' @param pc_wide Per-sample PC scores for the eligible samples (`sample_id`
#'   + `PC*` columns).
#' @param df_eligible The long-format input rows for the eligible samples
#'   (supplies each sample's own observed target values).
#' @inheritParams .predict_factor_long
#' @keywords internal
.predict_factor_conditional <- function(group, pc_wide, df_eligible, return,
                                        ndraws, batch_size) {
  target_analytes <- group$analytes
  safe_vec  <- unname(group$safe_names[target_analytes])
  J         <- length(target_analytes)
  log_floors <- group$log_floors

  # mu_j(X_i) for every eligible sample x target analyte, from the Stage-1 GAMs.
  mu_wide <- pc_wide["sample_id"]
  for (i in seq_len(J)) {
    mu_wide[[safe_vec[i]]] <- as.numeric(
      stats::predict(group$gams[[target_analytes[i]]], newdata = pc_wide, type = "response")
    )
  }

  # Observed residual (detected) / censoring bound (BDL) per (sample, analyte);
  # NA (unobserved -> to predict) for cells with no detected value.
  obs_wide   <- mu_wide["sample_id"]
  upper_wide <- mu_wide["sample_id"]
  for (s in safe_vec) { obs_wide[[s]] <- NA_real_; upper_wide[[s]] <- NA_real_ }

  target_vals <- df_eligible |>
    dplyr::filter(.data$analyte %in% .env$target_analytes) |>
    dplyr::select("sample_id", "analyte", "value", "detected") |>
    dplyr::group_by(.data$sample_id, .data$analyte) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      safe = unname(group$safe_names[.data$analyte]),
      lv   = log(pmax(.data$value, log_floors[.data$analyte]))
    )

  for (i in seq_len(J)) {
    a <- target_analytes[i]; s <- safe_vec[i]
    rows_a <- target_vals[target_vals$analyte == a, ]
    idx <- match(rows_a$sample_id, mu_wide$sample_id)
    keep <- !is.na(idx)
    idx <- idx[keep]; rows_a <- rows_a[keep, ]
    mu_a <- mu_wide[[s]][idx]

    det <- rows_a$detected
    obs_wide[[s]][idx[det]]   <- rows_a$lv[det]  - mu_a[det]
    upper_wide[[s]][idx[!det]] <- rows_a$lv[!det] - mu_a[!det]
  }

  # -- Posterior draws of (Lambda, Psi) ---------------------------------------
  params <- .route_c_draw_params(group, ndraws, return)
  n_draws <- params$n

  n_new <- nrow(mu_wide)
  bs <- if (is.null(batch_size)) n_new else max(1L, as.integer(batch_size))
  batches <- split(seq_len(n_new), ceiling(seq_len(n_new) / bs))

  predict_rows <- function(rows) {
    purrr::map_dfr(rows, function(i) {
      y <- as.numeric(obs_wide[i, safe_vec])
      miss <- which(is.na(y))
      if (length(miss) == 0L) return(NULL)

      sid    <- mu_wide$sample_id[i]
      analyt <- target_analytes[miss]
      mu0    <- rep(0, J)
      mu_i   <- as.numeric(mu_wide[i, safe_vec])[miss]
      upper_i <- as.numeric(upper_wide[i, safe_vec])[miss]

      if (return == "point") {
        has_bound <- !is.na(upper_i)
        if (any(has_bound) && !requireNamespace("truncnorm", quietly = TRUE)) {
          cli::cli_abort(c(
            "Point-mode BDL bounding needs the {.pkg truncnorm} package.",
            "i" = "Install with {.code install.packages(\"truncnorm\")}."
          ))
        }
        acc <- numeric(length(miss))
        for (d in seq_len(n_draws)) {
          cond <- .factor_condition(y, mu0, params$Lambda[[d]], params$Psi[[d]])
          contrib <- as.numeric(cond$mean)
          if (any(has_bound)) {
            sd_d <- sqrt(diag(as.matrix(cond$cov)))
            contrib[has_bound] <- truncnorm::etruncnorm(
              a = -Inf, b = upper_i[has_bound],
              mean = cond$mean[has_bound], sd = sd_d[has_bound]
            )
          }
          acc <- acc + contrib
        }
        tibble::tibble(
          sample_id = sid, analyte = analyt,
          .post_mean = exp(mu_i + acc / n_draws)
        )
      } else {
        vals <- vector("list", n_draws)
        for (d in seq_len(n_draws)) {
          draw_d <- .factor_condition_draw(
            y, mu0, params$Lambda[[d]], params$Psi[[d]], ndraws = 1L, upper = upper_i
          )
          vals[[d]] <- exp(mu_i + as.numeric(draw_d))
        }
        tibble::tibble(
          sample_id   = sid,
          analyte     = rep(analyt, times = n_draws),
          draw_id     = rep(seq_len(n_draws), each = length(miss)),
          .post_value = unlist(vals, use.names = FALSE)
        )
      }
    })
  }

  out <- purrr::map_dfr(batches, predict_rows)

  # Every eligible cell was already detected (nothing to predict): return a
  # correctly-typed empty frame so the join in .predict_and_merge() still
  # works rather than erroring on a columnless tibble.
  if (nrow(out) == 0L) {
    out <- if (return == "point") {
      tibble::tibble(sample_id = character(), analyte = character(),
                     .post_mean = double())
    } else {
      tibble::tibble(sample_id = character(), analyte = character(),
                     draw_id = integer(), .post_value = double())
    }
  }
  out
}
