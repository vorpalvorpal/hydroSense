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

#' Compiled Stage-2 censored-factor Stan model (cached for the session)
#' @keywords internal
.route_c_stan_model <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    cli::cli_abort(c(
      "The Route C {.val factor} imputation method needs the {.pkg cmdstanr} \\
       package (plus a working CmdStan install).",
      "i" = "Install with {.code install.packages(\"cmdstanr\", repos = \\
             \"https://stan-dev.r-universe.dev\")} then {.code \\
             cmdstanr::install_cmdstan()}."
    ))
  }
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
#' @exportS3Method posterior::rhat
rhat.route_c_stanfit <- function(x, ...) {
  s <- x$summary()
  stats::setNames(s$rhat, s$variable)
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
                                     group_name = "group", k = NULL, ...) {
  J <- length(target_analytes)
  safe_vec <- unname(safe_analytes[target_analytes])

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
    log_floors      = log_floors
  )
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
.route_c_draw_params <- function(group, ndraws = NULL) {
  J <- length(group$analytes)
  k <- group$k
  safe_vec <- unname(group$safe_names[group$analytes])

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
  params <- .route_c_draw_params(group, ndraws)
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
        acc <- numeric(length(miss))
        for (d in seq_len(n_draws)) {
          cond <- .factor_condition(y, mu0, params$Lambda[[d]], params$Psi[[d]])
          acc  <- acc + cond$mean
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
