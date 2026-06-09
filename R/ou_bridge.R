## ============================================================================
## ou_bridge.R  --  OU bridge for daily residual-state uncertainty (issue #16)
## ============================================================================
##
## OU process: dS = -θ·S dt + σ dW, stationary,
## covariance kernel k(s,t) = γ·exp(-θ|t-s|), marginal variance γ = σ²/(2θ).
##
## These functions supply mean-zero OU-bridge fluctuations to amspaf_daily().
## The existing deterministic .interp_residual() blend stays the centre line;
## draws are pinned exactly to 0 at anchor (grab-sample) dates.
##
## As θ → 0, the OU bridge converges exactly to the Brownian bridge — the code
## detects this automatically; no separate fallback model is needed.
##
## Public surface (package-internal):
##   .estimate_ou_params    MoM (θ, σ², γ) from residual anchor observations
##   .ou_bridge_factors     precompute per-gap Cholesky cache over a daily grid
##   .ou_bridge_draw        draw one coherent mean-zero path from the cache


## ── MoM parameter estimation ─────────────────────────────────────────────────

#' Method-of-moments OU parameter estimation from an anchor residual series
#'
#' Estimates the Ornstein-Uhlenbeck parameters θ (mean-reversion rate in
#' units of 1/day), σ² (diffusion coefficient in [S²/day]), and γ (marginal
#' variance in [S²]) from an irregularly-spaced series of residual-state S
#' values at grab-sample anchor dates.
#'
#' **Estimators:**
#' * γ̂ = sample variance of anchor S values.
#' * σ̂² = median[(ΔS)²/Δt] over consecutive anchor pairs (robust quadratic
#'   variation; median is resistant to outlier jumps at distant anchor pairs).
#' * θ̂ = σ̂²/(2γ̂)  (from the OU marginal-variance identity γ = σ²/(2θ)).
#'
#' Returns `degenerate = TRUE` when fewer than two finite anchors are
#' available, when S is essentially constant (γ̂ ≈ 0), or when consecutive
#' anchors show no detectable diffusion (σ̂² ≈ 0).  A degenerate result means
#' [.ou_bridge_draw()] returns a zero path (no added uncertainty).
#'
#' @param anchor_dates Date (or Date-coercible) vector of anchor dates.
#' @param anchor_S Numeric; residual state S at each anchor, same length.
#' @param scale Positive multiplier applied to γ and σ² (envelope-width
#'   tuning knob, default 1).  θ is unchanged (correlation structure
#'   preserved).
#' @return Named list: `theta` (1/day), `sigma2` ([S²/day]),
#'   `gamma` ([S²]), `degenerate` (logical).
#' @keywords internal
.estimate_ou_params <- function(anchor_dates, anchor_S, scale = 1) {
  checkmate::assert_number(scale, lower = 0, finite = TRUE)
  eps <- sqrt(.Machine$double.eps)
  n   <- length(anchor_S)

  if (n < 2L || !any(is.finite(anchor_S))) {
    return(list(theta = 0, sigma2 = 0, gamma = 0, degenerate = TRUE))
  }

  ord     <- order(as.Date(anchor_dates))
  a_S     <- anchor_S[ord]
  a_dates <- as.Date(anchor_dates)[ord]

  gamma_hat <- stats::var(a_S, na.rm = TRUE)
  if (!is.finite(gamma_hat) || gamma_hat < eps) {
    return(list(theta = 0, sigma2 = 0, gamma = 0, degenerate = TRUE))
  }

  dt   <- as.numeric(diff(a_dates))
  dS   <- diff(a_S)
  qv   <- dS^2 / pmax(dt, 1)             # (ΔS)²/Δt; pmax guards same-day pairs
  sigma2_hat <- stats::median(qv, na.rm = TRUE)
  if (!is.finite(sigma2_hat) || sigma2_hat < eps) {
    return(list(theta = 0, sigma2 = 0, gamma = 0, degenerate = TRUE))
  }

  theta_hat <- sigma2_hat / (2 * gamma_hat)   # always ≥ 0

  list(
    theta      = theta_hat,
    sigma2     = sigma2_hat * scale,
    gamma      = gamma_hat  * scale,
    degenerate = FALSE
  )
}


## ── Per-gap Cholesky cache ────────────────────────────────────────────────────

#' Precompute Cholesky factors for OU bridge draws over a daily grid
#'
#' For each inter-anchor gap and for the leading/trailing extrapolation zones,
#' builds the conditional covariance of the latent residual state S given
#' zero values at the bracketing anchor dates, then Cholesky-factorises it.
#'
#' Call this **once** per (analyte) after fitting the target model.  Each of
#' the N per-draw calls to [.ou_bridge_draw()] then requires only one
#' matrix-vector product per gap.
#'
#' **OU conditional covariance** for interior days `{t₁,…,tₘ}` in gap `[a,b]`
#' given `S(a) = S(b) = 0`:
#' ```
#' Σ_cond[i,j] = γ·exp(-θ|tᵢ−tⱼ|) −
#'   γ·(uᵢuⱼ + vᵢvⱼ − s(uᵢvⱼ + vᵢuⱼ)) / (1−s²)
#' ```
#' where `uᵢ = exp(-θ(tᵢ−a))`, `vᵢ = exp(-θ(b−tᵢ))`, `s = exp(-θ(b−a))`.
#' As `θ·(b−a) → 0` this equals the Brownian bridge covariance
#' `σ²·min(tᵢ−a, tⱼ−a)·(b−max(tᵢ,tⱼ)) / (b−a)` — the code switches
#' automatically when `θ·L < 1e-8`.
#'
#' @param anchor_dates Date vector of anchor dates (need not be sorted; need
#'   not all appear in `target_dates`).
#' @param target_dates Date vector of the full daily grid.
#' @param theta,sigma2,gamma OU parameters from [.estimate_ou_params()].
#' @return A list (the "factors" object) consumed by [.ou_bridge_draw()].
#' @keywords internal
.ou_bridge_factors <- function(anchor_dates, target_dates, theta, sigma2, gamma) {
  a_dates <- sort(as.Date(anchor_dates))
  t_dates <- sort(as.Date(target_dates))
  na      <- length(a_dates)
  nt      <- length(t_dates)

  if (gamma <= 0 || nt == 0L) {
    return(list(
      n_target         = nt,
      anchor_in_target = rep(FALSE, nt),
      theta = theta, sigma2 = sigma2, gamma = gamma,
      gaps = list()
    ))
  }

  anchor_in_target <- t_dates %in% a_dates

  # ── Cholesky helper for one inter-anchor bridge ──────────────────────────
  # t_int: days-from-left-anchor for each interior target date (numeric, >0)
  # L_gap: gap length in days (numeric, >0)
  # Returns lower-triangular Cholesky L s.t. Sigma_cond = L L^T, or NULL.
  .chol_bridge <- function(t_int, L_gap) {
    m <- length(t_int)
    ta <- t_int
    tb <- L_gap - t_int

    if (theta * L_gap < 1e-8) {
      # Brownian-bridge limit: Cov(ti,tj) = σ²·min(ta_i,ta_j)·(L−max(ta_i,ta_j))/L
      min_ta <- outer(ta, ta, pmin)
      max_ta <- outer(ta, ta, pmax)
      Sigma  <- (sigma2 / L_gap) * min_ta * (L_gap - max_ta)
    } else {
      u <- exp(-theta * ta)
      v <- exp(-theta * tb)
      s <- exp(-theta * L_gap)
      d <- 1 - s^2
      Sigma_II <- outer(ta, ta, function(x, y) gamma * exp(-theta * abs(x - y)))
      Sigma <- Sigma_II -
        gamma * (outer(u, u) + outer(v, v) -
                   s * (outer(u, v) + outer(v, u))) / d
    }

    Sigma <- (Sigma + t(Sigma)) / 2          # symmetrise numerically

    tryCatch(
      t(chol(Sigma)),
      error = function(e) {
        jitter <- pmax(max(abs(diag(Sigma))) * 1e-8, .Machine$double.eps)
        diag(Sigma) <<- diag(Sigma) + jitter
        tryCatch(t(chol(Sigma)), error = function(e2) NULL)
      }
    )
  }

  # ── One-sided conditional SD (leading / trailing) ────────────────────────
  # Var[S(t) | S(anchor) = 0] = γ(1 − exp(−2θΔt))   [OU]
  #                           ≈ σ²Δt                  [Brownian limit]
  .one_sided_sd <- function(dt_from_anchor) {
    if (theta < 1e-8) {
      sqrt(sigma2 * dt_from_anchor)
    } else {
      sqrt(gamma * (1 - exp(-2 * theta * dt_from_anchor)))
    }
  }

  gaps <- list()

  if (na > 0L) {
    # Leading zone: target dates before the first anchor
    lead_idx <- which(t_dates < a_dates[1L])
    if (length(lead_idx) > 0L) {
      dt <- as.numeric(a_dates[1L] - t_dates[lead_idx])
      gaps[["leading"]] <- list(
        type   = "one_sided",
        idx    = lead_idx,
        sd_vec = .one_sided_sd(dt)
      )
    }

    # Inter-anchor bridge gaps
    if (na >= 2L) {
      for (k in seq_len(na - 1L)) {
        a_date  <- a_dates[k]
        b_date  <- a_dates[k + 1L]
        int_idx <- which(t_dates > a_date & t_dates < b_date)
        gname   <- paste0("gap_", k)

        if (length(int_idx) == 0L) {
          gaps[[gname]] <- list(type = "bridge", idx = integer(0), chol_L = NULL)
          next
        }

        L_gap  <- as.numeric(b_date - a_date)
        t_int  <- as.numeric(t_dates[int_idx] - a_date)
        L_ch   <- .chol_bridge(t_int, L_gap)

        if (is.null(L_ch)) {
          # Cholesky failed: fall back to marginal N(0, γ) draws per interior day
          gaps[[gname]] <- list(
            type   = "one_sided",
            idx    = int_idx,
            sd_vec = rep(sqrt(gamma), length(int_idx))
          )
        } else {
          gaps[[gname]] <- list(type = "bridge", idx = int_idx, chol_L = L_ch)
        }
      }
    }

    # Trailing zone: target dates after the last anchor
    trail_idx <- which(t_dates > a_dates[na])
    if (length(trail_idx) > 0L) {
      dt <- as.numeric(t_dates[trail_idx] - a_dates[na])
      gaps[["trailing"]] <- list(
        type   = "one_sided",
        idx    = trail_idx,
        sd_vec = .one_sided_sd(dt)
      )
    }
  }

  list(
    n_target         = nt,
    anchor_in_target = anchor_in_target,
    theta            = theta,
    sigma2           = sigma2,
    gamma            = gamma,
    gaps             = gaps
  )
}


## ── Per-draw path sampler ─────────────────────────────────────────────────────

#' Draw one coherent mean-zero OU-bridge fluctuation path
#'
#' Uses precomputed Cholesky factors from [.ou_bridge_factors()] to generate
#' one realisation of the latent residual-state fluctuation ε(t) over the
#' full daily grid.
#'
#' All anchor dates carry ε = 0 (pinch condition).  Within each inter-anchor
#' gap, interior days are drawn jointly from the conditional Gaussian
#' distribution via the precomputed lower-triangular Cholesky `L`:
#' `ε_gap = L %*% rnorm(m)`.  This is O(m) per draw (the O(m³) Cholesky
#' was paid once in [.ou_bridge_factors()]).
#'
#' For leading/trailing extrapolation zones, each day is drawn independently
#' from `N(0, SD²)` where SD grows from 0 at the anchor toward √γ.
#'
#' @param factors Output of [.ou_bridge_factors()].
#' @return Numeric vector of length `factors$n_target`: ε for each target
#'   date in the sorted order used when calling [.ou_bridge_factors()].
#'   Call [set.seed()] before this function for reproducibility.
#' @keywords internal
.ou_bridge_draw <- function(factors) {
  epsilon <- numeric(factors$n_target)    # 0 everywhere; anchors stay 0

  for (g in factors$gaps) {
    idx <- g$idx
    if (length(idx) == 0L) next

    if (g$type == "bridge") {
      epsilon[idx] <- as.numeric(g$chol_L %*% stats::rnorm(length(idx)))
    } else {                               # one_sided (leading / trailing / fallback)
      epsilon[idx] <- stats::rnorm(length(idx), 0, g$sd_vec)
    }
  }

  epsilon
}
