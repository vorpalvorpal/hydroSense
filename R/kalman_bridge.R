## ============================================================================
## kalman_bridge.R -- state-space (Kalman/OU) residual smoother (issue #16)
## ============================================================================
##
## Replaces the deterministic `.interp_residual()` centre line AND the OU bridge
## with ONE coherent model: a scalar continuous-time AR(1)/Ornstein-Uhlenbeck
## latent state on the per-analyte residual S, fitted to the grab anchors and
## evaluated on a daily grid by a Kalman filter + smoother (via KFAS). The
## posterior mean is the centre line; `simulateSSM()` gives coherent draws.
##
## Model (Delta t = 1 day on the grid):
##   state:        x_{t+1} = phi x_t + w_t,  w_t ~ N(0, q_t),  phi = exp(-theta)
##   observation:  y_i     = x_{t_i} + v_i,  v_i ~ N(0, r_i)   at grab days t_i
##   q_t = gamma (1 - phi^2) * q_mult_t      (q_mult carries hydro modulation)
##
## Locality enters through the time-varying q_mult (hydrology-modulated process
## variance); the residual mean is a pure temporal smooth (S is already
## hydrology-de-trended by beta.f(hydro), so re-leaning it on hydro would
## double-count). The grid is clipped per analyte to [first_anchor, last_anchor].
##
## Public surface (package-internal):
##   .estimate_ou_kalman_params  MoM gamma + 1-D MLE theta, with fallback ladder
##   .build_kalman_model         KFAS SSModel on a grid (time-varying Q and H)
##   .kalman_smooth              posterior mean + variance per grid day (KFS)
##   .kalman_draw                simulateSSM state draws (matrix [n_grid x nsim])
##   .residual_smoother          per-analyte convenience: clip, fit, smooth
##   .loo_coverage_series        leave-one-anchor-out coverage of one S series
##   .loo_anchor_coverage        per-analyte + pooled coverage over a target_model

#' @importFrom KFAS SSModel SSMcustom KFS simulateSSM
NULL


## ── Parameter estimation ─────────────────────────────────────────────────────

#' Estimate OU parameters (gamma by MoM, theta by MLE) with a fallback ladder
#'
#' `gamma` (the marginal/stationary variance, the well-identified moment) is the
#' sample variance of the anchor residuals, times `scale`. `theta` (the
#' mean-reversion rate / inverse correlation length, weakly identified on sparse
#' data) is fitted by 1-D MLE of the exact irregular-spacing OU likelihood,
#' bounded so the correlation length stays in [0.5x, 10x] the median anchor
#' spacing.
#'
#' Tiers: `"degenerate"` (< 2 finite anchors or ~constant S -> no uncertainty),
#' `"fallback"` (>= 2 but < `n_fit_min` anchors -> theta pinned at the long-
#' correlation bound, approximating a Brownian bridge without diffuse init),
#' `"mle"` (>= `n_fit_min` -> theta by MLE).
#'
#' @param anchor_dates Date (or coercible) vector of anchor dates.
#' @param anchor_S Numeric residual S at each anchor (same length).
#' @param n_fit_min Minimum anchors to attempt the theta MLE (default 8).
#' @param scale Positive multiplier on `gamma` (the envelope-width knob;
#'   correlation length / theta is unchanged). Default 1.
#' @return list(theta, gamma, tier).
#' @keywords internal
.estimate_ou_kalman_params <- function(anchor_dates, anchor_S, n_fit_min = 8L,
                                       scale = 1) {
  eps <- sqrt(.Machine$double.eps)
  d   <- as.Date(anchor_dates)
  ok  <- is.finite(anchor_S) & !is.na(d)
  d   <- d[ok]; s <- anchor_S[ok]
  ord <- order(d); d <- d[ord]; s <- s[ord]
  n   <- length(s)

  if (n < 2L) return(list(theta = 0, gamma = 0, tier = "degenerate"))
  gamma_raw <- stats::var(s)          # fit theta against the data's own variance
  if (!is.finite(gamma_raw) || gamma_raw < eps)
    return(list(theta = 0, gamma = 0, tier = "degenerate"))
  gamma <- gamma_raw * scale          # `scale` inflates the envelope only

  med <- stats::median(diff(as.numeric(d)))
  if (!is.finite(med) || med <= 0) med <- 1
  theta_min <- 1 / (10 * med)        # correlation length 10x spacing
  theta_max <- 2 / med               # correlation length 0.5x spacing

  if (n < n_fit_min)
    return(list(theta = theta_min, gamma = gamma, tier = "fallback"))

  ## Exact irregular-spacing mean-zero OU log-likelihood (uses the RAW variance,
  ## so `scale` does not move the fitted correlation length).
  dt <- as.numeric(diff(d))
  ll <- function(theta) {
    phi <- exp(-theta * dt)
    qv  <- gamma_raw * (1 - phi^2)
    qv[qv < eps] <- eps
    sum(stats::dnorm(s[-1L], phi * s[-n], sqrt(qv), log = TRUE)) +
      stats::dnorm(s[1L], 0, sqrt(gamma_raw), log = TRUE)
  }
  opt <- stats::optimize(ll, c(theta_min, theta_max), maximum = TRUE)
  theta_hat <- opt$maximum
  if (!is.finite(theta_hat)) theta_hat <- theta_min
  list(theta = theta_hat, gamma = gamma, tier = "mle")
}


## ── KFAS model construction ──────────────────────────────────────────────────

#' Build a KFAS state-space model for one residual series on a grid
#'
#' **Invariant:** `grid_dates` must be *consecutive days* (step = 1 day). The
#' transition `phi = exp(-theta)` assumes a unit time step between adjacent grid
#' points, so a sparse/irregular grid would mis-state the correlation length.
#'
#' @param grid_dates Sorted, consecutive-day Date vector: the grid to predict on.
#' @param anchor_dates,anchor_S Observed anchors (must fall within `grid_dates`).
#' @param theta,gamma OU parameters from [.estimate_ou_kalman_params()].
#' @param r_vec Per-anchor observation variance (length = #anchors), or `NULL`
#'   (-> a tiny value relative to gamma, so the smoother pins to the anchors).
#' @param q_mult Per-grid-day positive multiplier on the process variance
#'   (hydrology modulation), length `length(grid_dates)`, or `NULL` (-> 1).
#' @return A KFAS `SSModel`, or `NULL` if `gamma <= 0` or no anchors on the grid.
#'   The grid is stored in `attr(., "grid_dates")`.
#' @keywords internal
.build_kalman_model <- function(grid_dates, anchor_dates, anchor_S, theta, gamma,
                                r_vec = NULL, q_mult = NULL) {
  grid_dates <- as.Date(grid_dates)
  n <- length(grid_dates)
  if (gamma <= 0 || n == 0L) return(NULL)

  pos <- match(as.Date(anchor_dates), grid_dates)
  keep <- !is.na(pos)
  pos <- pos[keep]; aS <- anchor_S[keep]
  if (!is.null(r_vec)) r_vec <- rep_len(r_vec, length(keep))[keep]
  if (length(pos) == 0L) return(NULL)

  phi <- exp(-theta)
  q_base <- gamma * (1 - phi^2)
  if (q_base <= 0) q_base <- gamma * sqrt(.Machine$double.eps)
  if (is.null(q_mult)) q_mult <- rep(1, n)
  Qvec <- q_base * q_mult

  if (is.null(r_vec)) r_vec <- rep(max(gamma, 1) * 1e-6, length(pos))
  r_vec <- rep_len(r_vec, length(pos))
  Hvec <- rep(max(gamma, 1), n)        # arbitrary at NA days (ignored by KFAS)
  Hvec[pos] <- r_vec

  y <- rep(NA_real_, n)
  y[pos] <- aS

  ## NOTE: SSMcustom must be called *unqualified* inside the formula — KFAS's
  ## formula parser does not resolve the `KFAS::` form. It is imported via
  ## @importFrom in the NAMESPACE.
  mod <- KFAS::SSModel(
    y ~ -1 + SSMcustom(
      Z = matrix(1), T = matrix(phi), R = matrix(1),
      Q = array(Qvec, c(1, 1, n)), a1 = 0, P1 = matrix(gamma)),
    H = array(Hvec, c(1, 1, n)))
  attr(mod, "grid_dates") <- grid_dates
  mod
}


## ── Smoothing and draws ──────────────────────────────────────────────────────

#' Posterior (smoothed) mean and variance of the residual on the grid
#' @param model A KFAS model from [.build_kalman_model()].
#' @return list(mean, var), each length `length(grid_dates)`.
#' @keywords internal
.kalman_smooth <- function(model) {
  out  <- KFAS::KFS(model, filtering = "state", smoothing = "state")
  mean <- as.numeric(out$alphahat)
  v    <- as.numeric(out$V)
  v[!is.finite(v) | v < 0] <- 0
  list(mean = mean, var = v)
}

#' Draw coherent residual trajectories via the simulation smoother
#' @param model A KFAS model from [.build_kalman_model()].
#' @param nsim Number of draws.
#' @return Numeric matrix `[n_grid x nsim]`. Call [set.seed()] beforehand for
#'   reproducibility.
#' @keywords internal
.kalman_draw <- function(model, nsim = 1L) {
  dr <- KFAS::simulateSSM(model, type = "states", nsim = nsim)
  matrix(dr[, 1L, ], nrow = dim(dr)[1L], ncol = nsim)
}


## ── Per-analyte convenience ──────────────────────────────────────────────────

#' Fit + smooth the residual for one analyte over a (clipped) daily grid
#'
#' Clips `target_dates` to the analyte's grab span `[first_anchor, last_anchor]`
#' (no extrapolation beyond the grabs), fits the OU parameters, builds the KFAS
#' model with optional hydrology-modulated process variance, and returns the
#' smoothed posterior. For degenerate series returns a flat mean with zero
#' variance.
#'
#' @param anchor_dates,anchor_S Observed anchors.
#' @param target_dates Candidate daily grid (will be clipped to the grab span).
#' @param z_hydro Optional standardised hydro feature aligned to `target_dates`
#'   (same length); drives `q_mult = exp(kappa * z)`. `NULL` -> no modulation.
#' @param kappa Hydrology process-variance sensitivity (default 0.5).
#' @param n_fit_min,scale Passed to [.estimate_ou_kalman_params()].
#' @param r_vec Per-anchor observation variance (default tiny).
#' @return list(grid_dates, params, model, mean, var). `model` is `NULL` for the
#'   degenerate path.
#' @keywords internal
.residual_smoother <- function(anchor_dates, anchor_S, target_dates,
                               z_hydro = NULL, kappa = 0.5, n_fit_min = 8L,
                               scale = 1, r_vec = NULL) {
  d  <- as.Date(anchor_dates)
  ok <- is.finite(anchor_S) & !is.na(d)
  d  <- d[ok]; aS <- anchor_S[ok]
  if (!is.null(r_vec)) r_vec <- rep_len(r_vec, length(ok))[ok]
  tgt <- sort(unique(as.Date(target_dates)))

  if (length(d) == 0L)
    return(list(grid_dates = tgt[0], params = list(tier = "degenerate"),
                model = NULL, mean = numeric(0), var = numeric(0)))

  span <- range(d)
  grid <- tgt[tgt >= span[1L] & tgt <= span[2L]]
  grid <- sort(unique(c(grid, d[d >= span[1L] & d <= span[2L]])))
  if (length(grid) == 0L) grid <- sort(unique(d))

  params <- .estimate_ou_kalman_params(d, aS, n_fit_min = n_fit_min,
                                       scale = scale)

  if (params$tier == "degenerate" || params$gamma <= 0) {
    return(list(grid_dates = grid, params = params, model = NULL,
                mean = rep(mean(aS), length(grid)),
                var  = rep(0, length(grid))))
  }

  q_mult <- NULL
  if (!is.null(z_hydro)) {
    z_grid <- z_hydro[match(grid, tgt)]
    z_grid[!is.finite(z_grid)] <- 0
    q_mult <- exp(kappa * z_grid)
  }

  model <- .build_kalman_model(grid, d, aS, params$theta, params$gamma,
                               r_vec = r_vec, q_mult = q_mult)
  if (is.null(model))
    return(list(grid_dates = grid, params = params, model = NULL,
                mean = rep(mean(aS), length(grid)),
                var  = rep(0, length(grid))))

  sm <- .kalman_smooth(model)
  list(grid_dates = grid, params = params, model = model,
       mean = sm$mean, var = sm$var)
}


## ── Coverage validation harness ──────────────────────────────────────────────

#' Leave-one(-block)-out coverage of one residual series
#'
#' Holds out each anchor (or each run of `block` consecutive anchors), refits the
#' smoother on the rest, and checks whether the held-out value falls inside the
#' predictive interval (state posterior + observation noise). Plain LOO mainly
#' tests the pinch zone; `block >= 2` probes mid-gap width.
#'
#' @param anchor_dates,anchor_S The series.
#' @param interval Nominal coverage (default 0.9).
#' @param block Number of consecutive anchors held out per fold (default 1).
#' @param n_fit_min,scale Passed through to the smoother.
#' @return list(coverage, mean_width, n, n_held). `coverage`/`mean_width` are
#'   `NA` when too few anchors remain to fit.
#' @keywords internal
.loo_coverage_series <- function(anchor_dates, anchor_S, interval = 0.9,
                                 block = 1L, n_fit_min = 8L, scale = 1) {
  d  <- as.Date(anchor_dates)
  ok <- is.finite(anchor_S) & !is.na(d)
  d  <- d[ok]; s <- anchor_S[ok]
  ord <- order(d); d <- d[ord]; s <- s[ord]
  n <- length(s)
  z <- stats::qnorm(1 - (1 - interval) / 2)

  if (n < 4L) return(list(coverage = NA_real_, mean_width = NA_real_,
                          n = n, n_held = 0L))

  folds <- split(seq_len(n), ceiling(seq_len(n) / block))
  hits <- logical(0); widths <- numeric(0)

  for (f in folds) {
    keep <- setdiff(seq_len(n), f)
    if (length(keep) < 2L) next
    p <- .estimate_ou_kalman_params(d[keep], s[keep], n_fit_min = n_fit_min,
                                    scale = scale)
    if (p$tier == "degenerate" || p$gamma <= 0) next
    grid <- seq(min(d), max(d), by = "day")       # DAILY grid (phi assumes dt=1)
    mod  <- .build_kalman_model(grid, d[keep], s[keep], p$theta, p$gamma)
    if (is.null(mod)) next
    sm   <- .kalman_smooth(mod)
    r_obs <- max(p$gamma, 1) * 1e-6               # match the model's anchor noise
    for (i in f) {
      gi  <- match(d[i], grid)
      sdi <- sqrt(sm$var[gi] + r_obs)
      lo  <- sm$mean[gi] - z * sdi; hi <- sm$mean[gi] + z * sdi
      hits   <- c(hits, s[i] >= lo & s[i] <= hi)
      widths <- c(widths, hi - lo)
    }
  }

  if (length(hits) == 0L)
    return(list(coverage = NA_real_, mean_width = NA_real_, n = n, n_held = 0L))
  list(coverage = mean(hits), mean_width = mean(widths), n = n,
       n_held = length(hits))
}

#' Per-analyte and pooled LOO coverage over a fitted target model
#'
#' @param target_model A `target_model` (uses each model's `$anchors$S` /
#'   `$d_anchors$S` and `$tier`).
#' @param interval,block,n_fit_min,scale Passed to [.loo_coverage_series()].
#' @return A tibble: one row per analyte (`analyte, tier, n, coverage,
#'   mean_width`) plus a final `(pooled)` row (n-weighted coverage, mean width).
#' @keywords internal
.loo_anchor_coverage <- function(target_model, interval = 0.9, block = 1L,
                                 n_fit_min = 8L, scale = 1) {
  models <- target_model$models %||% list()
  rows <- lapply(names(models), function(nm) {
    m <- models[[nm]]
    anch <- if (!is.null(m$d_anchors) && nrow(m$d_anchors) >= 2L) m$d_anchors
            else m$anchors
    if (is.null(anch) || is.null(anch$S))
      return(tibble::tibble(analyte = nm, tier = m$tier %||% NA_character_,
                            n = 0L, coverage = NA_real_, mean_width = NA_real_))
    cov <- .loo_coverage_series(anch$date, anch$S, interval = interval,
                                block = block, n_fit_min = n_fit_min,
                                scale = scale)
    tibble::tibble(analyte = nm, tier = m$tier %||% NA_character_,
                   n = cov$n, coverage = cov$coverage,
                   mean_width = cov$mean_width)
  })
  res <- dplyr::bind_rows(rows)

  ok <- !is.na(res$coverage)
  pooled <- tibble::tibble(
    analyte = "(pooled)", tier = NA_character_,
    n = sum(res$n[ok]),
    coverage   = if (any(ok)) stats::weighted.mean(res$coverage[ok], res$n[ok]) else NA_real_,
    mean_width = if (any(ok)) stats::weighted.mean(res$mean_width[ok], res$n[ok]) else NA_real_)
  dplyr::bind_rows(res, pooled)
}
