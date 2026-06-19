## reference_model.R
##
## Contemporaneous (temporal) ARA reference model.  Issue #9.
##
## The package assumes reference and target share the same catchment hydrology;
## a single daily hydrology series (rainfall → API, or stage/discharge →
## antecedent mean) drives both the tier-1 matching gate and the tier-2 GAM.
##
## Public surface
## ──────────────
##   fit_reference_model()   fit per-analyte temporal models; returns reference_model
##   get_silo_rainfall()     fetch SILO daily rainfall (in watertemp.R)
##   ara_summary()           accessor for per-cell ARA diagnostics on add_amspaf() output
##
## Internal
## ────────
##   .compute_api()                  API — recursive linear reservoir (#49)
##   .compute_antecedent_mean()      rolling mean for stage/discharge
##   .compute_hydro_features()       dispatcher; returns (date, hydro_short, hydro_long)
##   .build_analyte_model_df()       assemble (y, doy, hydro_short, hydro_long) for one analyte
##   .fit_ref_gam()                  fit mgcv::gam with gamm fallback
##   .aic_for_tau()                  AIC for a (tau_short, tau_long) candidate pair
##   .select_api_tau()               continuous tau selection by profiled GAM AIC (#49)
##   .resolve_ref_norm()             three-tier resolver dispatcher
##   .resolve_ref_norm_instant()     per-sample (instantaneous) resolver
##   .resolve_ref_norm_chronic()     per-sample (chronic window-integrated) resolver
##   .predict_ref_at_date()          predict value_norm from a fitted analyte model


## ============================================================================
## Internal hydrology helpers
## ============================================================================

## Parsimonious fallback recession constants (days) for the rainfall reservoir:
## a ~1-week fast store and a ~2-month slow store.  Used when tau selection is
## disabled (`auto_select = FALSE`), as the overfitting-gate baseline in
## [.select_api_tau()], and for analytes that fall to the static tier.
.REF_TAU_DEFAULT_SHORT <- 7
.REF_TAU_DEFAULT_LONG <- 60

#' Antecedent Precipitation Index (API) — exact recursive linear reservoir
#'
#' Implements the exact discrete solution of the single linear store
#' `dS/dt = -S/tau + P(t)`:
#'
#' ```
#' S_t = k^{dt} * S_{t-1} + P_t,   k = exp(-1 / tau)
#' ```
#'
#' This is the convergent infinite-horizon form of the classical Antecedent
#' Precipitation Index (Kohler & Linsley 1951), whose theoretical basis is the
#' Maillet (1905) exponential recession.  Unlike a windowed sum it has no
#' truncation horizon: all prior rainfall contributes with exponentially
#' decaying weight `k^{dt}` set by the actual day gap `dt`, so irregular
#' spacing is handled exactly.
#'
#' Algorithm: lay the rainfall onto a complete daily grid (`NA` and gap days
#' contribute 0 input, so the store simply decays across them).  On a daily
#' grid every step has `dt = 1`, so the reservoir reduces to the first-order
#' linear recursion `S_t = k * S_{t-1} + P_t`, evaluated in C by
#' [stats::filter()] (`method = "recursive"`).  This is an exact reorganisation
#' of the per-event form `S_t = k^{dt} * S_{t-1} + P_t` (`k^{dt}` is just `dt`
#' unit steps with zero input between events), with none of the per-target
#' re-summation the old windowed form required.  Target dates index straight
#' into the grid; any falling before the first hydro day return 0 (empty
#' reservoir).
#'
#' @param hydro_values Numeric vector of daily rainfall values (`NA` is 0 input).
#' @param hydro_dates Date vector matching `hydro_values`.
#' @param target_dates Date vector; API is evaluated at each.
#' @param tau Positive numeric; reservoir memory constant (days).  The weight
#'   for a gap of `dt` days is `exp(-dt / tau)`.
#' @return Numeric vector, same length as `target_dates`.
#' @references
#' Maillet, E. (1905) *Essais d'hydraulique souterraine et fluviale.* Hermann,
#' Paris.  Kohler, M.A. & Linsley, R.K. (1951) *Predicting the runoff from
#' storm rainfall.* US Weather Bureau Research Paper 34.
#' @keywords internal
.compute_api <- function(hydro_values, hydro_dates, target_dates, tau) {
  stopifnot(
    "`tau` must be a single positive number" =
      length(tau) == 1L && is.numeric(tau) && is.finite(tau) && tau > 0
  )
  k <- exp(-1 / tau)

  hdates <- as.Date(hydro_dates)
  hvals <- hydro_values
  hvals[is.na(hvals)] <- 0 # NA rainfall is zero input; the store still drains

  n <- length(hdates)
  if (n == 0L) {
    return(rep(0, length(target_dates)))
  }

  # Complete daily grid from the first hydro day to the last day we must report
  # (targets may fall after the final hydro observation; the grid decays into
  # that tail with zero input).
  grid_start <- min(hdates)
  grid_end <- max(c(hdates, as.Date(target_dates)))
  span <- as.integer(grid_end - grid_start) + 1L

  # Accumulate rainfall onto the grid; duplicate-date rows sum, matching the
  # sequential dt = 0 additions of the per-event recursion.
  rain <- numeric(span)
  pos <- as.integer(hdates - grid_start) + 1L
  agg <- tapply(hvals, pos, sum)
  rain[as.integer(names(agg))] <- as.numeric(agg)

  # First-order recursive filter in C: state[t] = rain[t] + k * state[t-1],
  # seeded from an empty reservoir.
  state <- as.numeric(stats::filter(rain, filter = k, method = "recursive"))

  # Targets index directly into the grid; those before grid_start stay 0.
  out <- numeric(length(target_dates))
  tpos <- as.integer(as.Date(target_dates) - grid_start) + 1L
  inside <- tpos >= 1L & tpos <= span
  out[inside] <- state[tpos[inside]]
  out
}

#' Antecedent mean for stage or discharge
#'
#' Rolling mean over `[t - window_days, t]`.  Returns `NA` for target dates
#' with no hydro record in the window.
#'
#' @param hydro_values Numeric vector of daily stage/discharge values.
#' @param hydro_dates Date vector matching `hydro_values`.
#' @param target_dates Date vector; mean is evaluated at each.
#' @param window_days Integer memory length (days).
#' @return Numeric vector, same length as `target_dates`.
#' @keywords internal
.compute_antecedent_mean <- function(hydro_values, hydro_dates, target_dates, window_days) {
  vapply(target_dates, function(td) {
    mask <- hydro_dates >= (td - window_days) & hydro_dates <= td
    if (!any(mask)) {
      return(NA_real_)
    }
    mean(hydro_values[mask], na.rm = TRUE)
  }, numeric(1L))
}

#' Compute hydrology features (hydro_short, hydro_long) at target dates
#'
#' Dispatches to [.compute_api()] for rainfall or [.compute_antecedent_mean()]
#' for stage/discharge, at both the short and long tau values.
#'
#' @param hydro Data frame with columns `date` (Date) and `value` (numeric).
#' @param target_dates Date vector.
#' @param tau_short,tau_long Memory constants (days).
#' @param hydro_type `"rainfall"`, `"stage"`, or `"discharge"`.
#' @return Tibble `(date, hydro_short, hydro_long)`.
#' @keywords internal
.compute_hydro_features <- function(hydro, target_dates,
                                    tau_short, tau_long,
                                    hydro_type = "rainfall") {
  hdates <- as.Date(hydro$date)
  hvals <- hydro$value

  if (hydro_type == "rainfall") {
    short <- .compute_api(hvals, hdates, target_dates, tau_short)
    long <- .compute_api(hvals, hdates, target_dates, tau_long)
  } else {
    short <- .compute_antecedent_mean(hvals, hdates, target_dates, tau_short)
    long <- .compute_antecedent_mean(hvals, hdates, target_dates, tau_long)
  }

  tibble::tibble(date = target_dates, hydro_short = short, hydro_long = long)
}


## ============================================================================
## Internal GAM helpers
## ============================================================================

#' Build the model data frame for one analyte
#'
#' Computes hydro features at the observation dates and assembles
#' `(y, doy, hydro_short, hydro_long)` where `y = log(value_norm + eps)`.
#'
#' @param obs Tibble `(date, value_norm)` for one analyte.
#' @param hydro Daily hydro series.
#' @param tau_short,tau_long Reservoir memory constants (days).
#' @param hydro_type Hydro type string.
#' @param eps Guard for log transform.
#' @return Tibble `(y, doy, hydro_short, hydro_long)` sorted by date.
#' @keywords internal
.build_analyte_model_df <- function(obs, hydro, tau_short, tau_long,
                                    hydro_type, eps = 1e-9) {
  obs <- obs[order(obs$date), ]
  feats <- .compute_hydro_features(hydro, obs$date, tau_short, tau_long, hydro_type)

  tibble::tibble(
    date        = obs$date,
    y           = log(pmax(obs$value_norm, eps)),
    doy         = as.integer(format(obs$date, "%j")),
    hydro_short = feats$hydro_short,
    hydro_long  = feats$hydro_long
  )
}

#' Fit a GAM on reference observations for one analyte
#'
#' Tries `mgcv::gam()` with cyclic spline on day-of-year and thin-plate
#' splines on the two hydro features.  Returns the fitted `gam` object, or
#' `NULL` if fitting fails.
#'
#' @param df_model Tibble `(y, doy, hydro_short, hydro_long)`.
#' @param eps Guard already applied to `y` upstream (unused here, for docs).
#' @return A `gam` object or `NULL`.
#' @keywords internal
.fit_ref_gam <- function(df_model, eps = 1e-9) { # nolint: unused arg
  n <- nrow(df_model)
  if (n < 10L) {
    return(NULL)
  }

  k_doy <- min(6L, n - 2L)
  k_h <- min(4L, n - 2L)

  tryCatch(
    mgcv::gam(
      y ~ s(doy, bs = "cc", k = k_doy) +
        s(hydro_short, k = k_h) +
        s(hydro_long, k = k_h),
      data = df_model,
      knots = list(doy = c(1L, 365L)),
      method = "REML"
    ),
    error = function(e) NULL,
    warning = function(w) {
      tryCatch(
        mgcv::gam(
          y ~ s(doy, bs = "cc", k = k_doy) +
            s(hydro_short, k = k_h) +
            s(hydro_long, k = k_h),
          data = df_model,
          knots = list(doy = c(1L, 365L)),
          method = "REML"
        ),
        error = function(e) NULL
      )
    }
  )
}

#' AIC for a single (tau_short, tau_long) candidate
#'
#' Returns `Inf` if fitting fails.
#' @keywords internal
.aic_for_tau <- function(obs, hydro, tau_short, tau_long, hydro_type, eps) {
  df_m <- tryCatch(
    .build_analyte_model_df(obs, hydro, tau_short, tau_long, hydro_type, eps),
    error = function(e) NULL
  )
  if (is.null(df_m)) {
    return(Inf)
  }
  fit <- .fit_ref_gam(df_m, eps)
  if (is.null(fit)) {
    return(Inf)
  }
  stats::AIC(fit)
}

#' Minimise an AIC objective over a (tau_short, tau_long) pair
#'
#' Deterministic golden-section coordinate descent ([stats::optimize()]) over
#' the two reservoir recession constants, holding one fixed while the other is
#' optimised (two passes).  The fast/slow separation `tau_long >= 1.5*tau_short`
#' is imposed by capping tau_short's upper bound at `tau_long/1.5` and lifting
#' tau_long's lower bound to `1.5*tau_short`.  A degenerate bound `lo == hi`
#' pins that store.  Shared by the reference ([.select_api_tau()]) and impact
#' ([.fit_impact_response()]) hydrology fits.  No RNG.
#'
#' @param aic_fn Function `(tau_short, tau_long) -> AIC` (`Inf` on fit failure).
#' @param tau_bounds_short,tau_bounds_long Length-2 numeric `c(lo, hi)` (days).
#' @param default_short,default_long Starting tau (days).
#' @return List `(tau_short, tau_long)`.
#' @keywords internal
.optimise_tau_pair <- function(aic_fn, tau_bounds_short, tau_bounds_long,
                               default_short = .REF_TAU_DEFAULT_SHORT,
                               default_long = .REF_TAU_DEFAULT_LONG) {
  opt1 <- function(lo, hi, f) {
    if (lo >= hi) {
      return(lo)
    }
    stats::optimize(f, lower = lo, upper = hi)$minimum
  }
  ts <- min(max(default_short, tau_bounds_short[1L]), tau_bounds_short[2L])
  tl <- min(max(default_long, tau_bounds_long[1L]), tau_bounds_long[2L])
  for (pass in seq_len(2L)) {
    hi_s <- min(tau_bounds_short[2L], tl / 1.5)
    ts <- opt1(tau_bounds_short[1L], hi_s, function(x) aic_fn(x, tl))
    lo_l <- max(tau_bounds_long[1L], 1.5 * ts)
    tl <- opt1(lo_l, tau_bounds_long[2L], function(x) aic_fn(ts, x))
  }
  list(tau_short = ts, tau_long = tl)
}

#' Select reservoir recession constants (tau) per analyte by profiled AIC
#'
#' Chooses the short- and long-store recession constants `tau_short`,
#' `tau_long` (days) for the rainfall reservoir by minimising the fitted GAM's
#' AIC.  Given tau the GAM is conditionally linear, so tau is selected by an
#' *outer* 1-D optimisation of the profiled AIC — profile likelihood for a
#' nonlinear-in-parameter feature (Wood 2017).  The two stores are optimised by
#' deterministic golden-section search ([stats::optimize()]) in a short
#' coordinate descent (tau_short with tau_long held, then tau_long with
#' tau_short held); no RNG is used.
#'
#' Guards:
#' * **Separation** `tau_long >= 1.5 * tau_short` keeps the fast and slow stores
#'   distinct (otherwise the two smooths become collinear).
#' * **Overfitting gate** — the optimised tau is adopted only when it improves
#'   AIC by `>= 2` over a parsimonious default `(default_short, default_long)`;
#'   otherwise the default is returned.  ΔAIC ≥ 2 is the conventional threshold
#'   for a distinguishable model (Burnham & Anderson 2002), so a flat or
#'   uninformative AIC surface falls back to the simpler model rather than
#'   chasing noise.
#'
#' @param obs Tibble `(date, value_norm)` for one analyte.
#' @param hydro Daily hydro series.
#' @param hydro_type Hydro type string.
#' @param tau_bounds_short,tau_bounds_long Length-2 numeric `c(lo, hi)` search
#'   ranges (days).  A degenerate `lo == hi` fixes that store's tau.
#' @param default_short,default_long Parsimonious fallback tau (days).
#' @param eps Log guard.
#' @return List `(tau_short, tau_long, best_aic, null_aic)`.
#' @references
#' Wood, S.N. (2017) *Generalized Additive Models: An Introduction with R*, 2nd
#' ed. CRC Press.  Burnham, K.P. & Anderson, D.R. (2002) *Model Selection and
#' Multimodel Inference*, 2nd ed. Springer.
#' @keywords internal
.select_api_tau <- function(obs, hydro, hydro_type,
                            tau_bounds_short, tau_bounds_long,
                            default_short = 7, default_long = 60,
                            eps = 1e-9) {
  aic <- function(ts, tl) .aic_for_tau(obs, hydro, ts, tl, hydro_type, eps)

  # Null model AIC (intercept only), for the downstream tier decision.
  df0 <- tryCatch(
    .build_analyte_model_df(
      obs, hydro, default_short, default_long, hydro_type, eps
    ),
    error = function(e) NULL
  )
  null_aic <- if (!is.null(df0)) {
    fit0 <- tryCatch(mgcv::gam(y ~ 1, data = df0, method = "REML"),
      error = function(e) NULL
    )
    if (!is.null(fit0)) stats::AIC(fit0) else Inf
  } else {
    Inf
  }

  tp <- .optimise_tau_pair(
    aic, tau_bounds_short, tau_bounds_long, default_short, default_long
  )
  ts <- tp$tau_short
  tl <- tp$tau_long

  best_aic <- aic(ts, tl)
  default_aic <- aic(default_short, default_long)

  # Adopt the optimised tau only if it beats the parsimonious default by ≥ 2
  # AIC; otherwise keep the default (guards against overfitting a flat surface).
  if (is.finite(best_aic) && best_aic <= default_aic - 2) {
    list(
      tau_short = ts, tau_long = tl,
      best_aic = best_aic, null_aic = null_aic
    )
  } else {
    list(
      tau_short = default_short, tau_long = default_long,
      best_aic = default_aic, null_aic = null_aic
    )
  }
}


## ============================================================================
## fit_reference_model()
## ============================================================================

#' Fit a temporal reference model for contemporaneous ARA background subtraction
#'
#' Fits per-analyte temporal models on reference-site chemistry so that
#' [add_amspaf()] can subtract a *contemporaneous* background (what the
#' reference site would have shown at the same moment as the target sample)
#' rather than a static time-average.
#'
#' **Package invariant:** reference and target are assumed to share the same
#' catchment.  The hydrology series supplied here (or fetched automatically
#' from SILO) applies to both sites.  This invariant should hold for
#' near-field leachate assessments where a headwater reference is paired with
#' a downstream target in the same sub-catchment.
#'
#' **Three-tier resolver** (evaluated per analyte × target date):
#' \enumerate{
#'   \item **Tier 1 — direct match.** A reference observation within
#'     `±match_window_days` days whose event-API is within `match_hydro_tol`
#'     (default: `0.5 × IQR` of the reference API series).  This gate rejects
#'     time-close but hydrologically-mismatched grabs (e.g. a dry-weather
#'     reference next to a wet-event target).
#'   \item **Tier 2 — GAM prediction.** A per-analyte `mgcv::gam` with
#'     `s(doy, bs="cc") + s(hydro_short) + s(hydro_long)`.  The reservoir
#'     recession constants `tau_short`, `tau_long` are selected per analyte by
#'     profiled AIC over the `api_tau_bounds_short`/`api_tau_bounds_long`
#'     ranges (continuous, with a `tau_long >= 1.5*tau_short` separation and a
#'     ΔAIC ≥ 2 adoption gate over a parsimonious default).  The analyte falls
#'     back to tier 3 if the best model has higher AIC than the null
#'     (intercept-only) model, or if fewer than `min_obs_model` detected
#'     observations are available.
#'   \item **Tier 3 — static fallback.** Geometric mean of all normalised
#'     reference observations (identical to [prepare_reference()] with
#'     `summary = "geom_mean"`).
#' }
#'
#' **Hydrology input** — supply exactly one of:
#' \itemize{
#'   \item `hydro`: a data frame with columns `date` (Date) and `value`
#'     (numeric), plus a `type` column OR supply `hydro_type =` separately.
#'     Supported types: `"rainfall"` (daily mm; → API), `"stage"` (m;
#'     → antecedent mean), `"discharge"` (m³/s; → antecedent mean).
#'   \item `latitude` + `longitude`: if `hydro = NULL`, SILO daily rainfall is
#'     fetched automatically via [get_silo_rainfall()].  Requires the
#'     **weatherOz** package and a valid SILO API key.
#' }
#'
#' @param reference Long-format reference chemistry data frame.  Required
#'   columns: `sample_id`, `site_id`, `datetime`, `analyte`, `value`,
#'   `detected`.  Toxicant concentrations must be in µg/L; supply them via a
#'   `units.analyte` column or via `conc_units`.
#' @param hydro Daily hydrology data frame (`date`, `value`), or `NULL` to
#'   fetch SILO rainfall automatically.  See *Hydrology input* above.
#' @param hydro_type Character; one of `"rainfall"`, `"stage"`, `"discharge"`.
#'   Ignored if `hydro` already has a `type` column.  Default `"rainfall"`.
#' @param latitude,longitude Decimal-degree coordinates; required only when
#'   `hydro = NULL` (SILO auto-fetch path).
#' @param conc_units Unit string (e.g. `"mg/L"`, `"ug/L"`) for reference
#'   chemistry when no `units.analyte` column is present.
#' @param analyte_metadata Analyte metadata, or `NULL` for the bundled CSV.
#' @param imputation_model Optional `imputation_model` from
#'   [fit_imputation_model()] (fit on the reference site's own chemistry).  When
#'   supplied, missing analytes are imputed in raw concentration space *before*
#'   the per-analyte models are fitted, so a well-sampled analyte lifts a
#'   sparsely-sampled one into a richer spread of hydrological regimes.  Imputed
#'   rows (`detected = TRUE`) are used as model anchors alongside measured
#'   observations.  Requires **brms**.  Default `NULL` (measured observations
#'   only).
#' @param match_window_days Integer; tier-1 time tolerance in days (default 5).
#' @param match_hydro_tol Numeric; tier-1 hydro tolerance (default `NULL` →
#'   `0.5 × IQR` of the reference event-API series).
#' @param api_tau_bounds_short Length-2 numeric `c(lo, hi)` search range (days)
#'   for the short-store recession constant `tau_short` (default `c(1, 30)`).
#'   A degenerate `c(x, x)` fixes `tau_short = x`.
#' @param api_tau_bounds_long Length-2 numeric `c(lo, hi)` search range (days)
#'   for the long-store recession constant `tau_long` (default `c(20, 365)`).
#' @param auto_select Logical; if `TRUE` (default), select `tau_short`,
#'   `tau_long` per analyte by profiled AIC.  If `FALSE`, use the parsimonious
#'   defaults (`tau_short = 7`, `tau_long = 60`) for all analytes.
#' @param min_obs_model Integer; minimum detected observations required to
#'   attempt tier-2 modelling (default `20L`).  Analytes below this fall
#'   directly to tier 3.
#' @param summary Static-fallback summary statistic (tier 3): one of
#'   `"geom_mean"` (default), `"median"`, `"arith_mean"`, `"p80"`, `"p90"`,
#'   `"p95"`.
#' @param silo_start,silo_end Start/end dates for SILO auto-fetch.  Default
#'   `NULL`: derived from the reference chemistry date range, padded on the left
#'   by `5 × max(api_tau_bounds_long)` days so the recursive reservoir has
#'   enough burn-in (about 5 tau) to converge before the first observation.
#' @param silo_api_key Passed to [get_silo_rainfall()] when `hydro = NULL`.
#' @param eps Small positive guard for log transform (default `1e-9`).
#'
#' @return An object of class `reference_model`:
#'   \describe{
#'     \item{`$models`}{Named list (one per analyte) carrying `gamm_fit`,
#'       `tau_short`, `tau_long`, `best_aic`, `null_aic`, `tier`,
#'       `n_obs`, `static_ref`, and `obs` (normalised observations with
#'       hydro features — used for tier-1 matching).}
#'     \item{`$hydro`}{Daily hydro series used for index computation.}
#'     \item{`$hydro_type`}{`"rainfall"`, `"stage"`, or `"discharge"`.}
#'     \item{`$match_window_days`}{Tier-1 time window.}
#'     \item{`$match_hydro_tol`}{Tier-1 hydro tolerance (computed or supplied).}
#'     \item{`$static_ref`}{Tibble `(analyte, ref_norm)` — tier-3 fallback for
#'       all analytes.}
#'     \item{`$fit_date`}{Date the model was fitted.}
#'   }
#'
#' @seealso [add_amspaf()], [get_silo_rainfall()], [prepare_reference()],
#'   [ara_summary()]
#'
#' @examples
#' \dontrun{
#' # Fit using SILO auto-fetch
#' ref <- subset(leachate_demo(), site_id == "reference")
#' ref_model <- fit_reference_model(
#'   reference  = ref,
#'   latitude   = -33.87,
#'   longitude  = 151.21,
#'   conc_units = "ug/L"
#' )
#' ref_model
#'
#' # Or supply your own gauge record
#' ref_model2 <- fit_reference_model(
#'   reference  = ref,
#'   hydro      = my_stage_df, # data.frame(date, value)
#'   hydro_type = "stage",
#'   conc_units = "ug/L"
#' )
#' }
#' @export
fit_reference_model <- function(
  reference,
  hydro = NULL,
  hydro_type = "rainfall",
  latitude = NULL,
  longitude = NULL,
  conc_units = NULL,
  analyte_metadata = NULL,
  imputation_model = NULL,
  match_window_days = 5L,
  match_hydro_tol = NULL,
  api_tau_bounds_short = c(1, 30),
  api_tau_bounds_long = c(20, 365),
  auto_select = TRUE,
  min_obs_model = 20L,
  summary = "geom_mean",
  silo_start = NULL,
  silo_end = NULL,
  silo_api_key = NULL,
  eps = 1e-9
) {
  ## ── Input validation ───────────────────────────────────────────────────────
  checkmate::assert_data_frame(reference)
  checkmate::assert_names(names(reference),
    must.include = c("sample_id", "datetime", "analyte", "value", "detected")
  )
  checkmate::assert_int(match_window_days, lower = 0L)
  checkmate::assert_flag(auto_select)
  checkmate::assert_int(min_obs_model, lower = 5L)
  summary <- match.arg(
    summary,
    c("geom_mean", "median", "arith_mean", "p80", "p90", "p95")
  )
  checkmate::assert_number(eps, lower = 0)
  checkmate::assert_numeric(
    api_tau_bounds_short,
    len = 2L, lower = 0, any.missing = FALSE, sorted = TRUE
  )
  checkmate::assert_numeric(
    api_tau_bounds_long,
    len = 2L, lower = 0, any.missing = FALSE, sorted = TRUE
  )
  if (!is.null(imputation_model) &&
    !inherits(imputation_model, "imputation_model")) {
    cli::cli_abort(
      "{.arg imputation_model} must be an {.cls imputation_model} from \\
       {.fn fit_imputation_model}, or {.val NULL}."
    )
  }

  ## ── Hydrology ──────────────────────────────────────────────────────────────
  if (is.null(hydro)) {
    if (is.null(latitude) || is.null(longitude)) {
      cli::cli_abort(c(
        "Supply either {.arg hydro} (a daily hydrology data frame) or \\
         both {.arg latitude} and {.arg longitude} to fetch SILO rainfall.",
        "i" = "See {.fn get_silo_rainfall} for the expected data format."
      ))
    }

    ref_dates <- as.Date(reference$datetime)
    s_start <- if (!is.null(silo_start)) {
      as.Date(silo_start)
    } else {
      # ~5\u03c4 burn-in so the recursive reservoir converges before obs start.
      min(ref_dates, na.rm = TRUE) - 5 * max(api_tau_bounds_long)
    }
    s_end <- if (!is.null(silo_end)) {
      as.Date(silo_end)
    } else {
      max(ref_dates, na.rm = TRUE)
    }

    cli::cli_inform("Fetching SILO rainfall for ({latitude}, {longitude}) \\
                    {s_start}\u2013{s_end}\u2026")
    rain_df <- get_silo_rainfall(
      latitude   = latitude,
      longitude  = longitude,
      start_date = s_start,
      end_date   = s_end,
      api_key    = silo_api_key
    )
    hydro <- dplyr::rename(rain_df, value = "rainfall_mm")
    hydro_type <- "rainfall"
  } else {
    checkmate::assert_data_frame(hydro)
    checkmate::assert_names(names(hydro), must.include = c("date", "value"))
    if ("type" %in% names(hydro)) {
      hydro_type <- unique(hydro$type)[1L]
    }
    hydro_type <- match.arg(hydro_type, c("rainfall", "stage", "discharge"))
    hydro <- dplyr::select(hydro, date = "date", value = "value")
  }

  hydro <- dplyr::mutate(hydro,
    date = as.Date(.data$date),
    value = as.numeric(.data$value)
  )
  hydro <- hydro[order(hydro$date), ]

  ## ── Unit conversion + normalisation ───────────────────────────────────────
  meta <- .load_analyte_metadata(analyte_metadata)
  ssd_analytes <- meta$analyte[
    !is.na(meta$ssd_available) & meta$ssd_available == TRUE &
      !meta$analyte %in% .AMSPAF_EXCLUDED_ANALYTES
  ]
  reference <- .convert_df_tox_to_ugL(reference, ssd_analytes, conc_units, "reference")

  ## ── Optional impute-first ──────────────────────────────────────────────────
  ## Complete missing analytes in raw concentration space before modelling, so a
  ## well-sampled analyte lifts a sparsely-sampled one into a richer spread of
  ## hydrological regimes (e.g. 100 Zn obs → 100 Cu anchors via the Cu–Zn
  ## relationship).  Imputed metal rows are marked detected = TRUE and flow into
  ## the per-analyte GAMs exactly like measured observations.
  if (!is.null(imputation_model) &&
    length(imputation_model$groups %||% list()) > 0L) {
    if (!"site_id" %in% names(reference)) {
      reference$site_id <- "reference"
    }
    reference <- impute_chemistry(reference, imputation_model)
  }

  # Build working frame: BDL rows contribute 0 to the static summary but keep
  # their DL value for normalisation (so value_norm is not zero-inflated before
  # the GAM sees it).  We model only detected observations in tier 2.
  ref_detected <- dplyr::filter(reference, .data$detected)
  ref_for_norm <- dplyr::mutate(
    reference,
    .date = as.Date(.data$datetime)
  )
  # BDL → 0 for the static summary path (mirrors prepare_reference())
  ref_for_static <- dplyr::mutate(
    ref_for_norm,
    value = dplyr::if_else(.data$detected, .data$value, 0)
  )

  ## Normalise (detected obs keep their value; BDL rows get value_norm based
  ## on DL for GAM fitting, but the GAM only uses detected rows anyway)
  norm_all <- .normalise_ref_observations(ref_for_norm, reference, meta)

  ## Detected-only normalised observations for GAM fitting
  norm_det <- dplyr::filter(
    norm_all, .data$detected, !is.na(.data$value_norm),
    .data$value_norm > 0, .data$analyte %in% ssd_analytes
  )

  ## Static fallback: geometric mean of normalised values (detected+BDL=0)
  norm_static <- .normalise_ref_observations(ref_for_static, reference, meta)
  static_ref <- norm_static |>
    dplyr::filter(.data$analyte %in% ssd_analytes) |>
    dplyr::group_by(.data$analyte) |>
    dplyr::summarise(
      ref_norm = .ref_summary(.data$value_norm, summary, eps),
      n_obs    = sum(!is.na(.data$value_norm)),
      .groups  = "drop"
    ) |>
    dplyr::filter(.data$n_obs > 0L)

  ## ── Per-analyte model fitting ──────────────────────────────────────────────
  target_analytes <- intersect(
    unique(norm_det$analyte),
    ssd_analytes
  )

  models <- vector("list", length(target_analytes))
  names(models) <- target_analytes

  for (nm in target_analytes) {
    obs <- norm_det |>
      dplyr::filter(.data$analyte == .env$nm) |>
      dplyr::select(date = ".date", "value_norm") |>
      dplyr::distinct(.data$date, .keep_all = TRUE)

    n_obs <- nrow(obs)
    static_val <- static_ref$ref_norm[static_ref$analyte == nm]
    static_val <- if (length(static_val)) static_val[1L] else NA_real_

    # Add hydro features at observation dates (longest-memory tau in each store)
    # for tier-1 matching; stored regardless of tier
    max_short <- max(api_tau_bounds_short)
    max_long <- max(api_tau_bounds_long)
    obs_feats <- .compute_hydro_features(
      hydro, obs$date, max_short, max_long, hydro_type
    )
    obs_with_feats <- dplyr::bind_cols(obs, dplyr::select(obs_feats, -"date"))

    if (n_obs < min_obs_model) {
      models[[nm]] <- list(
        gamm_fit   = NULL,
        tau_short  = .REF_TAU_DEFAULT_SHORT,
        tau_long   = .REF_TAU_DEFAULT_LONG,
        best_aic   = NA_real_,
        null_aic   = NA_real_,
        tier       = "static",
        n_obs      = n_obs,
        static_ref = static_val,
        obs        = obs_with_feats
      )
      next
    }

    if (auto_select) {
      sel <- .select_api_tau(
        obs, hydro, hydro_type, api_tau_bounds_short, api_tau_bounds_long,
        default_short = .REF_TAU_DEFAULT_SHORT,
        default_long = .REF_TAU_DEFAULT_LONG, eps = eps
      )
    } else {
      df0 <- .build_analyte_model_df(
        obs, hydro, .REF_TAU_DEFAULT_SHORT, .REF_TAU_DEFAULT_LONG,
        hydro_type, eps
      )
      fit0 <- tryCatch(
        mgcv::gam(y ~ 1, data = df0, method = "REML"),
        error = function(e) NULL
      )
      sel <- list(
        tau_short = .REF_TAU_DEFAULT_SHORT,
        tau_long  = .REF_TAU_DEFAULT_LONG,
        best_aic  = if (!is.null(fit0)) stats::AIC(fit0) else Inf,
        null_aic  = if (!is.null(fit0)) stats::AIC(fit0) else Inf
      )
    }

    df_m <- .build_analyte_model_df(
      obs, hydro, sel$tau_short, sel$tau_long, hydro_type, eps
    )
    gam_fit <- .fit_ref_gam(df_m, eps)

    # Model is useful only if it beats the null
    tier <- if (
      is.null(gam_fit) ||
        is.infinite(sel$best_aic) ||
        sel$best_aic >= sel$null_aic
    ) {
      "static"
    } else {
      "model"
    }

    models[[nm]] <- list(
      gamm_fit   = gam_fit,
      tau_short  = sel$tau_short,
      tau_long   = sel$tau_long,
      best_aic   = sel$best_aic,
      null_aic   = sel$null_aic,
      tier       = tier,
      n_obs      = n_obs,
      static_ref = static_val,
      obs        = obs_with_feats
    )
  }

  ## ── Tier-1 hydro tolerance ────────────────────────────────────────────────
  # Compute using the shortest-memory short tau on the reference observations
  if (is.null(match_hydro_tol)) {
    all_obs_dates <- as.Date(reference$datetime)
    event_api <- .compute_api(
      hydro$value, as.Date(hydro$date),
      all_obs_dates, api_tau_bounds_short[1L]
    )
    match_hydro_tol <- 0.5 * stats::IQR(event_api, na.rm = TRUE)
    if (match_hydro_tol < 1e-6) match_hydro_tol <- Inf # no variability → no gate
  }

  structure(
    list(
      models            = models,
      hydro             = hydro,
      hydro_type        = hydro_type,
      match_window_days = as.integer(match_window_days),
      match_hydro_tol   = match_hydro_tol,
      static_ref        = static_ref,
      fit_date          = Sys.Date(),
      summary           = summary
    ),
    class = "reference_model"
  )
}

#' @export
print.reference_model <- function(x, ...) {
  n_analytes <- length(x$models)
  n_model <- sum(vapply(x$models, function(m) m$tier == "model", logical(1L)))
  n_static <- n_analytes - n_model

  hydro_range <- if (!is.null(x$hydro) && nrow(x$hydro) > 0L) {
    sprintf("%s \u2013 %s", min(x$hydro$date), max(x$hydro$date))
  } else {
    "unknown"
  }

  cat(sprintf(
    "<reference_model>  fitted %s | %d analyte%s | hydro: %s (%s)\n",
    x$fit_date, n_analytes,
    if (n_analytes == 1L) "" else "s",
    x$hydro_type, hydro_range
  ))
  cat(sprintf(
    "  tier-1:  \u00b1%d d + API \u00b1%.2g\n",
    x$match_window_days, x$match_hydro_tol
  ))

  if (n_model > 0L) {
    model_nms <- names(Filter(function(m) m$tier == "model", x$models))
    detail <- vapply(model_nms, function(nm) {
      m <- x$models[[nm]]
      sprintf(
        "%s (\u03c4=%.0f/%.0fd, n=%d)", nm, m$tau_short, m$tau_long, m$n_obs
      )
    }, character(1L))
    cat(sprintf("  tier-2 (%d):  %s\n", n_model, paste(detail, collapse = ", ")))
  }
  if (n_static > 0L) {
    static_nms <- names(Filter(function(m) m$tier == "static", x$models))
    cat(sprintf("  tier-3 (%d):  %s\n", n_static, paste(static_nms, collapse = ", ")))
  }
  invisible(x)
}


## ============================================================================
## Internal resolver
## ============================================================================

#' Predict reference value_norm at a single date for one analyte model
#'
#' @param m Analyte model object from `reference_model$models[[nm]]`.
#' @param target_date Date scalar.
#' @param hydro Daily hydro series.
#' @param hydro_type Hydro type string.
#' @param eps Log guard.
#' @return List `(ref_norm, ref_tier)`.
#' @keywords internal
.predict_ref_at_date <- function(m, target_date, hydro, hydro_type, eps = 1e-9) {
  if (m$tier == "static" || is.null(m$gamm_fit)) {
    return(list(ref_norm = m$static_ref %||% NA_real_, ref_tier = "static"))
  }

  feats <- .compute_hydro_features(
    hydro, target_date, m$tau_short, m$tau_long, hydro_type
  )
  newdata <- tibble::tibble(
    doy         = as.integer(format(target_date, "%j")),
    hydro_short = feats$hydro_short,
    hydro_long  = feats$hydro_long
  )
  log_pred <- tryCatch(
    as.numeric(stats::predict(m$gamm_fit, newdata = newdata)),
    error = function(e) NA_real_
  )
  if (is.na(log_pred)) {
    return(list(ref_norm = m$static_ref %||% NA_real_, ref_tier = "static"))
  }
  list(ref_norm = exp(log_pred), ref_tier = "model")
}

#' Resolve ref_norm for instantaneous target samples (one date per sample)
#'
#' @param ref_model A `reference_model` object.
#' @param target_df Tibble with columns `sample_id`, `datetime` (Date).
#' @return Tibble `(sample_id, analyte, ref_norm, ref_tier)`.
#' @keywords internal
.resolve_ref_norm_instant <- function(ref_model, target_df) {
  target_df <- dplyr::mutate(target_df, .date = as.Date(.data$datetime))
  sample_dates <- dplyr::distinct(target_df, .data$sample_id, .data$.date)

  results <- purrr::pmap_dfr(sample_dates, function(sample_id, .date) {
    target_date <- .date
    target_event_api <- .compute_api(
      ref_model$hydro$value, as.Date(ref_model$hydro$date),
      target_date, ref_model$match_window_days
    )

    purrr::map_dfr(names(ref_model$models), function(nm) {
      m <- ref_model$models[[nm]]

      # Tier 1: direct contemporaneous match
      obs <- m$obs
      if (!is.null(obs) && nrow(obs) > 0L) {
        dt_diff <- abs(as.numeric(obs$date - target_date))
        time_ok <- dt_diff <= ref_model$match_window_days

        # Use short-window hydro for the gate; stored as hydro_short on the obs
        hydro_diff <- abs(obs$hydro_short - target_event_api)
        hydro_ok <- hydro_diff <= ref_model$match_hydro_tol

        candidate_mask <- time_ok & hydro_ok
        if (any(candidate_mask, na.rm = TRUE)) {
          best_idx <- which(candidate_mask)[which.min(dt_diff[candidate_mask])]
          return(tibble::tibble(
            sample_id = sample_id,
            analyte   = nm,
            ref_norm  = obs$value_norm[best_idx],
            ref_tier  = "matched"
          ))
        }
      }

      # Tier 2 / 3
      pred <- .predict_ref_at_date(
        m, target_date, ref_model$hydro, ref_model$hydro_type
      )
      tibble::tibble(
        sample_id = sample_id,
        analyte   = nm,
        ref_norm  = pred$ref_norm,
        ref_tier  = pred$ref_tier
      )
    })
  })

  results
}

#' Resolve ref_norm for chronic target samples (window-integrated)
#'
#' For each chronic sample (`focal_date` column), integrates the predicted
#' reference at daily resolution over `[focal_date - window_days, focal_date]`
#' using the same exponential-decay kernel as `time_weighted_aggregate()`.
#' Tier 1 does not apply to chronic targets (no single reference grab is a
#' proxy for an entire integration window).
#'
#' @param ref_model A `reference_model` object.
#' @param target_df Tibble with columns `sample_id`, `focal_date`.
#' @param tau_days Exponential-decay parameter (days).
#' @param window_days Look-back window (days).
#' @return Tibble `(sample_id, analyte, ref_norm, ref_tier)`.
#' @keywords internal
.resolve_ref_norm_chronic <- function(ref_model, target_df, tau_days, window_days) {
  focal_df <- dplyr::distinct(target_df, .data$sample_id, .data$focal_date) |>
    dplyr::mutate(focal_date = as.Date(.data$focal_date))

  purrr::pmap_dfr(focal_df, function(sample_id, focal_date) {
    dates_in_window <- seq(focal_date - window_days, focal_date, by = "day")
    # Kernel weights: pure exponential decay (delta_t = 1 for daily series)
    w <- exp(-as.numeric(focal_date - dates_in_window) / tau_days)
    w <- w / sum(w)

    purrr::map_dfr(names(ref_model$models), function(nm) {
      m <- ref_model$models[[nm]]

      # Tier 3 fallback: integrate static scalar (result = same scalar)
      if (m$tier == "static" || is.null(m$gamm_fit)) {
        return(tibble::tibble(
          sample_id = sample_id,
          analyte   = nm,
          ref_norm  = m$static_ref %||% NA_real_,
          ref_tier  = "static"
        ))
      }

      # Tier 2: daily GAM predictions over the window
      feats <- .compute_hydro_features(
        ref_model$hydro, dates_in_window,
        m$tau_short, m$tau_long, ref_model$hydro_type
      )
      newdata <- tibble::tibble(
        doy         = as.integer(format(dates_in_window, "%j")),
        hydro_short = feats$hydro_short,
        hydro_long  = feats$hydro_long
      )
      log_preds <- tryCatch(
        as.numeric(stats::predict(m$gamm_fit, newdata = newdata)),
        error = function(e) rep(NA_real_, length(dates_in_window))
      )

      if (all(is.na(log_preds))) {
        return(tibble::tibble(
          sample_id = sample_id,
          analyte   = nm,
          ref_norm  = m$static_ref %||% NA_real_,
          ref_tier  = "static"
        ))
      }

      # Replace NAs with static ref on log scale; weight accordingly
      eps <- 1e-9
      log_static <- log(max(m$static_ref %||% eps, eps))
      log_preds[is.na(log_preds)] <- log_static

      # Kernel-weighted geometric mean (matches time_weighted_aggregate geom_mean)
      ref_norm_val <- exp(sum(w * log_preds))

      tibble::tibble(
        sample_id = sample_id,
        analyte   = nm,
        ref_norm  = ref_norm_val,
        ref_tier  = "model_integrated"
      )
    })
  })
}

#' Resolve reference norms for a block of target samples
#'
#' Dispatches to [.resolve_ref_norm_instant()] or [.resolve_ref_norm_chronic()]
#' depending on whether `df` contains a `focal_date` column.  Returns a tibble
#' `(sample_id, analyte, ref_norm, ref_tier)` that [add_amspaf()] uses instead
#' of the static `(analyte, ref_norm)` table produced by [prepare_reference()].
#'
#' @param ref_model A `reference_model` object.
#' @param df Target chemistry data frame (as passed to [add_amspaf()]).
#' @param tau_days Exponential-decay parameter for chronic integration.
#' @param window_days Window length for chronic integration.
#' @return Tibble `(sample_id, analyte, ref_norm, ref_tier)`.
#' @keywords internal
.resolve_ref_norm <- function(ref_model, df, tau_days = 90, window_days = 365) {
  is_chronic <- "focal_date" %in% names(df)
  if (is_chronic) {
    .resolve_ref_norm_chronic(ref_model, df, tau_days, window_days)
  } else {
    .resolve_ref_norm_instant(ref_model, df)
  }
}


## ============================================================================
## ara_summary() accessor
## ============================================================================

#' Retrieve per-cell ARA diagnostics from an `add_amspaf()` result
#'
#' After calling [add_amspaf()] with a `reference_model` (or any reference),
#' this accessor returns a tibble describing what happened in the ARA
#' subtraction for every (sample × analyte) that was assessed.  This is the
#' primary tool for auditing the "reference higher than target" case (floored
#' to zero) and for understanding which tier was used per cell.
#'
#' The attribute is stored by [add_amspaf()] and is dropped by most dplyr
#' verbs, so read the summary before further wrangling.
#'
#' @param x A data frame returned by [add_amspaf()].
#' @return A tibble with columns:
#'   \describe{
#'     \item{`sample_id`}{Sample identifier.}
#'     \item{`analyte`}{Analyte name.}
#'     \item{`ref_norm`}{Normalised reference concentration subtracted.}
#'     \item{`C_norm`}{Normalised target concentration (before ARA).}
#'     \item{`C_adj`}{ARA-adjusted concentration (`max(C_norm - ref_norm, 0)`).}
#'     \item{`C_excess`}{Unfloored difference `C_norm - ref_norm`; negative
#'       values indicate the reference exceeded the target — possibly a
#'       geogenic artefact (e.g. low-pH upstream metal mobilisation).}
#'     \item{`floor_fired`}{Logical; `TRUE` when `C_norm < ref_norm`.}
#'     \item{`ref_source`}{`"disabled"`, `"matched"`, or `"unmatched"`.}
#'     \item{`ref_tier`}{`"matched"`, `"model"`, `"model_integrated"`, or
#'       `"static"`.  `NA` for non-temporal reference.}
#'   }
#'   Returns `NULL` (with a message) if the attribute is absent.
#'
#' @seealso [add_amspaf()], [fit_reference_model()]
#' @export
ara_summary <- function(x) {
  s <- attr(x, "ara_summary")
  if (is.null(s)) {
    cli::cli_inform("No {.field ara_summary} attribute found. \\
                    Call {.fn ara_summary} on the direct output of \\
                    {.fn add_amspaf} before any further wrangling.")
    return(NULL)
  }
  s
}

#' Per-analyte PAF breakdown from add_amspaf()
#'
#' After calling [add_amspaf()], this accessor returns the per-analyte
#' contribution breakdown behind each AmsPAF value: the ARA-adjusted
#' concentration, single-substance PAF, MOA group and reference source for every
#' assessed (sample × draw × analyte). It replaces the former `analyte_pafs`
#' list-column with a flat, tidy frame (filter/join directly).
#'
#' The attribute is stored by [add_amspaf()] and is dropped by most dplyr verbs,
#' so read it before further wrangling.
#'
#' @param x A data frame returned by [add_amspaf()].
#' @return A tibble with columns `site_id`, `sample_id`, `draw_id` (draws mode
#'   only), `analyte`, `C_adj`, `PAF`, `moa_group`, `ref_source`. Returns `NULL`
#'   (with a message) if the attribute is absent.
#' @seealso [add_amspaf()], [ara_summary()]
#' @export
analyte_pafs <- function(x) {
  s <- attr(x, "analyte_pafs")
  if (is.null(s)) {
    cli::cli_inform("No {.field analyte_pafs} attribute found. \\
                    Call {.fn analyte_pafs} on the direct output of \\
                    {.fn add_amspaf} before any further wrangling.")
    return(NULL)
  }
  s
}
