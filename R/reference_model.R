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
##   .compute_api()                  API (Antecedent Precipitation Index)
##   .compute_antecedent_mean()      rolling mean for stage/discharge
##   .compute_hydro_features()       dispatcher; returns (date, hydro_short, hydro_long)
##   .build_analyte_model_df()       assemble (y, doy, hydro_short, hydro_long) for one analyte
##   .fit_ref_gam()                  fit mgcv::gam with gamm fallback
##   .aic_for_windows()              AIC for a (window_short, window_long) candidate pair
##   .select_api_windows()           AIC-based window selection per analyte
##   .resolve_ref_norm()             three-tier resolver dispatcher
##   .resolve_ref_norm_instant()     per-sample (instantaneous) resolver
##   .resolve_ref_norm_chronic()     per-sample (chronic window-integrated) resolver
##   .predict_ref_at_date()          predict value_norm from a fitted analyte model


## ============================================================================
## Internal hydrology helpers
## ============================================================================

#' Antecedent Precipitation Index (API) for rainfall
#'
#' Computes `API(t) = sum_{i=0}^{window} P(t-i) * k^i` where the decay
#' constant `k = exp(-1 / window_days)` — meaning rainfall `window_days` days
#' ago carries weight `1/e ≈ 0.37`.  Returns 0 for target dates with no hydro
#' record in the preceding window.
#'
#' @param hydro_values Numeric vector of daily rainfall values.
#' @param hydro_dates Date vector matching `hydro_values`.
#' @param target_dates Date vector; API is evaluated at each.
#' @param window_days Integer memory length (days).
#' @return Numeric vector, same length as `target_dates`.
#' @keywords internal
.compute_api <- function(hydro_values, hydro_dates, target_dates, window_days) {
  k <- exp(-1 / window_days)
  vapply(target_dates, function(td) {
    start <- td - window_days
    mask  <- hydro_dates >= start & hydro_dates <= td
    if (!any(mask)) return(0)
    lags <- as.numeric(td - hydro_dates[mask])
    sum(hydro_values[mask] * k ^ lags, na.rm = TRUE)
  }, numeric(1L))
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
    if (!any(mask)) return(NA_real_)
    mean(hydro_values[mask], na.rm = TRUE)
  }, numeric(1L))
}

#' Compute hydrology features (hydro_short, hydro_long) at target dates
#'
#' Dispatches to [.compute_api()] for rainfall or [.compute_antecedent_mean()]
#' for stage/discharge, at both the short and long window lengths.
#'
#' @param hydro Data frame with columns `date` (Date) and `value` (numeric).
#' @param target_dates Date vector.
#' @param window_short,window_long Memory lengths (days).
#' @param hydro_type `"rainfall"`, `"stage"`, or `"discharge"`.
#' @return Tibble `(date, hydro_short, hydro_long)`.
#' @keywords internal
.compute_hydro_features <- function(hydro, target_dates,
                                     window_short, window_long,
                                     hydro_type = "rainfall") {
  hdates <- as.Date(hydro$date)
  hvals  <- hydro$value

  if (hydro_type == "rainfall") {
    short <- .compute_api(hvals, hdates, target_dates, window_short)
    long  <- .compute_api(hvals, hdates, target_dates, window_long)
  } else {
    short <- .compute_antecedent_mean(hvals, hdates, target_dates, window_short)
    long  <- .compute_antecedent_mean(hvals, hdates, target_dates, window_long)
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
#' @param window_short,window_long Memory window lengths (days).
#' @param hydro_type Hydro type string.
#' @param eps Guard for log transform.
#' @return Tibble `(y, doy, hydro_short, hydro_long)` sorted by date.
#' @keywords internal
.build_analyte_model_df <- function(obs, hydro, window_short, window_long,
                                     hydro_type, eps = 1e-9) {
  obs <- obs[order(obs$date), ]
  feats <- .compute_hydro_features(hydro, obs$date, window_short, window_long, hydro_type)

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
.fit_ref_gam <- function(df_model, eps = 1e-9) {  # nolint: unused arg
  n <- nrow(df_model)
  if (n < 10L) return(NULL)

  k_doy <- min(6L, n - 2L)
  k_h   <- min(4L, n - 2L)

  tryCatch(
    mgcv::gam(
      y ~ s(doy, bs = "cc", k = k_doy) +
          s(hydro_short, k = k_h) +
          s(hydro_long,  k = k_h),
      data   = df_model,
      knots  = list(doy = c(1L, 365L)),
      method = "REML"
    ),
    error   = function(e) NULL,
    warning = function(w) {
      tryCatch(
        mgcv::gam(
          y ~ s(doy, bs = "cc", k = k_doy) +
              s(hydro_short, k = k_h) +
              s(hydro_long,  k = k_h),
          data   = df_model,
          knots  = list(doy = c(1L, 365L)),
          method = "REML"
        ),
        error = function(e) NULL
      )
    }
  )
}

#' AIC for a single (window_short, window_long) candidate
#'
#' Returns `Inf` if fitting fails.
#' @keywords internal
.aic_for_windows <- function(obs, hydro, window_short, window_long, hydro_type, eps) {
  df_m <- tryCatch(
    .build_analyte_model_df(obs, hydro, window_short, window_long, hydro_type, eps),
    error = function(e) NULL
  )
  if (is.null(df_m)) return(Inf)
  fit <- .fit_ref_gam(df_m, eps)
  if (is.null(fit)) return(Inf)
  stats::AIC(fit)
}

#' Select API memory windows per analyte using AIC
#'
#' Evaluates every `(window_short, window_long)` candidate pair and picks the
#' pair with the lowest AIC.  Among pairs within 2 ΔAIC units of the best,
#' the most parsimonious (smallest `window_short + window_long`) is chosen.
#'
#' Returns a list `(window_short, window_long, best_aic, null_aic)`.  If all
#' candidate fits fail, returns the first candidate pair with `best_aic = Inf`.
#'
#' @keywords internal
.select_api_windows <- function(obs, hydro, hydro_type,
                                 api_windows_short, api_windows_long,
                                 eps = 1e-9) {
  grid <- expand.grid(
    window_short = api_windows_short,
    window_long  = api_windows_long,
    stringsAsFactors = FALSE
  )
  grid <- grid[grid$window_short < grid$window_long, , drop = FALSE]
  if (nrow(grid) == 0L) {
    grid <- data.frame(
      window_short = api_windows_short[1L],
      window_long  = api_windows_long[length(api_windows_long)]
    )
  }

  aics <- mapply(
    function(ws, wl) .aic_for_windows(obs, hydro, ws, wl, hydro_type, eps),
    grid$window_short, grid$window_long,
    SIMPLIFY = TRUE
  )

  # Null model AIC (intercept only)
  df0 <- tryCatch(
    .build_analyte_model_df(
      obs, hydro, api_windows_short[1L], api_windows_long[1L], hydro_type, eps
    ),
    error = function(e) NULL
  )
  null_aic <- if (!is.null(df0)) {
    fit0 <- tryCatch(mgcv::gam(y ~ 1, data = df0, method = "REML"),
                     error = function(e) NULL)
    if (!is.null(fit0)) stats::AIC(fit0) else Inf
  } else Inf

  best_aic <- min(aics, na.rm = TRUE)
  if (is.infinite(best_aic)) {
    return(list(
      window_short = grid$window_short[1L],
      window_long  = grid$window_long[1L],
      best_aic     = Inf,
      null_aic     = null_aic
    ))
  }

  # Among rows within 2 ΔAIC of best, pick the most parsimonious
  within2 <- which(aics <= best_aic + 2)
  candidate_rows <- grid[within2, , drop = FALSE]
  candidate_rows$.total <- candidate_rows$window_short + candidate_rows$window_long
  chosen <- candidate_rows[which.min(candidate_rows$.total), , drop = FALSE]

  list(
    window_short = chosen$window_short,
    window_long  = chosen$window_long,
    best_aic     = best_aic,
    null_aic     = null_aic
  )
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
#'     `s(doy, bs="cc") + s(hydro_short) + s(hydro_long)`.  Window lengths are
#'     auto-selected by AIC over the `api_windows_short × api_windows_long`
#'     candidate grid.  The analyte falls back to tier 3 if the best model has
#'     higher AIC than the null (intercept-only) model, or if fewer than
#'     `min_obs_model` detected observations are available.
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
#' @param match_window_days Integer; tier-1 time tolerance in days (default 5).
#' @param match_hydro_tol Numeric; tier-1 hydro tolerance (default `NULL` →
#'   `0.5 × IQR` of the reference event-API series).
#' @param api_windows_short Integer vector of candidate short-memory window
#'   lengths in days (default `c(3L, 7L, 14L)`).
#' @param api_windows_long Integer vector of candidate long-memory window
#'   lengths in days (default `c(30L, 60L, 90L, 180L)`).
#' @param auto_select Logical; if `TRUE` (default), select window lengths per
#'   analyte by AIC.  If `FALSE`, use `api_windows_short[1]` and
#'   `api_windows_long[1]` for all analytes without CV.
#' @param min_obs_model Integer; minimum detected observations required to
#'   attempt tier-2 modelling (default `20L`).  Analytes below this fall
#'   directly to tier 3.
#' @param summary Static-fallback summary statistic (tier 3): one of
#'   `"geom_mean"` (default), `"median"`, `"arith_mean"`, `"p80"`, `"p90"`,
#'   `"p95"`.
#' @param silo_start,silo_end Start/end dates for SILO auto-fetch.  Default
#'   `NULL`: derived from the reference chemistry date range (padded by 365
#'   days on the left for the maximum API window).
#' @param silo_api_key Passed to [get_silo_rainfall()] when `hydro = NULL`.
#' @param eps Small positive guard for log transform (default `1e-9`).
#'
#' @return An object of class `reference_model`:
#'   \describe{
#'     \item{`$models`}{Named list (one per analyte) carrying `gamm_fit`,
#'       `window_short`, `window_long`, `best_aic`, `null_aic`, `tier`,
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
#'   hydro      = my_stage_df,  # data.frame(date, value)
#'   hydro_type = "stage",
#'   conc_units = "ug/L"
#' )
#' }
#' @export
fit_reference_model <- function(
    reference,
    hydro              = NULL,
    hydro_type         = "rainfall",
    latitude           = NULL,
    longitude          = NULL,
    conc_units         = NULL,
    analyte_metadata   = NULL,
    match_window_days  = 5L,
    match_hydro_tol    = NULL,
    api_windows_short  = c(3L, 7L, 14L),
    api_windows_long   = c(30L, 60L, 90L, 180L),
    auto_select        = TRUE,
    min_obs_model      = 20L,
    summary            = "geom_mean",
    silo_start         = NULL,
    silo_end           = NULL,
    silo_api_key       = NULL,
    eps                = 1e-9
) {
  ## ── Input validation ───────────────────────────────────────────────────────
  checkmate::assert_data_frame(reference)
  checkmate::assert_names(names(reference),
    must.include = c("sample_id", "datetime", "analyte", "value", "detected"))
  checkmate::assert_int(match_window_days, lower = 0L)
  checkmate::assert_flag(auto_select)
  checkmate::assert_int(min_obs_model, lower = 5L)
  summary <- match.arg(summary,
    c("geom_mean", "median", "arith_mean", "p80", "p90", "p95"))
  checkmate::assert_number(eps, lower = 0)

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
    s_start <- if (!is.null(silo_start)) as.Date(silo_start) else
      min(ref_dates, na.rm = TRUE) - max(api_windows_long)
    s_end <- if (!is.null(silo_end)) as.Date(silo_end) else
      max(ref_dates, na.rm = TRUE)

    cli::cli_inform("Fetching SILO rainfall for ({latitude}, {longitude}) \\
                    {s_start}\u2013{s_end}\u2026")
    rain_df <- get_silo_rainfall(
      latitude   = latitude,
      longitude  = longitude,
      start_date = s_start,
      end_date   = s_end,
      api_key    = silo_api_key
    )
    hydro      <- dplyr::rename(rain_df, value = "rainfall_mm")
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

  hydro <- dplyr::mutate(hydro, date = as.Date(.data$date),
                         value = as.numeric(.data$value))
  hydro <- hydro[order(hydro$date), ]

  ## ── Unit conversion + normalisation ───────────────────────────────────────
  meta <- .load_analyte_metadata(analyte_metadata)
  ssd_analytes <- meta$analyte[
    !is.na(meta$ssd_available) & meta$ssd_available == TRUE &
    !meta$analyte %in% .AMSPAF_EXCLUDED_ANALYTES
  ]
  reference <- .convert_df_tox_to_ugL(reference, ssd_analytes, conc_units, "reference")

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
  norm_det <- dplyr::filter(norm_all, .data$detected, !is.na(.data$value_norm),
                            .data$value_norm > 0, .data$analyte %in% ssd_analytes)

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

    # Add hydro features at observation dates (longest candidate windows)
    # for tier-1 matching; stored regardless of tier
    max_short <- max(api_windows_short)
    max_long  <- max(api_windows_long)
    obs_feats <- .compute_hydro_features(
      hydro, obs$date, max_short, max_long, hydro_type
    )
    obs_with_feats <- dplyr::bind_cols(obs, dplyr::select(obs_feats, -"date"))

    if (n_obs < min_obs_model) {
      models[[nm]] <- list(
        gamm_fit     = NULL,
        window_short = api_windows_short[1L],
        window_long  = api_windows_long[length(api_windows_long)],
        best_aic     = NA_real_,
        null_aic     = NA_real_,
        tier         = "static",
        n_obs        = n_obs,
        static_ref   = static_val,
        obs          = obs_with_feats
      )
      next
    }

    if (auto_select) {
      sel <- .select_api_windows(
        obs, hydro, hydro_type, api_windows_short, api_windows_long, eps
      )
    } else {
      df0 <- .build_analyte_model_df(
        obs, hydro, api_windows_short[1L], api_windows_long[1L], hydro_type, eps
      )
      fit0 <- tryCatch(
        mgcv::gam(y ~ 1, data = df0, method = "REML"),
        error = function(e) NULL
      )
      sel <- list(
        window_short = api_windows_short[1L],
        window_long  = api_windows_long[1L],
        best_aic     = if (!is.null(fit0)) stats::AIC(fit0) else Inf,
        null_aic     = if (!is.null(fit0)) stats::AIC(fit0) else Inf
      )
    }

    df_m <- .build_analyte_model_df(
      obs, hydro, sel$window_short, sel$window_long, hydro_type, eps
    )
    gam_fit <- .fit_ref_gam(df_m, eps)

    # Model is useful only if it beats the null
    tier <- if (
      is.null(gam_fit) ||
      is.infinite(sel$best_aic) ||
      sel$best_aic >= sel$null_aic
    ) "static" else "model"

    models[[nm]] <- list(
      gamm_fit     = gam_fit,
      window_short = sel$window_short,
      window_long  = sel$window_long,
      best_aic     = sel$best_aic,
      null_aic     = sel$null_aic,
      tier         = tier,
      n_obs        = n_obs,
      static_ref   = static_val,
      obs          = obs_with_feats
    )
  }

  ## ── Tier-1 hydro tolerance ────────────────────────────────────────────────
  # Compute using the shortest short-window API on the reference observations
  if (is.null(match_hydro_tol)) {
    all_obs_dates <- as.Date(reference$datetime)
    event_api <- .compute_api(
      hydro$value, as.Date(hydro$date),
      all_obs_dates, api_windows_short[1L]
    )
    match_hydro_tol <- 0.5 * stats::IQR(event_api, na.rm = TRUE)
    if (match_hydro_tol < 1e-6) match_hydro_tol <- Inf  # no variability → no gate
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
  n_model    <- sum(vapply(x$models, function(m) m$tier == "model", logical(1L)))
  n_static   <- n_analytes - n_model

  hydro_range <- if (!is.null(x$hydro) && nrow(x$hydro) > 0L) {
    sprintf("%s \u2013 %s", min(x$hydro$date), max(x$hydro$date))
  } else "unknown"

  cat(sprintf(
    "<reference_model>  fitted %s | %d analyte%s | hydro: %s (%s)\n",
    x$fit_date, n_analytes,
    if (n_analytes == 1L) "" else "s",
    x$hydro_type, hydro_range
  ))
  cat(sprintf("  tier-1:  \u00b1%d d + API \u00b1%.2g\n",
              x$match_window_days, x$match_hydro_tol))

  if (n_model > 0L) {
    model_nms <- names(Filter(function(m) m$tier == "model", x$models))
    detail <- vapply(model_nms, function(nm) {
      m <- x$models[[nm]]
      sprintf("%s (w=%d/%dd, n=%d)", nm, m$window_short, m$window_long, m$n_obs)
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
    hydro, target_date, m$window_short, m$window_long, hydro_type
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
    target_date    <- .date
    target_event_api <- .compute_api(
      ref_model$hydro$value, as.Date(ref_model$hydro$date),
      target_date, ref_model$match_window_days
    )

    purrr::map_dfr(names(ref_model$models), function(nm) {
      m <- ref_model$models[[nm]]

      # Tier 1: direct contemporaneous match
      obs <- m$obs
      if (!is.null(obs) && nrow(obs) > 0L) {
        dt_diff  <- abs(as.numeric(obs$date - target_date))
        time_ok  <- dt_diff <= ref_model$match_window_days

        # Use short-window hydro for the gate; stored as hydro_short on the obs
        hydro_diff <- abs(obs$hydro_short - target_event_api)
        hydro_ok   <- hydro_diff <= ref_model$match_hydro_tol

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
        m$window_short, m$window_long, ref_model$hydro_type
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
