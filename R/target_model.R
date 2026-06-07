## target_model.R
##
## Season-blind, hydrology-modulated predictive model of the *site impact* at a
## target (impacted) site.  Issue #14.
##
## Where fit_reference_model() (#9) predicts the natural background, this model
## predicts the anthropogenic increment
##
##     I(t) = C_norm(t) - ref_norm(t)            (== ara_summary()'s C_excess)
##
## so that days without a full toxicant suite get a chemistry-grounded estimate
## rather than a forward-filled one.
##
## The critical asymmetry vs. the reference model: day-of-year may NEVER generate
## impact.  Site impacts are driven by management failure (a breach), not the
## calendar; rain enables an existing breach to escape but does not cause one.
## So the impact model is season-blind — hydrology enters only as a *modulator*
## of how an already-present impact expresses itself (first-flush mobilisation,
## sorption, antecedent memory), never as a generator.
##
##     I_metal(t) = beta_metal . f(hydro_t)      (non-conservative expression)
##                + S_metal(t)                   (latent state, bridge-interpolated)
##
## This is a season-blind, hydrology-modulated variant of WRTDS-Kalman
## (Zhang & Hirsch 2019, WRR 55(11):9705-9723) applied to the ARA impact
## residual; the hydrological response f(hydro) follows concentration-discharge
## theory (Godsey, Kirchner & Clow 2009, Hydrol. Process. 23:1844-1864).
##
## Public surface
## ──────────────
##   fit_target_model()      fit per-analyte impact-residual models
##   print.target_model()    S3 print method
##
## Internal
## ────────
##   .fit_impact_response()  season-blind GAM I ~ s(hydro_short) + s(hydro_long)
##   .resolve_target_impact() predict impact (+ implied C_norm) at query dates
##   .interp_residual()      hydrology-weighted bracketing-anchor bridge


## ============================================================================
## fit_target_model()
## ============================================================================

#' Fit a season-blind predictive model of the site impact
#'
#' Models the anthropogenic increment `I = C_norm - ref_norm` (the ARA
#' "added risk", i.e. `ara_summary()`'s `C_excess`) at a target site as a
#' function of hydrology and a persistent latent state — **never** of
#' day-of-year.  Used by [amspaf_daily()] (`interpolation = "model"`) to fill
#' the gaps between grab samples with a chemistry-grounded impact estimate
#' instead of a forward-filled concentration.
#'
#' **Why season-blind.** Site impacts are driven by management failure (a
#' leachate breach), not the calendar.  Heavy summer rain *enables* an existing
#' breach to escape; it does not *cause* one.  A perfect-management year has
#' zero impact under identical rainfall.  Day-of-year is therefore a confounded
#' covariate for impact and is excluded; hydrology enters only as a modulator of
#' how an already-present impact expresses itself (first-flush mobilisation,
#' sorption, antecedent memory).
#'
#' **Model (per analyte).**
#' \deqn{I(t) = \beta \cdot f(\mathrm{hydro}_t) + S(t)}
#' where `f(hydro)` is a thin-plate GAM on the short- and long-window
#' antecedent indices (no cyclic seasonal term), fitted only when it beats an
#' intercept-only null by AIC and at least `min_obs_model` anchors are
#' available; otherwise `beta = 0` and the analyte falls to a pure
#' state-interpolation ("bridge") tier.  The residual state `S(t)` is
#' interpolated between observation anchors by a hydrology-weighted bridge that
#' pinches to the observed residual at each anchor and leans toward the
#' hydrologically more similar bracketing anchor in between.
#'
#' This is a season-blind, hydrology-modulated variant of WRTDS-Kalman
#' (Zhang & Hirsch 2019); `f(hydro)` follows concentration–discharge theory
#' (Godsey, Kirchner & Clow 2009).
#'
#' @param target Long-format target chemistry. Required columns: `sample_id`,
#'   `datetime`, `analyte`, `value`, `detected`. Toxicants must be in µg/L;
#'   supply via a `units.analyte` column or `conc_units`. Co-analyte rows
#'   (pH, DOC, hardness, temperature) should be present for normalisation.
#' @param reference_model A `reference_model` from [fit_reference_model()].
#'   Supplies `ref_norm(t)` and (by default) the catchment hydrology series.
#' @param hydro Optional target-specific daily hydrology data frame
#'   (`date`, `value`); when `NULL` (default) the reference model's hydrology is
#'   reused (shared-catchment assumption). A target-local stage/discharge gauge
#'   can capture breach-mobilising flow that catchment rainfall misses.
#' @param hydro_type `"rainfall"`, `"stage"`, or `"discharge"`; used only when
#'   `hydro` is supplied. Default `"rainfall"`.
#' @param imputation_model Optional `imputation_model` from
#'   [fit_imputation_model()] (fit on the target's chemistry). When supplied,
#'   missing analytes are imputed in raw concentration space before the impact
#'   is computed, adding more anchor days (tier 2). Requires **brms**.
#' @param conc_units Unit string for target chemistry when no `units.analyte`
#'   column is present.
#' @param analyte_metadata Analyte metadata, or `NULL` for the bundled CSV.
#' @param api_windows_short,api_windows_long Candidate short/long antecedent
#'   memory windows (days) for `f(hydro)`, selected by AIC.
#' @param auto_select Logical; AIC window selection per analyte (default `TRUE`).
#' @param min_obs_model Integer; minimum impact anchors required to attempt the
#'   `f(hydro)` GAM. Below this, the analyte uses the bridge tier. Default `12L`.
#' @param eps Small positive guard. Default `1e-9`.
#'
#' @return An object of class `target_model`:
#'   \describe{
#'     \item{`$models`}{Named per-analyte list: `impact_fit` (gam or `NULL`),
#'       `window_short`, `window_long`, `tier` (`"model"` or `"bridge"`),
#'       `n_obs`, and `anchors` (tibble `date`, `I`, `S`, `hydro_short`,
#'       `hydro_long`).}
#'     \item{`$reference_model`}{The supplied reference model.}
#'     \item{`$hydro`,`$hydro_type`}{Hydrology series used for `f(hydro)`.}
#'     \item{`$fit_date`}{Date fitted.}
#'   }
#'
#' @seealso [amspaf_daily()], [fit_reference_model()], [add_amspaf()],
#'   [ara_summary()]
#'
#' @references
#' Zhang Q, Hirsch RM (2019) Water Resources Research 55(11):9705–9723.
#' Godsey SE, Kirchner JW, Clow DW (2009) Hydrological Processes 23:1844–1864.
#'
#' @examples
#' \dontrun{
#' ref_model <- fit_reference_model(reference_chem, latitude = -33.8,
#'                                  longitude = 151.2, conc_units = "ug/L")
#' tgt_model <- fit_target_model(target_chem, ref_model, conc_units = "ug/L")
#' tgt_model
#' }
#' @export
fit_target_model <- function(
    target,
    reference_model,
    hydro              = NULL,
    hydro_type         = "rainfall",
    imputation_model   = NULL,
    conc_units         = NULL,
    analyte_metadata   = NULL,
    api_windows_short  = c(3L, 7L, 14L),
    api_windows_long   = c(30L, 60L, 90L, 180L),
    auto_select        = TRUE,
    min_obs_model      = 12L,
    eps                = 1e-9
) {
  ## ── Validation ─────────────────────────────────────────────────────────────
  checkmate::assert_data_frame(target)
  checkmate::assert_names(names(target),
    must.include = c("sample_id", "datetime", "analyte", "value", "detected"))
  if (!inherits(reference_model, "reference_model")) {
    cli::cli_abort(
      "{.arg reference_model} must be a {.cls reference_model} from \\
       {.fn fit_reference_model}."
    )
  }
  if (!is.null(imputation_model) &&
      !inherits(imputation_model, "imputation_model")) {
    cli::cli_abort(
      "{.arg imputation_model} must be an {.cls imputation_model} or {.val NULL}."
    )
  }
  checkmate::assert_int(min_obs_model, lower = 4L)

  ## ── Hydrology (default: reuse the reference model's series) ─────────────────
  if (is.null(hydro)) {
    hydro      <- reference_model$hydro
    hydro_type <- reference_model$hydro_type
  } else {
    checkmate::assert_names(names(hydro), must.include = c("date", "value"))
    if ("type" %in% names(hydro)) hydro_type <- unique(hydro$type)[1L]
    hydro_type <- match.arg(hydro_type, c("rainfall", "stage", "discharge"))
    hydro <- dplyr::select(hydro, date = "date", value = "value")
  }
  hydro <- dplyr::mutate(hydro, date = as.Date(.data$date),
                         value = as.numeric(.data$value))
  hydro <- hydro[order(hydro$date), ]

  ## ── Units, optional impute-first, normalisation ────────────────────────────
  meta <- .load_analyte_metadata(analyte_metadata)
  ssd_analytes <- meta$analyte[
    !is.na(meta$ssd_available) & meta$ssd_available == TRUE &
    !meta$analyte %in% .AMSPAF_EXCLUDED_ANALYTES
  ]
  target <- .convert_df_tox_to_ugL(target, ssd_analytes, conc_units, "target")

  if (!is.null(imputation_model)) {
    if (!"site_id" %in% names(target)) target$site_id <- "target"
    target <- impute_chemistry(target, imputation_model)
  }

  target <- dplyr::mutate(target, .date = as.Date(.data$datetime))
  norm_all <- .normalise_ref_observations(target, target, meta)
  norm_det <- dplyr::filter(
    norm_all, .data$detected, !is.na(.data$value_norm),
    .data$value_norm > 0, .data$analyte %in% ssd_analytes
  )

  ## ── Reference norm at each target sample (instant resolver) ────────────────
  ref_resolved <- .resolve_ref_norm_instant(
    reference_model,
    dplyr::distinct(target, .data$sample_id, .data$datetime)
  )

  ## ── Impact anchors: I = value_norm - ref_norm ──────────────────────────────
  anchors_all <- norm_det |>
    dplyr::select("sample_id", date = ".date", "analyte", "value_norm") |>
    dplyr::inner_join(
      dplyr::select(ref_resolved, "sample_id", "analyte", "ref_norm"),
      by = c("sample_id", "analyte")
    ) |>
    dplyr::mutate(I = .data$value_norm - .data$ref_norm) |>
    dplyr::filter(is.finite(.data$I)) |>
    ## one anchor per (date, analyte) — average duplicate same-day grabs
    dplyr::summarise(
      I = mean(.data$I, na.rm = TRUE),
      .by = c("date", "analyte")
    )

  ## ── Per-analyte fit ────────────────────────────────────────────────────────
  target_analytes <- intersect(unique(anchors_all$analyte), ssd_analytes)
  models <- vector("list", length(target_analytes))
  names(models) <- target_analytes

  for (nm in target_analytes) {
    obs <- anchors_all |>
      dplyr::filter(.data$analyte == .env$nm) |>
      dplyr::arrange(.data$date)

    feats <- .compute_hydro_features(
      hydro, obs$date, max(api_windows_short), max(api_windows_long), hydro_type
    )
    obs <- dplyr::bind_cols(obs, dplyr::select(feats, -"date"))

    fit_res <- .fit_impact_response(
      obs, hydro, hydro_type, api_windows_short, api_windows_long,
      auto_select, min_obs_model, eps
    )
    models[[nm]] <- fit_res
  }

  structure(
    list(
      models       = models,
      reference_model = reference_model,
      hydro        = hydro,
      hydro_type   = hydro_type,
      fit_date     = Sys.Date()
    ),
    class = "target_model"
  )
}


## ============================================================================
## Internal: season-blind hydrological-response fit
## ============================================================================

#' Fit the season-blind impact response for one analyte
#'
#' Selects short/long antecedent windows by AIC, fits `I ~ s(hydro_short) +
#' s(hydro_long)` (no seasonal term), and keeps it only if it beats the
#' intercept-only null.  Stores the de-trended residual `S = I - fitted` as the
#' bridge-interpolation state.
#'
#' @return List with `impact_fit`, `window_short`, `window_long`, `tier`,
#'   `n_obs`, `anchors` (`date`, `I`, `S`, `hydro_short`, `hydro_long`).
#' @keywords internal
.fit_impact_response <- function(obs, hydro, hydro_type,
                                  api_windows_short, api_windows_long,
                                  auto_select, min_obs_model, eps) {
  n_obs <- nrow(obs)

  flat <- function(ws, wl) {
    list(
      impact_fit   = NULL,
      window_short = ws,
      window_long  = wl,
      tier         = "bridge",
      n_obs        = n_obs,
      anchors      = dplyr::mutate(obs, S = .data$I)
    )
  }

  ws0 <- api_windows_short[1L]
  wl0 <- api_windows_long[length(api_windows_long)]
  if (n_obs < min_obs_model) return(flat(ws0, wl0))

  ## AIC window selection (season-blind: 2-smooth model on the two indices)
  fit_one <- function(ws, wl) {
    df_m <- tibble::tibble(
      I           = obs$I,
      hydro_short = .compute_hydro_features(hydro, obs$date, ws, wl, hydro_type)$hydro_short,
      hydro_long  = .compute_hydro_features(hydro, obs$date, ws, wl, hydro_type)$hydro_long
    )
    k_h <- min(4L, n_obs - 2L)
    fit <- tryCatch(
      mgcv::gam(I ~ s(hydro_short, k = k_h) + s(hydro_long, k = k_h),
                data = df_m, method = "REML"),
      error = function(e) NULL, warning = function(w) NULL
    )
    list(fit = fit, df = df_m)
  }

  if (auto_select) {
    grid <- expand.grid(ws = api_windows_short, wl = api_windows_long)
    grid <- grid[grid$ws < grid$wl, , drop = FALSE]
    if (nrow(grid) == 0L) grid <- data.frame(ws = ws0, wl = wl0)
  } else {
    grid <- data.frame(ws = ws0, wl = wl0)
  }

  best <- NULL; best_aic <- Inf; best_ws <- ws0; best_wl <- wl0
  for (i in seq_len(nrow(grid))) {
    r <- fit_one(grid$ws[i], grid$wl[i])
    if (is.null(r$fit)) next
    a <- stats::AIC(r$fit)
    if (is.finite(a) && a < best_aic) {
      best_aic <- a; best <- r; best_ws <- grid$ws[i]; best_wl <- grid$wl[i]
    }
  }

  if (is.null(best)) return(flat(ws0, wl0))

  ## Null (intercept-only) — the model must beat it to earn the "model" tier
  null_fit <- tryCatch(
    mgcv::gam(I ~ 1, data = best$df, method = "REML"),
    error = function(e) NULL
  )
  null_aic <- if (!is.null(null_fit)) stats::AIC(null_fit) else Inf
  if (best_aic >= null_aic) return(flat(best_ws, best_wl))

  ## Recompute anchors' hydro features at the chosen windows; store residual S
  feats <- .compute_hydro_features(hydro, obs$date, best_ws, best_wl, hydro_type)
  anchors <- obs
  anchors$hydro_short <- feats$hydro_short
  anchors$hydro_long  <- feats$hydro_long
  fitted_I <- as.numeric(stats::predict(best$fit, newdata = dplyr::tibble(
    hydro_short = feats$hydro_short, hydro_long = feats$hydro_long
  )))
  anchors$S <- anchors$I - fitted_I

  list(
    impact_fit   = best$fit,
    window_short = best_ws,
    window_long  = best_wl,
    tier         = "model",
    n_obs        = n_obs,
    anchors      = anchors
  )
}


## ============================================================================
## Internal: resolve the impact at query dates
## ============================================================================

#' Hydrology-weighted bracketing-anchor interpolation of the residual state
#'
#' Pinches to the observed residual at each anchor; between two bracketing
#' anchors it leans toward the one whose antecedent hydrology more closely
#' matches the query day (and, secondarily, the nearer in time).  Outside the
#' anchor span it carries the nearest anchor's residual forward/backward (flat).
#'
#' @param anchors Tibble with `date`, `S`, `hydro_short`, `hydro_long`
#'   (≥1 row).
#' @param qdate Query Date scalar.
#' @param qshort,qlong Query-day hydro features.
#' @return Numeric interpolated residual `S`.
#' @keywords internal
.interp_residual <- function(anchors, qdate, qshort, qlong) {
  if (nrow(anchors) == 0L) return(0)
  if (nrow(anchors) == 1L) return(anchors$S[1L])

  ## Pinch exactly to an observation on its own date (residual bridge anchors).
  exact <- which(anchors$date == qdate)
  if (length(exact)) return(anchors$S[exact[1L]])

  anchors <- anchors[order(anchors$date), ]
  ## Standardise hydro features by anchor spread for a scale-free distance
  s_sd <- stats::sd(anchors$hydro_short, na.rm = TRUE); if (!is.finite(s_sd) || s_sd == 0) s_sd <- 1
  l_sd <- stats::sd(anchors$hydro_long,  na.rm = TRUE); if (!is.finite(l_sd) || l_sd == 0) l_sd <- 1
  tau_t <- stats::median(diff(as.numeric(anchors$date)), na.rm = TRUE)
  if (!is.finite(tau_t) || tau_t <= 0) tau_t <- 1

  prev_idx <- which(anchors$date <= qdate)
  next_idx <- which(anchors$date >  qdate)

  if (length(prev_idx) == 0L) return(anchors$S[next_idx[1L]])         # before first
  if (length(next_idx) == 0L) return(anchors$S[prev_idx[length(prev_idx)]])  # after last

  ip <- prev_idx[length(prev_idx)]
  inx <- next_idx[1L]

  hydro_dist <- function(i) {
    sqrt(((qshort - anchors$hydro_short[i]) / s_sd)^2 +
         ((qlong  - anchors$hydro_long[i])  / l_sd)^2)
  }
  wt <- function(i) {
    dt <- abs(as.numeric(qdate - anchors$date[i])) / tau_t
    exp(-dt - 0.5 * hydro_dist(i)^2)
  }
  wp <- wt(ip); wn <- wt(inx)
  if (wp + wn <= 0) return(0.5 * (anchors$S[ip] + anchors$S[inx]))
  (wp * anchors$S[ip] + wn * anchors$S[inx]) / (wp + wn)
}

#' Predict the site impact (and implied normalised concentration) at query dates
#'
#' For each query (date × analyte): `I_hat = beta·f(hydro) + S_interp`, with
#' `beta·f(hydro)` from the fitted response (0 for bridge-tier analytes) and
#' `S_interp` from [.interp_residual()].  Also returns `ref_norm` (from the
#' reference model) and the implied `C_norm = max(ref_norm + I_hat, 0)`.
#'
#' @param target_model A `target_model`.
#' @param query Tibble with `date` (Date) — the days to predict.
#' @param analytes Character; analytes to predict (default: all modelled).
#' @return Tibble `(date, analyte, ref_norm, impact, C_norm, impact_tier)`.
#' @keywords internal
.resolve_target_impact <- function(target_model, query, analytes = NULL) {
  qdates <- sort(unique(as.Date(query$date)))
  if (is.null(analytes)) analytes <- names(target_model$models)
  analytes <- intersect(analytes, names(target_model$models))
  if (length(analytes) == 0L || length(qdates) == 0L) {
    return(tibble::tibble(
      date = as.Date(character()), analyte = character(),
      ref_norm = numeric(), impact = numeric(), C_norm = numeric(),
      impact_tier = character()
    ))
  }

  ## ref_norm at the query dates (synthetic per-date "samples" for the resolver)
  ref_q <- .resolve_ref_norm_instant(
    target_model$reference_model,
    tibble::tibble(sample_id = as.character(qdates), datetime = qdates)
  ) |>
    dplyr::mutate(date = as.Date(.data$sample_id))

  out <- vector("list", length(analytes))
  for (j in seq_along(analytes)) {
    nm <- analytes[j]
    m  <- target_model$models[[nm]]

    feats <- .compute_hydro_features(
      target_model$hydro, qdates, m$window_short, m$window_long,
      target_model$hydro_type
    )

    hydro_pred <- if (m$tier == "model" && !is.null(m$impact_fit)) {
      as.numeric(stats::predict(m$impact_fit, newdata = dplyr::tibble(
        hydro_short = feats$hydro_short, hydro_long = feats$hydro_long
      )))
    } else rep(0, length(qdates))

    s_interp <- vapply(seq_along(qdates), function(i) {
      .interp_residual(m$anchors, qdates[i], feats$hydro_short[i], feats$hydro_long[i])
    }, numeric(1L))

    impact <- hydro_pred + s_interp

    ref_norm_nm <- ref_q$ref_norm[ref_q$analyte == nm]
    ref_dates_nm <- ref_q$date[ref_q$analyte == nm]
    ref_lookup <- stats::setNames(ref_norm_nm, as.character(ref_dates_nm))
    ref_vec <- as.numeric(ref_lookup[as.character(qdates)])
    ref_vec[is.na(ref_vec)] <- 0

    out[[j]] <- tibble::tibble(
      date        = qdates,
      analyte     = nm,
      ref_norm    = ref_vec,
      impact      = impact,
      C_norm      = pmax(ref_vec + impact, 0),
      impact_tier = m$tier
    )
  }
  dplyr::bind_rows(out)
}


## ============================================================================
## print method
## ============================================================================

#' @export
print.target_model <- function(x, ...) {
  n_an    <- length(x$models)
  n_model <- sum(vapply(x$models, function(m) m$tier == "model", logical(1L)))
  n_bridge <- n_an - n_model

  hr <- if (!is.null(x$hydro) && nrow(x$hydro) > 0L) {
    sprintf("%s – %s", min(x$hydro$date), max(x$hydro$date))
  } else "unknown"

  cat(sprintf(
    "<target_model>  fitted %s | %d analyte%s | season-blind | hydro: %s (%s)\n",
    x$fit_date, n_an, if (n_an == 1L) "" else "s", x$hydro_type, hr
  ))
  if (n_model > 0L) {
    nms <- names(Filter(function(m) m$tier == "model", x$models))
    detail <- vapply(nms, function(nm) {
      m <- x$models[[nm]]
      sprintf("%s (w=%d/%dd, n=%d)", nm, m$window_short, m$window_long, m$n_obs)
    }, character(1L))
    cat(sprintf("  hydro-response (%d):  %s\n", n_model, paste(detail, collapse = ", ")))
  }
  if (n_bridge > 0L) {
    nms <- names(Filter(function(m) m$tier == "bridge", x$models))
    cat(sprintf("  bridge-only (%d):  %s\n", n_bridge, paste(nms, collapse = ", ")))
  }
  invisible(x)
}
