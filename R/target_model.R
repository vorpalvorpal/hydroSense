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
## So the impact model is season-blind -- hydrology enters only as a *modulator*
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
## --------------
##   fit_target_model()      fit per-analyte impact-residual models
##   print.target_model()    S3 print method
##
## Internal
## --------
##   .fit_impact_response()  season-blind GAM I ~ s(hydro_short) + s(hydro_long)
##   .resolve_target_impact() predict impact (+ implied C_norm) at query dates
##   .resolve_target_impact() combines ref + beta.f(hydro) + smoothed residual


## ============================================================================
## fit_target_model()
## ============================================================================

#' Fit a season-blind predictive model of the site impact
#'
#' Models the anthropogenic increment `I = C_norm - ref_norm` (the ARA
#' "added risk", i.e. `ara_summary()`'s `C_excess`) at a target site as a
#' function of hydrology and a persistent latent state -- **never** of
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
#' (Zhang & Hirsch 2019); `f(hydro)` follows concentration--discharge theory
#' (Godsey, Kirchner & Clow 2009).
#'
#' @param target Long-format target chemistry. Required columns: `sample_id`,
#'   `datetime`, `analyte`, `value`, `detected`. Toxicants must be in ug/L;
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
#' @param method SSD method (`"multi"` or `"anzecc"`) used to derive each
#'   analyte's HC5 transform scale when `analyte_c` is not supplied.
#' @param guideline_dir Path to the ANZG guideline data folder (for the SSD
#'   fits); falls back to `getOption("leachatetools.guideline_dir")`.
#' @param transform `"pseudo_log"` (default) or `"additive"`. Controls the
#'   variance-stabilising transform applied to the impact residual before
#'   smoothing. `"pseudo_log"` uses `g = asinh(I / c)` with per-analyte scale
#'   `c = HC5` (issue #15), which compresses the dynamic range and prevents
#'   event spikes from inflating the baseline draw spread. `"additive"` keeps
#'   `g = I` (pre-#15 behaviour): the smoother operates in the original additive
#'   impact space. Ignored when `analyte_c` is supplied directly.
#' @param analyte_c Optional named numeric vector of per-analyte transform
#'   scales `c` (SSD HC5; issue #15). When `NULL` (default) it is computed from
#'   the fitted SSDs. The impact residual is smoothed on the variance-stabilising
#'   scale `g = asinh(I / c)`; an analyte with `NA`/absent `c` keeps the additive
#'   model.
#' @param api_tau_bounds_short,api_tau_bounds_long Length-2 numeric `c(lo, hi)`
#'   search ranges (days) for the short/long reservoir recession constants of
#'   `f(hydro)`, selected per analyte by profiled AIC (with a
#'   `tau_long >= 1.5*tau_short` separation).  A degenerate `c(x, x)` fixes that
#'   store's tau.  Defaults `c(1, 30)` and `c(20, 365)`.
#' @param auto_select Logical; profiled-AIC tau selection per analyte (default
#'   `TRUE`).  If `FALSE`, use the parsimonious defaults (tau 7 / 60 days).
#' @param min_obs_model Integer; minimum impact anchors required to attempt the
#'   `f(hydro)` GAM. Below this, the analyte uses the bridge tier. Default `12L`.
#' @param pool Logical (default `TRUE`). When `TRUE`, the per-analyte hydro
#'   responses are **partially pooled**: a single factor-smooth GAM
#'   (`bs = "fs"`) is fitted across all sufficiently-sampled analytes at one
#'   common AIC-selected window, shrinking each analyte's response toward a
#'   shared shape. This *regularises* noisy, low-signal analytes (it does not
#'   add hydrological coverage -- co-sampled analytes already share the same
#'   regimes), and falls back to independent fits if it fails or doesn't beat an
#'   analyte-intercept null. Set `pool = FALSE` to force independent per-analyte
#'   fits (appropriate only when all analytes are densely sampled).
#' @param eps Small positive guard. Default `1e-9`.
#'
#' @return An object of class `target_model`:
#'   \describe{
#'     \item{`$models`}{Named per-analyte list: `impact_fit` (gam or `NULL`),
#'       `tau_short`, `tau_long`, `tier` (`"model"` or `"bridge"`),
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
#' Zhang Q, Hirsch RM (2019) Water Resources Research 55(11):9705--9723.
#' Godsey SE, Kirchner JW, Clow DW (2009) Hydrological Processes 23:1844--1864.
#'
#' @examples
#' \dontrun{
#' ref_model <- fit_reference_model(reference_chem,
#'   latitude = -33.8,
#'   longitude = 151.2, conc_units = "ug/L"
#' )
#' tgt_model <- fit_target_model(target_chem, ref_model, conc_units = "ug/L")
#' tgt_model
#' }
#' @export
fit_target_model <- function(
  target,
  reference_model,
  hydro = NULL,
  hydro_type = "rainfall",
  imputation_model = NULL,
  conc_units = NULL,
  analyte_metadata = NULL,
  method = c("multi", "anzecc"),
  guideline_dir = getOption("leachatetools.guideline_dir"),
  transform = c("pseudo_log", "additive"),
  analyte_c = NULL,
  api_tau_bounds_short = c(1, 30),
  api_tau_bounds_long = c(20, 365),
  auto_select = TRUE,
  min_obs_model = 12L,
  pool = TRUE,
  eps = 1e-9
) {
  ## -- Validation -------------------------------------------------------------
  checkmate::assert_data_frame(target)
  checkmate::assert_flag(pool)
  checkmate::assert_names(names(target),
    must.include = c("sample_id", "datetime", "analyte", "value", "detected")
  )
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

  ## -- Hydrology (default: reuse the reference model's series) -----------------
  if (is.null(hydro)) {
    hydro <- reference_model$hydro
    hydro_type <- reference_model$hydro_type
  } else {
    checkmate::assert_names(names(hydro), must.include = c("date", "value"))
    if ("type" %in% names(hydro)) hydro_type <- unique(hydro$type)[1L]
    hydro_type <- match.arg(hydro_type, c("rainfall", "stage", "discharge"))
    hydro <- dplyr::select(hydro, date = "date", value = "value")
  }
  hydro <- dplyr::mutate(hydro,
    date = as.Date(.data$date),
    value = as.numeric(.data$value)
  )
  hydro <- hydro[order(hydro$date), ]

  ## -- Units, optional impute-first, normalisation ----------------------------
  meta <- .load_analyte_metadata(analyte_metadata)
  ssd_analytes <- meta$analyte[
    !is.na(meta$ssd_available) & meta$ssd_available == TRUE &
      !meta$analyte %in% .AMSPAF_EXCLUDED_ANALYTES
  ]
  target <- .convert_df_tox_to_ugL(target, ssd_analytes, conc_units, "target")

  ## Per-analyte variance-stabilising scale c = SSD HC5 (issue #15). The impact
  ## residual is smoothed on the g = asinh(I / c) scale so the process variance
  ## is proportional to concentration; c is the additive->proportional crossover.
  ## Computed once from the fitted SSDs (NA where an HC5 is unavailable -> that
  ## analyte keeps the additive model). The caller may supply `analyte_c`
  ## directly (e.g. amspaf_daily, which already derives the SSD params).
  transform <- match.arg(transform)
  method <- match.arg(method)
  if (is.null(analyte_c)) {
    if (transform == "pseudo_log") {
      ssd_p <- suppressMessages(
        derive_ssd_params(meta, method = method, guideline_dir = guideline_dir)
      )
      analyte_c <- stats::setNames(
        vapply(ssd_p$fit, function(f) {
          tryCatch(.analyte_c(f), error = function(e) NA_real_)
        }, numeric(1L)),
        ssd_p$analyte
      )
    } else {
      ## additive: g = I for every analyte. An empty named numeric means every
      ## analyte_c[nm] lookup returns NA_real_, routing through the existing
      ## additive branch (the pre-#15 model) throughout.
      analyte_c <- numeric(0L)
    }
  }

  ## Impute-first only when the model can actually impute (has fitted groups).
  ## A PCA-only imputation_model still feeds the WQ layer below via its $pca.
  if (!is.null(imputation_model) &&
    length(imputation_model$groups %||% list()) > 0L) {
    if (!"site_id" %in% names(target)) target$site_id <- "target"
    target <- impute_chemistry(target, imputation_model)
  }

  target <- dplyr::mutate(target, .date = as.Date(.data$datetime))
  norm_all <- .normalise_ref_observations(target, target, meta)
  norm_det <- dplyr::filter(
    norm_all, .data$detected, !is.na(.data$value_norm),
    .data$value_norm > 0, .data$analyte %in% ssd_analytes
  )

  ## -- Reference norm at each target sample (instant resolver) ----------------
  ref_resolved <- .resolve_ref_norm_instant(
    reference_model,
    dplyr::distinct(target, .data$sample_id, .data$datetime)
  )

  ## -- Per-sample impact: I = value_norm - ref_norm ---------------------------
  obs_samples <- norm_det |>
    dplyr::select("sample_id", date = ".date", "analyte", "value_norm") |>
    dplyr::inner_join(
      dplyr::select(ref_resolved, "sample_id", "analyte", "ref_norm"),
      by = c("sample_id", "analyte")
    ) |>
    dplyr::mutate(I = .data$value_norm - .data$ref_norm) |>
    dplyr::filter(is.finite(.data$I))

  ## WQ-layer predictor scores (issue #14 item B). The WQ->metal prediction is a
  ## regression on the imputation model's chemistry PCA -- not Bayesian (the
  ## cross-metal coupling only acts when sibling metals are observed). We reuse
  ## the PCA to predict each metal from water quality and interpolate only the
  ## residual `d`, which lets WQ-only days beat pure impact interpolation.
  pca <- if (!is.null(imputation_model)) imputation_model$pca else NULL
  pc_cols <- character(0)
  if (!is.null(pca)) {
    pc_scores <- .compute_pca_scores(target, pca)
    pc_cols <- grep("^PC", names(pc_scores), value = TRUE)
    obs_samples <- dplyr::left_join(obs_samples, pc_scores, by = "sample_id")
  }

  ## Date-aggregated impact anchors (average duplicate same-day grabs)
  anchors_all <- obs_samples |>
    dplyr::summarise(I = mean(.data$I, na.rm = TRUE), .by = c("date", "analyte"))

  ## Variance-stabilising transform of the impact, per analyte (issue #15):
  ## g = asinh(I / c). The smoother bridges `g`; analytes without an HC5 keep
  ## g = I (additive). Same-day grabs are averaged in impact space first (the
  ## anchor is the day's impact level), then transformed.
  anchors_all <- anchors_all |>
    dplyr::mutate(scale_c = unname(analyte_c[.data$analyte])) |>
    dplyr::group_by(.data$analyte) |>
    dplyr::mutate(g = {
      cc <- dplyr::first(.data$scale_c)
      if (is.na(cc)) .data$I else .g_transform(.data$I, cc)
    }) |>
    dplyr::ungroup()

  ## -- Hydro response: pooled (factor-smooth shrinkage) or per-analyte ---------
  target_analytes <- intersect(unique(anchors_all$analyte), ssd_analytes)

  if (pool && length(target_analytes) >= 2L) {
    base_models <- .fit_pooled_impact_response(
      anchors_all, target_analytes, hydro, hydro_type,
      api_tau_bounds_short, api_tau_bounds_long, auto_select, min_obs_model, eps
    )
  } else {
    base_models <- stats::setNames(
      lapply(target_analytes, function(nm) {
        obs <- anchors_all |>
          dplyr::filter(.data$analyte == .env$nm) |>
          dplyr::arrange(.data$date)
        feats <- .compute_hydro_features(
          hydro, obs$date, max(api_tau_bounds_short), max(api_tau_bounds_long),
          hydro_type
        )
        obs <- dplyr::bind_cols(obs, dplyr::select(feats, -"date"))
        .fit_impact_response(
          obs, hydro, hydro_type, api_tau_bounds_short,
          api_tau_bounds_long, auto_select, min_obs_model, eps
        )
      }),
      target_analytes
    )
  }

  ## -- WQ layer + residual d per analyte (only when a PCA is available) --------
  models <- vector("list", length(target_analytes))
  names(models) <- target_analytes
  for (nm in target_analytes) {
    fit_res <- base_models[[nm]]
    fit_res$wq_fit <- NULL
    fit_res$d_anchors <- NULL
    fit_res$scale_c <- unname(analyte_c[nm]) # transform scale (NA -> additive)
    if (length(pc_cols) > 0L) {
      os_nm <- dplyr::filter(obs_samples, .data$analyte == .env$nm)
      fit_res <- c(fit_res, .fit_wq_layer(
        os_nm, pc_cols, hydro, hydro_type,
        fit_res$tau_short, fit_res$tau_long, min_obs_model,
        scale_c = fit_res$scale_c
      ))
    }
    models[[nm]] <- fit_res
  }

  structure(
    list(
      models          = models,
      reference_model = reference_model,
      hydro           = hydro,
      hydro_type      = hydro_type,
      pca             = pca,
      pc_cols         = pc_cols,
      fit_date        = Sys.Date()
    ),
    class = "target_model"
  )
}

#' Fit the WQ->metal layer and its residual `d` for one analyte
#'
#' A GAM of normalised concentration on the chemistry-PCA scores (the
#' non-Bayesian WQ prediction), kept only if it beats an intercept-only null by
#' AIC. The residual `d = value_norm - WQ-prediction` is date-aggregated with
#' hydro features for bracketing-bridge interpolation, exactly like the impact
#' state `S`.
#'
#' @return List with `wq_fit` (gam or `NULL`) and `d_anchors` (or `NULL`).
#' @keywords internal
.fit_wq_layer <- function(os_nm, pc_cols, hydro, hydro_type,
                          tau_short, tau_long, min_obs_model,
                          scale_c = NA_real_) {
  os_nm <- os_nm[stats::complete.cases(os_nm[, pc_cols, drop = FALSE]), , drop = FALSE]
  if (nrow(os_nm) < min_obs_model) {
    return(list(wq_fit = NULL, d_anchors = NULL))
  }

  ## Model the WQ->concentration prediction on the variance-stabilising scale
  ## (issue #15): gvn = asinh(value_norm / c) (or value_norm when c is NA). The
  ## residual `d` is then in g-space, like the impact residual `S`.
  os_nm$gvn <- if (is.na(scale_c)) {
    os_nm$value_norm
  } else {
    .g_transform(os_nm$value_norm, scale_c)
  }

  k_pc <- min(4L, nrow(os_nm) - 2L)
  form <- stats::as.formula(paste(
    "gvn ~", paste(sprintf("s(%s, k = %d)", pc_cols, k_pc), collapse = " + ")
  ))
  wq_gam <- tryCatch(mgcv::gam(form, data = os_nm, method = "REML"),
    error = function(e) NULL, warning = function(w) NULL
  )
  if (is.null(wq_gam)) {
    return(list(wq_fit = NULL, d_anchors = NULL))
  }

  null_gam <- tryCatch(mgcv::gam(gvn ~ 1, data = os_nm, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(null_gam) || stats::AIC(wq_gam) >= stats::AIC(null_gam)) {
    return(list(wq_fit = NULL, d_anchors = NULL))
  }

  os_nm$d <- os_nm$gvn - as.numeric(stats::predict(wq_gam))
  d_anch <- os_nm |>
    dplyr::summarise(S = mean(.data$d, na.rm = TRUE), .by = "date") |>
    dplyr::arrange(.data$date)
  feats <- .compute_hydro_features(hydro, d_anch$date, tau_short, tau_long, hydro_type)
  d_anch$hydro_short <- feats$hydro_short
  d_anch$hydro_long <- feats$hydro_long

  list(wq_fit = wq_gam, d_anchors = d_anch)
}


## ============================================================================
## Internal: season-blind hydrological-response fit
## ============================================================================

#' Fit the season-blind impact response for one analyte
#'
#' Selects the short/long reservoir recession constants (tau) by profiled AIC,
#' fits `I ~ s(hydro_short) + s(hydro_long)` (no seasonal term), and keeps it
#' only if it beats the intercept-only null.  Stores the de-trended residual
#' `S = I - fitted` as the bridge-interpolation state.
#'
#' @return List with `impact_fit`, `tau_short`, `tau_long`, `tier`,
#'   `n_obs`, `anchors` (`date`, `I`, `S`, `hydro_short`, `hydro_long`).
#' @keywords internal
.fit_impact_response <- function(obs, hydro, hydro_type,
                                 tau_bounds_short, tau_bounds_long,
                                 auto_select, min_obs_model, eps) {
  n_obs <- nrow(obs)

  flat <- function(ts, tl) {
    list(
      impact_fit = NULL,
      tau_short  = ts,
      tau_long   = tl,
      tier       = "bridge",
      n_obs      = n_obs,
      anchors    = dplyr::mutate(obs, S = .data$g)
    )
  }

  ts0 <- .REF_TAU_DEFAULT_SHORT
  tl0 <- .REF_TAU_DEFAULT_LONG
  if (n_obs < min_obs_model) {
    return(flat(ts0, tl0))
  }

  ## Profiled-AIC tau selection (season-blind 2-smooth model on the indices)
  fit_one <- function(ts, tl) {
    feats <- .compute_hydro_features(hydro, obs$date, ts, tl, hydro_type)
    df_m <- tibble::tibble(
      g           = obs$g, # variance-stabilised impact (issue #15)
      hydro_short = feats$hydro_short,
      hydro_long  = feats$hydro_long
    )
    k_h <- min(4L, n_obs - 2L)
    fit <- tryCatch(
      mgcv::gam(g ~ s(hydro_short, k = k_h) + s(hydro_long, k = k_h),
        data = df_m, method = "REML"
      ),
      error = function(e) NULL, warning = function(w) NULL
    )
    list(fit = fit, df = df_m)
  }
  aic_fn <- function(ts, tl) {
    f <- fit_one(ts, tl)$fit
    if (is.null(f)) Inf else stats::AIC(f)
  }

  if (auto_select) {
    tp <- .optimise_tau_pair(aic_fn, tau_bounds_short, tau_bounds_long)
    best_ts <- tp$tau_short
    best_tl <- tp$tau_long
  } else {
    best_ts <- ts0
    best_tl <- tl0
  }
  best <- fit_one(best_ts, best_tl)
  if (is.null(best$fit)) {
    return(flat(ts0, tl0))
  }
  best_aic <- stats::AIC(best$fit)

  ## Null (intercept-only) -- the model must beat it to earn the "model" tier
  null_fit <- tryCatch(
    mgcv::gam(g ~ 1, data = best$df, method = "REML"),
    error = function(e) NULL
  )
  null_aic <- if (!is.null(null_fit)) stats::AIC(null_fit) else Inf
  if (best_aic >= null_aic) {
    return(flat(best_ts, best_tl))
  }

  ## Recompute anchors' hydro features at the chosen tau; store residual S
  ## in the transformed (g) space (issue #15).
  feats <- .compute_hydro_features(hydro, obs$date, best_ts, best_tl, hydro_type)
  anchors <- obs
  anchors$hydro_short <- feats$hydro_short
  anchors$hydro_long <- feats$hydro_long
  fitted_g <- as.numeric(stats::predict(best$fit, newdata = dplyr::tibble(
    hydro_short = feats$hydro_short, hydro_long = feats$hydro_long
  )))
  anchors$S <- anchors$g - fitted_g

  list(
    impact_fit = best$fit,
    tau_short  = best_ts,
    tau_long   = best_tl,
    tier       = "model",
    n_obs      = n_obs,
    anchors    = anchors
  )
}

#' Pooled (hierarchical) season-blind impact response across analytes
#'
#' Fits a single factor-smooth GAM
#' `z ~ s(hydro_short, analyte, bs="fs") + s(hydro_long, analyte, bs="fs")`
#' at one AIC-selected common window pair, where `z = (I - mu_a) / sd_a` is the
#' impact standardised **per analyte**. Pooling on this common unit scale shares
#' the response *shape* across analytes while leaving each analyte's *magnitude*
#' untouched (it is restored by de-standardising the fitted shape,
#' `fitted_I = mu_a + sd_a * z_hat`). This is essential: the `bs = "fs"` penalty
#' shrinks each analyte's level as well as its wiggliness, so pooling the raw
#' `I` would drag a large-signal analyte (e.g. Cu) toward a population dominated
#' by near-zero ones and inflate the near-zero ones in turn.
#'
#' Pooling regularises noisy, low-SNR analytes by borrowing a response shape
#' from co-varying ones -- it does **not** add hydrological coverage (co-sampled
#' analytes share the same regimes). Analytes with fewer than `min_obs_model`
#' anchors, or with ~no impact variance (no shape to share), get a per-analyte
#' flat bridge; if pooling fails or doesn't beat the no-shared-shape null
#' (`z ~ 1`), the poolable analytes fall back to independent fits.
#'
#' @return Named list (per analyte) of `.fit_impact_response()`-shaped objects;
#'   pooled analytes additionally carry `pooled = TRUE`, `analyte`,
#'   `pool_levels`, and the de-standardisation pair `pool_center` (`mu_a`) and
#'   `pool_scale` (`sd_a`) for prediction.
#' @keywords internal
.fit_pooled_impact_response <- function(anchors_all, target_analytes, hydro,
                                        hydro_type, tau_bounds_short,
                                        tau_bounds_long, auto_select,
                                        min_obs_model, eps) {
  anchors_for <- function(nm, ts, tl) {
    obs <- anchors_all |>
      dplyr::filter(.data$analyte == .env$nm) |>
      dplyr::arrange(.data$date)
    f <- .compute_hydro_features(hydro, obs$date, ts, tl, hydro_type)
    obs$hydro_short <- f$hydro_short
    obs$hydro_long <- f$hydro_long
    obs
  }
  per_analyte_fit <- function(nm) {
    obs <- anchors_for(nm, max(tau_bounds_short), max(tau_bounds_long))
    .fit_impact_response(
      obs, hydro, hydro_type, tau_bounds_short,
      tau_bounds_long, auto_select, min_obs_model, eps
    )
  }

  bridge_for <- function(nm, ts, tl) {
    obs <- anchors_for(nm, ts, tl)
    list(
      impact_fit = NULL, tau_short = ts, tau_long = tl, tier = "bridge",
      n_obs = nrow(obs), anchors = dplyr::mutate(obs, S = .data$g)
    )
  }

  counts <- anchors_all |> dplyr::count(.data$analyte)
  poolable <- intersect(
    target_analytes,
    counts$analyte[counts$n >= min_obs_model]
  )
  small <- setdiff(target_analytes, poolable)

  out <- list()
  for (nm in small) { # too few anchors -> flat bridge
    out[[nm]] <- bridge_for(nm, .REF_TAU_DEFAULT_SHORT, .REF_TAU_DEFAULT_LONG)
  }

  ## Per-analyte standardisation of the impact `I` BEFORE pooling.  `I` is
  ## window-independent (windows change only the hydro features), so the
  ## centre/scale are computed once.  Pooling the SHAPE on a common z-scale is
  ## what keeps a large-magnitude analyte (e.g. Cu) from being shrunk toward a
  ## population dominated by near-zero analytes -- and a near-zero analyte
  ## (e.g. Ni) from being inflated toward one carrying the big signal.  Without
  ## this, the shared `bs = "fs"` penalty (which shrinks each analyte's level,
  ## not just its wiggliness) cross-contaminates impact magnitudes across
  ## chemically unrelated analytes.
  ## Standardise on the transformed scale g (issue #15), not raw impact I.
  zstats <- anchors_all |>
    dplyr::filter(.data$analyte %in% .env$poolable) |>
    dplyr::summarise(
      mu = mean(.data$g, na.rm = TRUE),
      sd = stats::sd(.data$g, na.rm = TRUE), .by = "analyte"
    )
  ## Analytes with ~no impact variance carry no shape to share -> flat bridge.
  flat_nm <- zstats$analyte[!is.finite(zstats$sd) | zstats$sd < eps]
  poolable <- setdiff(poolable, flat_nm)
  for (nm in flat_nm) {
    out[[nm]] <- bridge_for(nm, .REF_TAU_DEFAULT_SHORT, .REF_TAU_DEFAULT_LONG)
  }

  if (length(poolable) < 2L) { # nothing to pool
    for (nm in poolable) out[[nm]] <- per_analyte_fit(nm)
    return(out[target_analytes])
  }

  zsd <- function(nm) zstats$sd[zstats$analyte == nm]
  zmu <- function(nm) zstats$mu[zstats$analyte == nm]

  ## Profiled-AIC select a common (short, long) tau for the pooled factor-smooth
  ## fit, fitted on the standardised response z = (I - mu_a) / sd_a.
  build_df <- function(ts, tl) {
    purrr::map_dfr(poolable, function(nm) {
      o <- anchors_for(nm, ts, tl)
      tibble::tibble(
        z = (o$g - zmu(nm)) / zsd(nm), hydro_short = o$hydro_short,
        hydro_long = o$hydro_long, analyte = nm
      )
    }) |> dplyr::mutate(analyte = factor(.data$analyte))
  }
  fit_pool <- function(df_m) {
    nlev <- nlevels(df_m$analyte)
    k_h <- max(3L, min(4L, floor(nrow(df_m) / nlev) - 1L))
    ## Factor-smooth fits routinely emit benign convergence warnings -- suppress
    ## them but keep the fit (only a hard error means "unusable").
    tryCatch(
      suppressWarnings(
        mgcv::gam(
          z ~ s(hydro_short, analyte, bs = "fs", k = k_h) +
            s(hydro_long, analyte, bs = "fs", k = k_h),
          data = df_m, method = "REML"
        )
      ),
      error = function(e) NULL
    )
  }
  pool_aic <- function(ts, tl) {
    f <- fit_pool(build_df(ts, tl))
    if (is.null(f)) Inf else stats::AIC(f)
  }
  if (auto_select) {
    tp <- .optimise_tau_pair(pool_aic, tau_bounds_short, tau_bounds_long)
    best_ws <- tp$tau_short
    best_wl <- tp$tau_long
  } else {
    best_ws <- .REF_TAU_DEFAULT_SHORT
    best_wl <- .REF_TAU_DEFAULT_LONG
  }
  best_df <- build_df(best_ws, best_wl)
  best <- fit_pool(best_df)
  if (is.null(best)) { # pooled fit failed -> independent
    for (nm in poolable) out[[nm]] <- per_analyte_fit(nm)
    return(out[target_analytes])
  }
  best_aic <- stats::AIC(best)

  ## Null on the SAME (z) scale: "no shared hydro shape, every analyte flat at
  ## its own mean" (z is de-meaned per analyte, so the null is z ~ 1).
  null_fit <- tryCatch(mgcv::gam(z ~ 1, data = best_df, method = "REML"),
    error = function(e) NULL
  )
  pooled_useful <- !is.null(null_fit) && best_aic < stats::AIC(null_fit)
  lev <- levels(best_df$analyte)

  for (nm in poolable) {
    obs <- anchors_for(nm, best_ws, best_wl)
    if (pooled_useful) {
      nd <- tibble::tibble(
        hydro_short = obs$hydro_short,
        hydro_long = obs$hydro_long,
        analyte = factor(nm, levels = lev)
      )
      ## De-standardise the shared shape back to this analyte's own magnitude.
      z_hat <- as.numeric(stats::predict(best, newdata = nd))
      fitted_g <- zmu(nm) + zsd(nm) * z_hat # de-standardise in g-space
      out[[nm]] <- list(
        impact_fit = best, tau_short = best_ws, tau_long = best_wl,
        tier = "model", n_obs = nrow(obs),
        anchors = dplyr::mutate(obs, S = .data$g - fitted_g),
        pooled = TRUE, analyte = nm, pool_levels = lev,
        pool_center = zmu(nm), pool_scale = zsd(nm)
      )
    } else {
      out[[nm]] <- bridge_for(nm, best_ws, best_wl)
    }
  }
  out[target_analytes]
}


## ============================================================================
## Internal: resolve the impact at query dates
## ============================================================================

#' Standardised short-window hydro feature over query dates (q-modulation input)
#'
#' Drives the hydrology-modulated process variance of the residual smoother
#' (`q_mult = exp(kappa * z)`).
#' @keywords internal
.hydro_zscore <- function(target_model, qdates, m) {
  feats <- .compute_hydro_features(
    target_model$hydro, qdates, m$tau_short,
    m$tau_long, target_model$hydro_type
  )
  zs <- feats$hydro_short
  sdv <- stats::sd(zs, na.rm = TRUE)
  if (!is.finite(sdv) || sdv == 0) sdv <- 1
  (zs - mean(zs, na.rm = TRUE)) / sdv
}

#' Deterministic season-blind hydrology prediction `beta.f(hydro)` at query dates
#' @keywords internal
.hydro_pred <- function(m, feats, qdates) {
  if (m$tier == "model" && !is.null(m$impact_fit)) {
    nd <- dplyr::tibble(
      hydro_short = feats$hydro_short,
      hydro_long = feats$hydro_long
    )
    if (isTRUE(m$pooled)) {
      nd$analyte <- factor(m$analyte, levels = m$pool_levels)
      z_hat <- as.numeric(stats::predict(m$impact_fit, newdata = nd))
      m$pool_center + m$pool_scale * z_hat
    } else {
      as.numeric(stats::predict(m$impact_fit, newdata = nd))
    }
  } else {
    rep(0, length(qdates))
  }
}

#' Build the state-space residual smoother for one analyte
#'
#' Uses the WQ residual `d` (`m$d_anchors`, WQ tier) or the impact residual `S`
#' (`m$anchors`, impact/bridge tier). The daily grid is clipped to the analyte's
#' grab span; hydrology modulates the process variance. See
#' [.residual_smoother()].
#' @keywords internal
.analyte_residual_smoother <- function(m, target_model, qdates, kappa = 0.5,
                                       scale = 1, r_vec = NULL) {
  use_wq <- !is.null(m$wq_fit) && !is.null(m$d_anchors) && nrow(m$d_anchors) >= 2L
  anch <- if (use_wq) m$d_anchors else m$anchors
  if (is.null(anch) || nrow(anch) < 1L) {
    return(NULL)
  }
  z_hydro <- .hydro_zscore(target_model, qdates, m)
  .residual_smoother(anch$date, anch$S, qdates,
    z_hydro = z_hydro,
    kappa = kappa, scale = scale, r_vec = r_vec
  )
}

#' Align a residual path (values on the smoother's clipped grid) to query dates
#'
#' Returns `NA` for query dates outside the analyte's clipped grab span -- those
#' rows are dropped by [.resolve_target_impact()] (per-analyte clipping).
#' @keywords internal
.residual_on_qdates <- function(grid_dates, path, qdates) {
  lut <- stats::setNames(path, as.character(grid_dates))
  as.numeric(lut[as.character(qdates)])
}

#' Predict the site impact (and implied normalised concentration) at query dates
#'
#' For each query (date x analyte): `I_hat = beta*f(hydro) + S`, with
#' `beta*f(hydro)` from the fitted response (0 for bridge-tier analytes) and the
#' residual `S` from the state-space smoother (`residual_paths`, or the smoother
#' posterior mean when `NULL`).  Also returns `ref_norm` (from the reference
#' model) and the implied `C_norm = max(ref_norm + I_hat, 0)`.
#'
#' @param target_model A `target_model`.
#' @param query Tibble with `date` (Date) -- the days to predict.
#' @param analytes Character; analytes to predict (default: all modelled).
#' @param wq Optional long-format water-quality data frame (`sample_id`,
#'   `analyte`, `value`) giving each query day's WQ, with `sample_id` equal to
#'   the query date as a character string. When supplied and the analyte has a
#'   fitted WQ layer, the day's concentration is predicted as
#'   `WQ-prediction + d_interp` (tier `"wq"`) instead of `ref + impact`.
#' @return Tibble `(date, analyte, ref_norm, impact, C_norm, impact_tier)`.
#' @keywords internal
.resolve_target_impact <- function(target_model, query, analytes = NULL,
                                   wq = NULL, residual_paths = NULL,
                                   kappa = 0.5, scale = 1, ref_q = NULL,
                                   static = NULL) {
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

  ## ref_norm at the query dates (synthetic per-date "samples" for the resolver).
  ## Static across draws (ARA cancels it), so it can be precomputed once in
  ## .fit_daily_target() and passed in via `ref_q`; otherwise compute it here.
  if (is.null(ref_q)) {
    ref_q <- .resolve_ref_norm_instant(
      target_model$reference_model,
      tibble::tibble(sample_id = as.character(qdates), datetime = qdates)
    ) |>
      dplyr::mutate(date = as.Date(.data$sample_id))
  }

  ## WQ-layer PC scores at the query dates (sample_id == date string)
  pc_q <- NULL
  if (!is.null(wq) && !is.null(target_model$pca)) {
    pc_q <- .compute_pca_scores(wq, target_model$pca)
  }

  out <- vector("list", length(analytes))
  for (j in seq_along(analytes)) {
    nm <- analytes[j]
    m <- target_model$models[[nm]]
    sc <- static[[nm]]
    c_a <- m$scale_c %||% NA_real_ # transform scale c (NA -> additive model)

    ## ref_vec and beta.f(hydro) are static across draws -- use the precomputed
    ## context when supplied (and compute beta.f as lpmatrix %*% coef, which
    ## equals predict.gam() but avoids rebuilding the basis every draw).
    if (!is.null(sc)) {
      feats <- list(hydro_short = sc$hydro_short, hydro_long = sc$hydro_long)
      ref_vec <- sc$ref_vec
      hp <- if (!is.null(sc$lp)) {
        z <- as.numeric(sc$lp %*% stats::coef(m$impact_fit))
        if (isTRUE(sc$pooled)) sc$pool_center + sc$pool_scale * z else z
      } else {
        rep(0, length(qdates))
      }
    } else {
      feats <- .compute_hydro_features(
        target_model$hydro, qdates, m$tau_short, m$tau_long,
        target_model$hydro_type
      )
      ref_norm_nm <- ref_q$ref_norm[ref_q$analyte == nm]
      ref_dates_nm <- ref_q$date[ref_q$analyte == nm]
      ref_lookup <- stats::setNames(ref_norm_nm, as.character(ref_dates_nm))
      ref_vec <- as.numeric(ref_lookup[as.character(qdates)])
      ref_vec[is.na(ref_vec)] <- 0
      hp <- .hydro_pred(m, feats, qdates) # beta.f(hydro), season-blind
    }

    ## Residual path on qdates (NA outside the analyte's clipped grab span).
    ## Provided by the caller (draw column or centre mean) or, if NULL, built
    ## here from the smoother posterior mean (standalone / point use).
    rp <- residual_paths[[nm]]
    if (is.null(rp)) {
      sm <- .analyte_residual_smoother(m, target_model, qdates,
        kappa = kappa,
        scale = scale
      )
      rp <- if (is.null(sm)) {
        rep(NA_real_, length(qdates))
      } else {
        .residual_on_qdates(sm$grid_dates, sm$mean, qdates)
      }
    }

    use_wq <- !is.null(pc_q) && !is.null(m$wq_fit) && !is.null(m$d_anchors)
    if (use_wq) {
      ## Tier "wq": WQ-prediction + smoothed residual d (rp). Where the WQ
      ## prediction is unavailable (NA PCs), fall back to ref + beta.f(hydro).
      pc_lookup <- dplyr::left_join(
        tibble::tibble(sample_id = as.character(qdates)), pc_q,
        by = "sample_id"
      )
      ## c_wq + rp live on the g scale (issue #15); invert to concentration.
      c_wq <- as.numeric(stats::predict(m$wq_fit, newdata = pc_lookup))
      g_c <- c_wq + rp
      c_norm <- pmax(if (is.na(c_a)) g_c else .g_inverse(g_c, c_a), 0)
      tier_vec <- rep("wq", length(qdates))
      bad <- !is.finite(c_norm) & !is.na(rp) # within span but WQ missing
      if (any(bad)) {
        hp_imp <- if (is.na(c_a)) hp[bad] else .g_inverse(hp[bad], c_a)
        c_norm[bad] <- pmax(ref_vec[bad] + hp_imp, 0)
        tier_vec[bad] <- m$tier
      }
      impact <- c_norm - ref_vec
    } else {
      ## hp (beta.f(hydro)) + rp (smoothed residual) live on the g scale;
      ## invert to recover the impact in concentration units (issue #15).
      g_hat <- hp + rp
      impact <- if (is.na(c_a)) g_hat else .g_inverse(g_hat, c_a)
      c_norm <- pmax(ref_vec + impact, 0)
      tier_vec <- rep(m$tier, length(qdates))
    }

    ## Per-analyte clip: drop query dates outside the smoother grid (rp == NA).
    keep <- !is.na(rp)
    if (!any(keep)) next
    out[[j]] <- tibble::tibble(
      date        = qdates[keep],
      analyte     = nm,
      ref_norm    = ref_vec[keep],
      impact      = impact[keep],
      C_norm      = c_norm[keep],
      impact_tier = tier_vec[keep]
    )
  }
  dplyr::bind_rows(out)
}


## ============================================================================
## print method
## ============================================================================

#' @export
print.target_model <- function(x, ...) {
  n_an <- length(x$models)
  n_model <- sum(vapply(x$models, function(m) m$tier == "model", logical(1L)))
  n_bridge <- n_an - n_model

  hr <- if (!is.null(x$hydro) && nrow(x$hydro) > 0L) {
    sprintf("%s -- %s", min(x$hydro$date), max(x$hydro$date))
  } else {
    "unknown"
  }

  cat(sprintf(
    "<target_model>  fitted %s | %d analyte%s | season-blind | hydro: %s (%s)\n",
    x$fit_date, n_an, if (n_an == 1L) "" else "s", x$hydro_type, hr
  ))
  if (n_model > 0L) {
    nms <- names(Filter(function(m) m$tier == "model", x$models))
    detail <- vapply(nms, function(nm) {
      m <- x$models[[nm]]
      sprintf("%s (τ=%.0f/%.0fd, n=%d)", nm, m$tau_short, m$tau_long, m$n_obs)
    }, character(1L))
    cat(sprintf("  hydro-response (%d):  %s\n", n_model, paste(detail, collapse = ", ")))
  }
  if (n_bridge > 0L) {
    nms <- names(Filter(function(m) m$tier == "bridge", x$models))
    cat(sprintf("  bridge-only (%d):  %s\n", n_bridge, paste(nms, collapse = ", ")))
  }
  n_wq <- sum(vapply(x$models, function(m) !is.null(m$wq_fit), logical(1L)))
  if (n_wq > 0L) {
    nms <- names(Filter(function(m) !is.null(m$wq_fit), x$models))
    cat(sprintf("  WQ layer (%d):  %s\n", n_wq, paste(nms, collapse = ", ")))
  }
  invisible(x)
}
