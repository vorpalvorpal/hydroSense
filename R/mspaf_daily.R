## ============================================================================
## mspaf_daily -- continuous daily-resolved msPAF time series
## ============================================================================
##
## The core idea: grab chemistry is sparse (bi-monthly, weekly, etc.); daily
## msPAF requires daily chemistry.  This function interpolates each analyte
## onto a fine date grid, constructs synthetic "one-per-day" samples, and
## runs the existing add_mspaf() engine on them.
##
## Two interpolation styles are offered:
##   forward_fill  -- step function; each day inherits the most recently
##                    observed value.  Conservative; never exceeds observations.
##   linear        -- piecewise-linear in concentration space (log-space for
##                    SSD-eligible toxicants; linear for co-analytes like pH
##                    and temperature).  Smoother trajectories between grabs.
##
## Critical design constraint: synthetic samples must NOT carry a focal_date
## column.  That column triggers .resolve_ref_norm_chronic() inside
## add_mspaf(), which would double-integrate the ARA reference.  We want
## .resolve_ref_norm_instant() -- pointwise matching per daily sample.

## ============================================================================
## mspaf_daily
## ============================================================================

#' Continuous daily msPAF time series from interpolated grab chemistry
#'
#' Interpolates per-analyte grab chemistry onto a daily date grid and computes
#' msPAF for every day within the requested date range.  The result is a tidy
#' tibble with one row per (site \eqn{\times} day), suitable for trend
#' analysis, visualisation, and downstream [time_weighted_aggregate()] calls.
#'
#' @section Interpolation:
#' Each analyte is interpolated independently:
#' \itemize{
#'   \item `"forward_fill"` carries the last directly observed value forward
#'     until the next grab.  Produces a step function.
#'   \item `"linear"` linearly interpolates between consecutive observations.
#'     For SSD-eligible toxicants (metals, ammonia), interpolation is performed
#'     in log-concentration space, which is more appropriate for log-normally
#'     distributed data and avoids negative intermediate values.  Co-analytes
#'     (pH, temperature, DOC, hardness) are interpolated linearly.
#'   \item `"model"` fits a season-blind [fit_target_model()] on the grab
#'     chemistry and the supplied `reference_model`, and predicts each
#'     toxicant's concentration between grabs as
#'     `reference + impact`, where the impact (the leachate-attributable
#'     increment) is modelled from hydrology and a persistent latent state but
#'     **never** from day-of-year. Co-analytes are forward-filled. Requires
#'     `reference_model`; see [fit_target_model()] for the method.
#' }
#' Below-detection values are treated as their detection-limit value for
#' interpolation purposes, matching the treatment in [add_mspaf()].
#'
#' @section Leading edge:
#' Days before the first grab sample for an analyte are outside the observation
#' record:
#' \itemize{
#'   \item `"drop"` (default) excludes such days from the output.
#'   \item `"backfill"` assigns the first observed value to all earlier days.
#'     Use cautiously -- it assumes the analyte was at its first observed level
#'     before sampling commenced.
#' }
#'
#' @section Temperature for ammonia:
#' When `NH3-N` is in `df`, water temperature is required per sample for the
#' un-ionised fraction normalisation.  Two sources are accepted:
#' \itemize{
#'   \item Temperature rows already in `df` (measured on grab-sample days) are
#'     interpolated along with other analytes.
#'   \item An external daily temperature series supplied via the `temperature`
#'     argument (e.g. from [estimate_water_temp()]) fills in days without a
#'     grab temperature.  Grab-sample-day measurements take priority.
#' }
#' Set `require_temperature = FALSE` only for datasets that do not contain
#' ammonia.
#'
#' @param df Long-format grab chemistry data frame. Required columns:
#'   `sample_id`, `site_id`, `datetime` (Date or POSIXct), `analyte`, `value`,
#'   `detected`. Optional but propagated: `units.analyte`, `imputed`.
#'   Chemistry for multiple sites may be stacked; interpolation and msPAF are
#'   computed per site.
#' @param temperature Optional daily water temperature data frame for days
#'   without a grab temperature measurement. Required columns: `datetime` (Date
#'   or POSIXct) and `value` (temperature in \eqn{{}^\circ}C). The output of
#'   [estimate_water_temp()] is accepted directly (extra columns are ignored).
#'   When both this argument and grab-sample temperature rows in `df` are
#'   present for the same day, the grab measurement takes priority.
#'   `NULL` (default) means temperature must come from `df` rows alone.
#' @param reference Background reference chemistry for ARA adjustment. Accepts
#'   the same four forms as [add_mspaf()]. Controls **only** whether background
#'   is subtracted; it is independent of `interpolation`. With
#'   `interpolation = "model"`, pass the same `reference_model` here to assess
#'   the leachate-attributable increment, or `NULL` to assess total
#'   concentration.
#' @param reference_model A `reference_model` from [fit_reference_model()].
#'   **Required** when `interpolation = "model"` — it supplies the background
#'   and catchment hydrology used by the season-blind target impact model
#'   ([fit_target_model()]) that interpolates toxicants between grabs. Ignored
#'   for the `"forward_fill"` and `"linear"` paths.
#' @param imputation_model Optional `imputation_model` from
#'   [fit_imputation_model()], passed to [fit_target_model()] for tier-2
#'   enrichment under `interpolation = "model"`. Requires **brms**.
#' @param start,end Date boundaries for the daily grid. Default: earliest and
#'   latest `datetime` values in `df`. Coerced to `Date`.
#' @param by Temporal resolution string passed to [seq.Date()]. Default
#'   `"day"`.
#' @param interpolation How to fill gaps between grab samples. One of
#'   `"forward_fill"` (default) or `"linear"`. See the *Interpolation* section.
#' @param leading_edge What to do with days before the first grab sample for
#'   an analyte. `"drop"` (default) or `"backfill"`. See the *Leading edge*
#'   section.
#' @param analyte_metadata Data frame of analyte metadata, or `NULL` to use
#'   the bundled metadata. Passed to [add_mspaf()].
#' @param method SSD method. One of `"multi"` (default, model-averaged) or
#'   `"anzecc"`. Passed to [add_mspaf()].
#' @param guideline_dir Path to the ANZG guideline data folder. Falls back to
#'   `getOption("hydroSense.guideline_dir")`.
#' @param min_analytes Minimum number of SSD-eligible analytes per day for
#'   msPAF to be computed. Default `3L`.
#' @param conc_units Character. Concentration units for SSD-eligible rows when
#'   `df` lacks a `units.analyte` column. Passed to [add_mspaf()].
#' @param require_temperature Logical (default `TRUE`). When `TRUE`, any daily
#'   sample with `NH3-N` must also carry a `temperature` value. Passed to
#'   [add_mspaf()]. Set `FALSE` for datasets without ammonia.
#' @param ndraws Positive integer or `NULL` (default). `NULL` returns the
#'   **deterministic** daily msPAF: the fast, grabs-exact point estimate (the
#'   default and recommended best guess).  When non-`NULL`, runs the full
#'   OU/GAM uncertainty propagation for `ndraws` draws and returns a **draws**
#'   product instead (see Details for the distinction).  Requires
#'   `interpolation = "model"`.
#' @param seed Integer or `NULL`. RNG seed for reproducibility of draws.
#' @param return `"summary"` (default) or `"draws"`; relevant only when
#'   `ndraws` is supplied.  `"summary"` collapses the draws to per-envelope
#'   central estimates (`median_*`, the draws' own central tendency — see
#'   `central`) plus credible-interval bounds (`lo_*`, `hi_*`) for the chosen
#'   `gap_uncertainty`, alongside the `deterministic` centre line.  `"draws"`
#'   returns one row per (site \eqn{\times} day \eqn{\times} draw) with a
#'   `draw_id` column and the per-draw envelope value(s).  See *Value*.  Ignored
#'   in point mode (`ndraws = NULL`).
#' @param interval Credible interval width for `return = "summary"` (default
#'   `0.9`).  The lower bound is the `(1 - interval)/2` quantile, the upper
#'   is the `1 - (1 - interval)/2` quantile.
#' @param central Central tendency for the `median_*` envelope columns when
#'   `return = "summary"`: `"median"` (default) or `"mean"` of the per-day
#'   draws.  Because it is a summary *of the draws*, each `median_*` is coherent
#'   with its band (it lies inside `[lo_*, hi_*]` by construction) and depends on
#'   `ndraws`/`seed`.  It is NOT the deterministic point estimate, which is
#'   reported separately in the `deterministic` column (and is the sole estimate
#'   when called with `ndraws = NULL`).
#' @param grab_cv Numeric scalar or named numeric vector of coefficients of
#'   variation for grab-sample measurement error.  A scalar applies the same
#'   CV to all analytes; a named vector (e.g. `c(Cu = 0.1, pH = 0.02)`)
#'   applies analyte-specific CVs.  Controls two uncertainty sources:
#'   S6 (anchor S-value spread at measured dates) and S7 (co-analyte
#'   normalisation spread).  `NULL` (default) disables S6 and S7.
#' @param ou_scale Positive numeric scale factor for the OU bridge envelope
#'   (default `1`).  Multiplies \eqn{\sigma^2} and \eqn{\gamma} (marginal
#'   variance) without changing \eqn{\theta} (correlation length).  Use values
#'   > 1 to widen the between-grab uncertainty bands.
#' @param kappa Non-negative numeric (default `0.5`).  Hydrological modulation
#'   exponent: the smoother's process variance is multiplied by
#'   \eqn{\exp(\kappa \cdot z_h)} where \eqn{z_h} is the standardised daily
#'   flow (capped at ±4 SD).  `kappa = 0` disables hydrological modulation;
#'   larger values widen the credible band more aggressively across high-flow
#'   gaps.  Requires a hydrology column in `df`.
#' @param parallel Logical (default `FALSE`).  When `TRUE`, the per-draw loop
#'   for each site is parallelised via [future.apply::future_lapply()], which
#'   honours whatever [future::plan()] the caller has established.  Requires
#'   the **future.apply** package.  Set a parallel plan before calling, e.g.
#'   `future::plan(future::multisession, workers = 4)`.  In parallel mode the
#'   RNG stream for each draw is managed by `future.apply` (L'Ecuyer-CMRG),
#'   so draws will differ from sequential mode even with the same `seed`, but
#'   are themselves reproducible.
#' @param couple_residuals Logical (default `TRUE`).  When `TRUE` and >= 2
#'   analytes have fitted residual smoothers, daily residual draws are
#'   correlated across analytes using the empirical anchor-residual correlation
#'   (see [.anchor_residual_cor()]).  This widens the combined msPAF interval
#'   to reflect co-movement of co-toxicants on breach events while leaving
#'   per-analyte marginals unchanged.  Set to `FALSE` to reproduce the pre-#32
#'   independent-draw path exactly.
#' @param gap_uncertainty One of `"bracket"` (default), `"ignorable"`, or
#'   `"informative"`; relevant only in draws mode.  In observation gaps the
#'   latent residual reverts to its marginal variance, widening the band — the
#'   honest posterior **only under ignorable (MAR) missingness**.  Field
#'   sampling is often **informative (MNAR)**: gaps exist *because* the system
#'   was judged quiescent, so that band over-states gap uncertainty.  The model
#'   cannot identify the mechanism from grabs, so the default **brackets** both:
#'   the *ignorable* (upper) envelope keeps the simulation-smoother draw; the
#'   *informative* (lower) envelope freezes the residual at its posterior mean on
#'   in-gap days (the fully-informative extreme — gaps perfectly predictable).
#'   The two are nested and coincide at observation days.  `"ignorable"` /
#'   `"informative"` return only that envelope.  See *Interpreting gap
#'   uncertainty* below.  Reference: Rubin (1976) \doi{10.1093/biomet/63.3.581}.
#' @param transform `"pseudo_log"` (default) or `"additive"`.  Controls the
#'   variance-stabilising transform for the daily impact residual smoother.
#'   `"pseudo_log"` applies `g = asinh(I / c)` with per-analyte scale `c = HC5`
#'   (issue #15); `"additive"` keeps `g = I` (pre-#15 behaviour, no HC5 scaling).
#'   Forwarded to [fit_target_model()]; ignored unless
#'   `interpolation = "model"`.
#'
#' @details
#' `mspaf_daily()` returns one of **two distinct products**, chosen by
#' `ndraws`; they answer different questions and should not be expected to
#' coincide:
#'
#' * **Deterministic** (`ndraws = NULL`, default).  A single grabs-exact daily
#'   line: the residual smoother is pinned to the measured grab samples, so the
#'   curve threads the reported lab values.  *Pro:* fast, reproducible, the best
#'   single-number estimate; ideal as a sanity check.  *Con:* carries no
#'   uncertainty and over-fits noisy grabs.
#' * **Draws** (`ndraws > 0`).  A Monte-Carlo posterior propagating trend
#'   (GAM) and between-grab (OU bridge) uncertainty, plus grab measurement error
#'   when `grab_cv` is set.  `return = "summary"` reports the central tendency
#'   plus a credible band; `return = "draws"` returns per-draw paths.  *Pro:*
#'   honest uncertainty quantification.  *Con:* the central estimate is
#'   seed/`ndraws`-dependent and, by Jensen's inequality on the bounded msPAF
#'   index, generally differs from the deterministic line.
#'
#' The summary centre is a summary *of the draws* (issue #42), so it always lies
#' within its own credible band; it is not the deterministic line overlaid on a
#' band built from a different posterior.
#'
#' @section Interpreting gap uncertainty:
#' In draws mode the width of the band across an observation gap is the honest
#' posterior **only under ignorable (MAR) missingness**.  Where you have
#' external grounds that a gap was quiescent (informative/MNAR sampling), the
#' ignorable band over-states uncertainty there.  The `"bracket"` output gives
#' both extremes: read the **informative** (lower) envelope where you vouch the
#' gap was quiet, the **ignorable** (upper) envelope otherwise.  The
#' `precautionary_lo`/`precautionary_hi` columns are the composite
#' `[lo_informative, hi_ignorable]` — a **decision bound, not a calibrated
#' credible interval** (its coverage exceeds nominal and is undefined).  Applied
#' blanket, the informative envelope under-covers genuinely eventful gaps, so use
#' it per-gap; automatic per-gap conditioning on a continuous proxy is future
#' work (#18).  The `deterministic` line is **not** a safe blanket alternative —
#' it under-states risk by Jensen's inequality (#39/#42).
#'
#' @return
#' **Point mode** (`ndraws = NULL`): a tibble with one row per (site
#'   \eqn{\times} day) for days with sufficient analyte coverage; `mspaf` is
#'   the deterministic daily estimate.
#'
#' **Draws mode** (`ndraws > 0`, `return = "summary"`): one row per (site
#'   \eqn{\times} day) with the `deterministic` centre line and the envelope
#'   columns for the chosen `gap_uncertainty`: `median_*`, `lo_*`, `hi_*` per
#'   envelope (`informative` and/or `ignorable`, the central tendency per
#'   `central` and the `interval` credible bounds), plus
#'   `precautionary_lo`/`precautionary_hi` in `"bracket"` mode.
#'
#' **Draws mode** (`ndraws > 0`, `return = "draws"`): one row per (site
#'   \eqn{\times} day \eqn{\times} draw), with a `draw_id` integer column and the
#'   per-draw msPAF value(s) `mspaf_ignorable` and/or `mspaf_informative` for
#'   the chosen `gap_uncertainty`.
#'
#' Common columns:
#'   \describe{
#'     \item{`date`}{Date of this daily estimate.}
#'     \item{`site_id`}{Site identifier.}
#'     \item{`n_analytes_used`}{SSD-eligible analytes contributing to msPAF.}
#'     \item{`dominant_analyte`}{Analyte with the highest individual PAF.}
#'     \item{`max_paf`}{PAF of the dominant analyte (proportion 0--1).}
#'     \item{`n_measured_analytes`}{SSD-eligible analytes with a direct grab
#'       sample on this day (not interpolated).}
#'     \item{`days_since_last_sample`}{Days since the most recent grab sample
#'       for any SSD-eligible analyte.}
#'   }
#'   `"analyte_pafs"` (per-analyte PAF breakdown, re-keyed by date) and
#'   `"ara_summary"` attributes are attached; retrieve them with [analyte_pafs()]
#'   and [ara_summary()]. (`analyte_pafs` is now a flat attribute, not a
#'   list-column — issue #30.)
#'
#' @seealso [add_mspaf()], [time_weighted_aggregate()],
#'   [estimate_water_temp()], [get_silo_air_temp()]
#'
#' @examples
#' \donttest{
#' demo <- leachate_demo()
#' ds <- subset(demo, site_id == "downstream")
#' out <- mspaf_daily(ds, require_temperature = FALSE)
#' head(out[, c(
#'   "date", "site_id", "mspaf", "n_measured_analytes",
#'   "days_since_last_sample"
#' )])
#' }
#' @export
mspaf_daily <- function(
  df,
  temperature = NULL,
  reference = NULL,
  reference_model = NULL,
  imputation_model = NULL,
  start = NULL,
  end = NULL,
  by = "day",
  interpolation = c("forward_fill", "linear", "model"),
  leading_edge = c("drop", "backfill"),
  analyte_metadata = NULL,
  method = c("multi", "anzecc"),
  guideline_dir = getOption("hydroSense.guideline_dir"),
  min_analytes = 3L,
  conc_units = NULL,
  require_temperature = TRUE,
  ndraws = NULL,
  seed = NULL,
  return = c("summary", "draws"),
  interval = 0.9,
  central = c("median", "mean"),
  grab_cv = NULL,
  ou_scale = 1,
  kappa = 0.5,
  parallel = FALSE,
  couple_residuals = TRUE,
  gap_uncertainty = c("bracket", "ignorable", "informative"),
  transform = c("pseudo_log", "additive")
) {
  ## --- Validate inputs -------------------------------------------------------
  checkmate::assert_data_frame(df)
  checkmate::assert_names(names(df),
    must.include = c("sample_id", "site_id", "datetime", "analyte", "value")
  )
  interpolation <- match.arg(interpolation)
  leading_edge <- match.arg(leading_edge)
  method <- match.arg(method)
  return <- match.arg(return)
  central <- match.arg(central)
  checkmate::assert_flag(require_temperature)
  checkmate::assert_int(min_analytes, lower = 1L)
  checkmate::assert_string(by, min.chars = 1L)
  checkmate::assert_number(interval, lower = 0, upper = 1)
  checkmate::assert_number(ou_scale, lower = 0, finite = TRUE)
  checkmate::assert_number(kappa, lower = 0, finite = TRUE)
  if (!is.null(ndraws)) checkmate::assert_count(ndraws, positive = TRUE)
  if (!is.null(grab_cv)) {
    checkmate::assert_numeric(grab_cv, lower = 0, finite = TRUE, min.len = 1L)
  }
  checkmate::assert_flag(parallel)
  checkmate::assert_flag(couple_residuals)
  gap_uncertainty <- match.arg(gap_uncertainty)
  transform <- match.arg(transform)
  if (parallel && !requireNamespace("future.apply", quietly = TRUE)) {
    cli::cli_abort(c(
      "{.arg parallel = TRUE} requires the {.pkg future.apply} package.",
      "i" = "Install it with {.code install.packages(\"future.apply\")} and set a \\
             parallel plan with {.code future::plan(future::multisession)}."
    ))
  }
  if (!is.null(seed) && !parallel) {
    ## Sequential mode: set a global seed and restore on exit.
    ## Parallel mode: future.apply manages per-draw seeds via future.seed.
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    set.seed(as.integer(seed))
    on.exit(
      {
        if (!is.null(old_seed)) {
          assign(".Random.seed", old_seed, envir = .GlobalEnv)
        }
      },
      add = TRUE
    )
  }

  draws_mode <- !is.null(ndraws)

  ## interpolation = "model" needs a fitted reference_model to predict the
  ## site impact; co-analytes are still forward-filled (toxicants come from the
  ## target model). ARA on/off is controlled separately by `reference`.
  if (interpolation == "model" && !inherits(reference_model, "reference_model")) {
    cli::cli_abort(c(
      "{.code interpolation = \"model\"} requires a fitted {.arg reference_model} \\
       (from {.fn fit_reference_model}).",
      "i" = "ARA is controlled separately by {.arg reference}: pass the same \\
             model (or a {.cls prepared_reference}) to subtract background, or \\
             {.val NULL} to assess total concentration."
    ))
  }

  if (draws_mode && interpolation != "model") {
    cli::cli_abort(c(
      "{.arg ndraws} requires {.code interpolation = \"model\"}.",
      "i" = "Draws-mode uncertainty propagation uses the season-blind target \\
             model to generate chemistry traces; supply a {.arg reference_model} \\
             and set {.code interpolation = \"model\"}."
    ))
  }

  ## #50: with rainfall hydrology there is no continuous in-stream state signal,
  ## so the smoother treats observation gaps as ignorable (MAR) and the residual
  ## reverts to its marginal variance there. Warn that the intervals over-state
  ## uncertainty in any gap the analyst knows was quiescent (informative/MNAR),
  ## without steering them to the deterministic line (which under-states risk by
  ## Jensen, #39/#42). Silent in point mode and for stage/discharge hydrology.
  if (draws_mode && inherits(reference_model, "reference_model") &&
    identical(reference_model$hydro_type, "rainfall")) {
    cli::cli_warn(c(
      "Daily gap uncertainty is treated as {.emph ignorable} (MAR) for \\
       rainfall hydrology.",
      "i" = "The residual reverts to its marginal variance across gaps, so the \\
             intervals {.emph over-state} uncertainty in any gap you know was \\
             quiescent. The {.field informative} envelope brackets the opposite \\
             (fully-informative) extreme; read the two together and choose \\
             per gap.",
      "i" = "This is not a reason to prefer the {.field deterministic} line: it \\
             {.emph under-states} risk by Jensen's inequality (#39/#42)."
    ))
  }

  if (!is.null(temperature)) {
    checkmate::assert_data_frame(temperature)
    checkmate::assert_names(names(temperature),
      must.include = c("datetime", "value")
    )
  }

  ## Draw-bearing input chemistry is not supported (mspaf_daily generates its
  ## own draws internally via ndraws; use add_mspaf() directly for existing
  ## draw-carrier frames).
  if ("draw_id" %in% names(df) && !all(is.na(df[["draw_id"]]))) {
    cli::cli_abort(c(
      "{.fn mspaf_daily} does not accept draw-bearing input {.arg df}.",
      "i" = "Use {.arg ndraws} to generate daily draws internally (requires \\
             {.code interpolation = \"model\"}), or call {.fn add_mspaf} on \\
             your draw-carrier frame directly."
    ))
  }

  ## --- Normalise datetime to Date; ensure detected column -------------------
  df <- dplyr::mutate(df, datetime = as.Date(.data$datetime))
  if (!"detected" %in% names(df)) {
    df <- dplyr::mutate(df, detected = TRUE)
  }

  ## --- Determine date range --------------------------------------------------
  if (is.null(start)) start <- min(df$datetime, na.rm = TRUE)
  if (is.null(end)) end <- max(df$datetime, na.rm = TRUE)
  start <- as.Date(start)
  end <- as.Date(end)
  all_dates <- seq(start, end, by = by)

  ## --- Load analyte metadata once (needed for log-space decisions) -----------
  meta <- .load_analyte_metadata(analyte_metadata)

  tox_analytes <- meta$analyte[!is.na(meta$ssd_available) &
    meta$ssd_available == TRUE]

  ## --- Per-site processing ---------------------------------------------------
  sites <- unique(df$site_id)

  ## Reference GAM perturbation (S1) is only meaningful in total-concentration
  ## mode (reference = NULL) where the reference GAM is NOT subtracted from
  ## both sides.  When reference = reference_model, ref cancels in C_excess.
  perturb_ref <- draws_mode && is.null(reference)

  site_results <- lapply(sites, function(site) {
    site_rows <- dplyr::filter(df, .data$site_id == .env$site)

    ## Step 1: Interpolate each analyte onto the daily grid.  For the "model"
    ## path, co-analytes are forward-filled here; toxicants are overwritten
    ## by the fitted target model (step 1b).
    base_interp <- if (interpolation == "model") "forward_fill" else interpolation
    daily_long <- .build_daily_chem(
      site_rows     = site_rows,
      dates         = all_dates,
      interpolation = base_interp,
      leading_edge  = leading_edge,
      tox_analytes  = tox_analytes
    )

    if (nrow(daily_long) == 0L) {
      return(NULL)
    }

    ## Step 2: Fill temperature from external series on non-grab days.
    if (!is.null(temperature)) {
      daily_long <- .fill_external_temperature(daily_long, temperature)
    }

    ## Step 1b: Model interpolation of toxicants (season-blind impact model).
    impact_tiers <- NULL
    ## #50 bracket: parallel synthetic frames for the informative envelope and
    ## the deterministic centre line (NULL unless built in draws mode).
    synth_inf <- NULL
    synth_det <- NULL

    if (interpolation == "model") {
      if (draws_mode) {
        ## ── Draws mode: fit once, iterate N draws ─────────────────────────────
        fdm <- .fit_daily_target(
          site_rows        = site_rows,
          reference_model  = reference_model,
          imputation_model = imputation_model,
          conc_units       = conc_units,
          meta             = meta,
          tox_analytes     = tox_analytes,
          daily_long       = daily_long,
          ou_scale         = ou_scale,
          grab_cv          = grab_cv,
          kappa            = kappa,
          method           = method,
          guideline_dir    = guideline_dir,
          transform        = transform
        )

        ## Diagnostics from forward-filled daily_long (counts grab dates,
        ## independent of draw perturbations).
        diag <- .compute_daily_diag(daily_long, tox_analytes, site)

        if (is.null(fdm)) {
          cli::cli_warn(c(
            "Target model fit failed for site {.val {site}} in draws mode.",
            "i" = "This site will contribute exact forward-fill rows with \\
                   zero-width credible intervals."
          ))
          synth <- .build_synthetic_samples(daily_long, site)
          ## No model -> the bracket collapses: all three frames coincide.
          return(list(
            synth = synth, synth_inf = synth, synth_det = synth,
            diag = diag, tiers = NULL, site = site
          ))
        }

        ## Deterministic prediction: supplies the per-analyte impact tier
        ## ("model"/"bridge") attached to the ARA diagnostics.  Since issue #42
        ## the summary centre line is the draws' own central tendency, so no
        ## separate point-mode msPAF frame is built here (the deterministic
        ## point estimate is its own product, via point mode / ndraws = NULL).
        pt_rows <- .predict_daily_tox(fdm)
        impact_tiers <- if (!is.null(pt_rows)) attr(pt_rows, "impact_tiers") else NULL

        ## Non-modelled + co-analyte rows (exact); shared by every draw's
        ## reconstruction and the co-analyte perturbation path below.
        daily_long_exact <- daily_long[!daily_long$analyte %in% fdm$modelled, ,
          drop = FALSE
        ]
        if (!"units.analyte" %in% names(daily_long_exact)) {
          daily_long_exact$units.analyte <- NA_character_
        }

        ## Compute the cross-analyte empirical correlation once per site, from
        ## the anchor residuals of the fitted target model.  Used by the coupled
        ## draw path below.  Computed unconditionally here so that the seed block
        ## that follows is the SINGLE chokepoint for all RNG state in draw mode.
        cor_res <- .anchor_residual_cor(fdm$tm, fdm$modelled)

        ## Precompute coherent residual trajectories once per analyte — the
        ## SINGLE innovation chokepoint.  S4 (mid-gap residual spread) and S6
        ## (grab measurement error, via the draw model's observation noise) both
        ## enter through these paths; theta/gamma are estimated once (the per-
        ## draw GAM perturbation supplies the S1-S3 trend uncertainty).
        ##
        ## Cross-analyte coupling (#32): when couple_residuals = TRUE and >= 2
        ## analytes have fitted draw models, innovations are correlated via the
        ## Cholesky of cor_res$R.  The DK identity guarantees each per-analyte
        ## marginal is unchanged; only the joint distribution (combined msPAF
        ## interval) widens to reflect co-movement of co-toxicants.
        if (!is.null(seed)) set.seed(as.integer(seed))
        n_couplable <- sum(vapply(
          fdm$smoothers[fdm$modelled],
          function(s) !is.null(s$draw_model), logical(1L)
        ))
        res_draws <- if (couple_residuals && n_couplable >= 2L) {
          .coupled_residual_draws(
            fdm$smoothers, fdm$modelled, ndraws,
            cor_res$R, cor_res$analytes
          )
        } else {
          stats::setNames(lapply(fdm$modelled, function(nm) {
            sm <- fdm$smoothers[[nm]]
            if (is.null(sm) || length(sm$grid_dates) == 0L) {
              return(list(grid_dates = as.Date(character()), draws = NULL))
            }
            dr <- if (is.null(sm$draw_model)) {
              matrix(sm$mean, nrow = length(sm$grid_dates), ncol = ndraws)
            } else {
              .kalman_draw(sm$draw_model, ndraws)
            }
            list(grid_dates = sm$grid_dates, draws = dr)
          }), fdm$modelled)
        }

        ## #50 informative envelope: per analyte, the gap mask (where the
        ## residual variance has ballooned) and the posterior mean, aligned to
        ## each draw's own grid. On in-gap days the informative path replaces the
        ## simulation-smoother draw by this mean (freezing the residual at its
        ## posterior expectation -> the fully-informative extreme); off-gap the
        ## draw is retained, so observed-day uncertainty (incl. S6 grab error) is
        ## preserved. Built once here so the per-draw closure stays RNG-neutral.
        want_inf <- gap_uncertainty %in% c("bracket", "informative")
        res_frozen <- if (!want_inf) {
          NULL
        } else {
          stats::setNames(lapply(fdm$modelled, function(nm) {
            sm <- fdm$smoothers[[nm]]
            rd <- res_draws[[nm]]
            if (is.null(rd$draws) || is.null(sm) ||
              length(sm$grid_dates) == 0L) {
              return(NULL)
            }
            mask <- .residual_gap_mask(
              list(grid_dates = sm$grid_dates, var = sm$var),
              sm$anchor_dates
            )
            ## Align mask + mean to the draw grid (defensive: identical here).
            mask_lut <- stats::setNames(mask, as.character(sm$grid_dates))
            mean_lut <- stats::setNames(sm$mean, as.character(sm$grid_dates))
            gd <- as.character(rd$grid_dates)
            m_rd <- mask_lut[gd]
            m_rd[is.na(m_rd)] <- FALSE
            list(mask = as.logical(m_rd), mean = as.numeric(mean_lut[gd]))
          }), fdm$modelled)
        }

        ## N stochastic draw iterations (draw_id = 1..N).
        ## G2: .predict_daily_tox uses fdm$co_split (exact) for clean C_raw
        ## reconstruction; S7 co-analyte perturbations enter add_mspaf's
        ## normalisation via draw-bearing co-analyte rows in the synthetic frame.
        ##
        ## H2: the draw loop is a closure so it can run sequentially (lapply) or
        ## in parallel (future.apply::future_lapply).
        .draw_fn <- local({
          .fdm <- fdm
          .perturb_ref <- perturb_ref
          .grab_cv <- grab_cv
          .daily_long_exact <- daily_long_exact
          .res_draws <- res_draws
          .res_frozen <- res_frozen
          .want_inf <- want_inf
          function(d_idx) {
            tm_p <- .perturb_target_model(.fdm$tm, perturb_reference = .perturb_ref)
            ## Residual path for this draw (S impact / d WQ), per analyte.
            residual_paths <- stats::setNames(lapply(.fdm$modelled, function(nm) {
              rd <- .res_draws[[nm]]
              if (is.null(rd$draws)) {
                return(rep(NA_real_, length(.fdm$qdates)))
              }
              .residual_on_qdates(rd$grid_dates, rd$draws[, d_idx], .fdm$qdates)
            }), .fdm$modelled)
            ## S7: perturbed wq_long shifts WQ-layer scores; co_split stays exact.
            co_p_wq <- if (!is.null(.grab_cv)) {
              .perturb_co_split(.fdm, .grab_cv)
            } else {
              list(co_split = .fdm$co_split, wq_long = .fdm$wq_long)
            }
            mr_d <- .predict_daily_tox(
              .fdm,
              tm_p           = tm_p,
              residual_paths = residual_paths,
              co_split       = .fdm$co_split, # exact → clean C_raw
              wq_long        = co_p_wq$wq_long # perturbed when S7 active
            )
            if (!is.null(mr_d)) mr_d$draw_id <- as.integer(d_idx)

            ## Informative envelope: same perturbations (tm_p, wq), residual
            ## frozen at its posterior mean on in-gap days. Deterministic given
            ## the already-drawn quantities -> consumes no RNG (preserves the
            ## ignorable draws exactly).
            mr_d_inf <- NULL
            if (.want_inf) {
              residual_paths_inf <- stats::setNames(lapply(.fdm$modelled, function(nm) {
                rd <- .res_draws[[nm]]
                if (is.null(rd$draws)) {
                  return(rep(NA_real_, length(.fdm$qdates)))
                }
                col <- rd$draws[, d_idx]
                fr <- .res_frozen[[nm]]
                if (!is.null(fr) && any(fr$mask)) {
                  col[fr$mask] <- fr$mean[fr$mask]
                }
                .residual_on_qdates(rd$grid_dates, col, .fdm$qdates)
              }), .fdm$modelled)
              mr_d_inf <- .predict_daily_tox(
                .fdm,
                tm_p           = tm_p,
                residual_paths = residual_paths_inf,
                co_split       = .fdm$co_split,
                wq_long        = co_p_wq$wq_long
              )
              if (!is.null(mr_d_inf)) mr_d_inf$draw_id <- as.integer(d_idx)
            }

            list(
              tox_rows = mr_d,
              tox_rows_inf = mr_d_inf,
              co_rows = if (!is.null(.grab_cv)) {
                .co_draw_rows(.daily_long_exact, co_p_wq$co_split, d_idx)
              } else {
                NULL
              }
            )
          }
        })

        draw_results <- if (parallel) {
          future.apply::future_lapply(
            seq_len(ndraws), .draw_fn,
            future.seed = if (!is.null(seed)) as.integer(seed) else TRUE
          )
        } else {
          lapply(seq_len(ndraws), .draw_fn)
        }

        tox_draw_rows <- lapply(draw_results, `[[`, "tox_rows")
        co_draw_rows <- if (!is.null(grab_cv)) {
          lapply(draw_results, `[[`, "co_rows")
        } else {
          NULL
        }

        tox_draw_long <- dplyr::bind_rows(Filter(Negate(is.null), tox_draw_rows))

        ## Assemble a draws-mode daily_long from a set of modelled-tox draw rows,
        ## reusing the (shared) exact co-analyte / non-modelled context. The
        ## ignorable and informative envelopes differ ONLY in their modelled-tox
        ## rows, so they share this assembler exactly.
        non_mod_tox <- daily_long_exact[
          daily_long_exact$analyte %in% tox_analytes, ,
          drop = FALSE
        ]
        co_draws_all <- if (!is.null(co_draw_rows)) {
          dplyr::bind_rows(Filter(Negate(is.null), co_draw_rows))
        } else {
          NULL
        }
        assemble_daily <- function(tox_long) {
          if (!is.null(co_draws_all)) {
            ## S7 active: co-analytes are draw-bearing; non-modelled tox exact.
            dplyr::bind_rows(non_mod_tox, co_draws_all, tox_long)
          } else {
            ## S7 inactive: co-analytes + non-modelled tox exact (draw_id = NA →
            ## broadcast to all stochastic draws by add_mspaf).
            dplyr::bind_rows(daily_long_exact, tox_long)
          }
        }
        daily_long <- assemble_daily(tox_draw_long)

        ## #50 informative envelope frame (same draws, residual frozen in gaps).
        if (want_inf) {
          tox_inf_rows <- lapply(draw_results, `[[`, "tox_rows_inf")
          tox_inf_long <- dplyr::bind_rows(Filter(Negate(is.null), tox_inf_rows))
          synth_inf <- .build_synthetic_samples(
            assemble_daily(tox_inf_long), site
          )
        }

        ## #50 deterministic centre frame: posterior-mean residual, unperturbed
        ## reference — identical to point mode (.daily_tox_from_model reuses the
        ## same fit + .predict_daily_tox(fdm)), so the column matches ndraws=NULL.
        if (!is.null(pt_rows)) {
          daily_det <- dplyr::bind_rows(daily_long_exact, pt_rows)
          synth_det <- .build_synthetic_samples(daily_det, site)
        }
      } else {
        ## ── Point mode: thin wrapper ──────────────────────────────────────────
        daily_long <- .daily_tox_from_model(
          daily_long       = daily_long,
          site_rows        = site_rows,
          reference_model  = reference_model,
          imputation_model = imputation_model,
          conc_units       = conc_units,
          meta             = meta,
          tox_analytes     = tox_analytes,
          method           = method,
          guideline_dir    = guideline_dir,
          transform        = transform
        )
        impact_tiers <- attr(daily_long, "impact_tiers")
        diag <- .compute_daily_diag(daily_long, tox_analytes, site)
      }
    } else {
      diag <- .compute_daily_diag(daily_long, tox_analytes, site)
    }

    ## Step 4: Build synthetic long-format samples (no focal_date!).
    synth <- .build_synthetic_samples(daily_long, site)

    list(
      synth = synth, synth_inf = synth_inf, synth_det = synth_det,
      diag = diag, tiers = impact_tiers, site = site
    )
  })

  site_results <- Filter(Negate(is.null), site_results)
  if (length(site_results) == 0L) {
    cli::cli_warn("No daily chemistry could be built. Returning empty tibble.")
    empty_mode <- if (!draws_mode) "point" else if (return == "draws") "draws" else "summary"
    return(.empty_daily_result(empty_mode, gap_uncertainty))
  }

  all_synth <- dplyr::bind_rows(lapply(site_results, `[[`, "synth"))
  all_diag <- dplyr::bind_rows(lapply(site_results, `[[`, "diag"))

  ## Build sample_id -> date lookup before passing to add_mspaf() (add_mspaf
  ## may rearrange rows but sample_id is stable throughout).
  id_date_map <- dplyr::distinct(
    dplyr::select(all_synth, "sample_id", "site_id", ".date")
  )

  ## Remove .date from the df passed to add_mspaf() so it stays unaware of it.
  all_synth_clean <- dplyr::select(all_synth, -".date")

  ## --- Run add_mspaf on the daily synthetic samples -------------------------
  ## In draws mode, request raw per-draw rows so we can extract draw_id and
  ## per-analyte diagnostics before collapsing ourselves.
  mspaf_out <- add_mspaf(
    df                  = all_synth_clean,
    reference           = reference,
    analyte_metadata    = analyte_metadata,
    method              = method,
    guideline_dir       = guideline_dir,
    min_analytes        = min_analytes,
    conc_units          = conc_units,
    require_temperature = require_temperature,
    return              = if (draws_mode) "draws" else "summary"
  )

  ara_summ <- attr(mspaf_out, "ara_summary")

  ## Attach the target model's per-analyte impact tier ("model" / "bridge") to
  ## the ARA diagnostics.
  all_tiers <- dplyr::bind_rows(lapply(site_results, function(z) {
    if (is.null(z$tiers) || nrow(z$tiers) == 0L) {
      return(NULL)
    }
    dplyr::mutate(z$tiers, site_id = z$site)
  }))
  if (!is.null(ara_summ) && nrow(all_tiers) > 0L) {
    site_lookup <- dplyr::distinct(id_date_map[, c("sample_id", "site_id")])
    ara_summ <- ara_summ |>
      dplyr::left_join(site_lookup, by = "sample_id") |>
      dplyr::left_join(all_tiers, by = c("site_id", "analyte"))
  }

  ## --- Extract and annotate msPAF rows -------------------------------------
  mspaf_rows <- dplyr::filter(mspaf_out, .data$analyte == "msPAF")

  if (nrow(mspaf_rows) == 0L) {
    cli::cli_warn(
      "No daily msPAF rows produced. \\
       Check {.arg min_analytes} ({min_analytes}) and data coverage."
    )
    empty_mode <- if (!draws_mode) "point" else if (return == "draws") "draws" else "summary"
    result <- .empty_daily_result(empty_mode, gap_uncertainty)
    attr(result, "ara_summary") <- ara_summ
    return(result)
  }

  mspaf_dated <- mspaf_rows |>
    dplyr::left_join(id_date_map, by = c("sample_id", "site_id")) |>
    dplyr::rename(date = ".date")

  if (draws_mode) {
    ## #50 bracket. The primary pass (mspaf_dated) is the IGNORABLE envelope —
    ## byte-identical to the pre-#50 behaviour. The informative envelope and the
    ## deterministic centre are additive secondary add_mspaf passes over the
    ## frozen-residual / posterior-mean synthetic frames (all share the same
    ## seeded draws, so they are deterministic given the primary pass).
    want_inf <- gap_uncertainty %in% c("bracket", "informative")

    run_pass <- function(col, ret) {
      s <- dplyr::bind_rows(lapply(site_results, `[[`, col))
      if (is.null(s) || nrow(s) == 0L) {
        return(NULL)
      }
      out <- add_mspaf(
        df                  = dplyr::select(s, -".date"),
        reference           = reference,
        analyte_metadata    = analyte_metadata,
        method              = method,
        guideline_dir       = guideline_dir,
        min_analytes        = min_analytes,
        conc_units          = conc_units,
        require_temperature = require_temperature,
        return              = ret
      )
      rows <- dplyr::filter(out, .data$analyte == "msPAF")
      if (nrow(rows) == 0L) {
        return(NULL)
      }
      rows |>
        dplyr::left_join(id_date_map, by = c("sample_id", "site_id")) |>
        dplyr::rename(date = ".date")
    }

    ig_draws <- mspaf_dated |>
      dplyr::select(
        "date", "site_id", "draw_id",
        mspaf_ignorable = "value",
        "n_analytes_used", "dominant_analyte", "max_paf"
      )
    inf_dated <- if (want_inf) run_pass("synth_inf", "draws") else NULL
    draws_wide <- ig_draws
    if (!is.null(inf_dated)) {
      draws_wide <- dplyr::full_join(
        draws_wide,
        dplyr::select(inf_dated, "date", "site_id", "draw_id",
          mspaf_informative = "value"
        ),
        by = c("date", "site_id", "draw_id")
      )
    }
    ## Informative requested but unproducible -> collapse onto ignorable.
    if (want_inf && !"mspaf_informative" %in% names(draws_wide)) {
      draws_wide$mspaf_informative <- draws_wide$mspaf_ignorable
    }

    if (return == "draws") {
      keep <- .bracket_draw_cols(gap_uncertainty)
      result <- draws_wide |>
        dplyr::left_join(all_diag, by = c("date", "site_id")) |>
        dplyr::select(
          "date", "site_id", "draw_id", dplyr::all_of(keep),
          "n_analytes_used", "dominant_analyte", "max_paf",
          "n_measured_analytes", "days_since_last_sample"
        ) |>
        dplyr::arrange(.data$site_id, .data$date, .data$draw_id)
    } else {
      ## Envelope summary (median + CI per #42), the deterministic centre line,
      ## and the composition diagnostics (envelope-invariant; from one draw).
      env_summ <- .summarise_bracket(
        draws_wide,
        interval = interval, central = central,
        gap_uncertainty = gap_uncertainty
      )
      det_dated <- run_pass("synth_det", "summary")
      det_col <- if (!is.null(det_dated)) {
        dplyr::select(det_dated, "date", "site_id", deterministic = "value")
      } else {
        tibble::tibble(
          date = env_summ$date, site_id = env_summ$site_id,
          deterministic = NA_real_
        )
      }
      comp <- mspaf_dated |>
        dplyr::group_by(.data$date, .data$site_id) |>
        dplyr::summarise(
          n_analytes_used = dplyr::first(.data$n_analytes_used),
          dominant_analyte = dplyr::first(.data$dominant_analyte),
          max_paf = dplyr::first(.data$max_paf),
          .groups = "drop"
        )
      result <- env_summ |>
        dplyr::left_join(det_col, by = c("date", "site_id")) |>
        dplyr::left_join(comp, by = c("date", "site_id")) |>
        dplyr::left_join(all_diag, by = c("date", "site_id")) |>
        dplyr::select(
          "date", "site_id", "deterministic",
          dplyr::all_of(.bracket_summary_cols(gap_uncertainty)),
          "n_analytes_used", "dominant_analyte", "max_paf",
          "n_measured_analytes", "days_since_last_sample"
        ) |>
        dplyr::arrange(.data$site_id, .data$date)
    }
  } else {
    ## Point mode: unchanged output schema.
    result <- mspaf_dated |>
      dplyr::rename(mspaf = "value") |>
      dplyr::left_join(all_diag, by = c("date", "site_id")) |>
      dplyr::select(
        "date", "site_id", "mspaf",
        "n_analytes_used", "dominant_analyte", "max_paf",
        "n_measured_analytes", "days_since_last_sample"
      ) |>
      dplyr::arrange(.data$site_id, .data$date)
  }

  ## Per-analyte PAF breakdown, re-keyed from the synthetic sample_id to date.
  apafs <- attr(mspaf_out, "analyte_pafs")
  if (!is.null(apafs) && nrow(apafs) > 0L && exists("id_date_map")) {
    apafs <- apafs |>
      dplyr::left_join(
        dplyr::distinct(id_date_map, .data$sample_id, .data$site_id, .data$.date),
        by = intersect(c("sample_id", "site_id"), names(apafs))
      ) |>
      dplyr::rename(date = ".date") |>
      dplyr::select(-dplyr::any_of("sample_id"))
  }
  attr(result, "analyte_pafs") <- apafs
  attr(result, "ara_summary") <- ara_summ
  result
}


## ============================================================================
## Internal helpers
## ============================================================================

#' Interpolate per-analyte grab chemistry onto a daily date grid
#'
#' For each analyte present in `site_rows`, produces one row per date in
#' `dates` using either forward-fill or linear interpolation.  Adds `.date`
#' and `.measured` columns for downstream diagnostics.
#'
#' @param site_rows Long-format chemistry for one site.
#' @param dates Date vector giving the target daily grid.
#' @param interpolation `"forward_fill"` or `"linear"`.
#' @param leading_edge `"drop"` or `"backfill"`.
#' @param tox_analytes SSD-eligible analyte names (for log-space decision).
#' @return Long-format tibble with `.date`, `analyte`, `value`, `detected`,
#'   `.measured`, and any pass-through columns (`units.analyte`, etc.).
#' @keywords internal
.build_daily_chem <- function(site_rows, dates, interpolation, leading_edge,
                              tox_analytes) {
  analytes <- unique(site_rows$analyte)
  passthru_cols <- intersect(
    names(site_rows),
    c(
      "units.analyte", "valence.analyte",
      "atomic_mass.analyte"
    )
  )

  rows_per_analyte <- lapply(analytes, function(a) {
    a_rows <- dplyr::filter(site_rows, .data$analyte == .env$a) |>
      dplyr::arrange(.data$datetime)

    ## Use detected rows as interpolation anchors; fall back to all rows if none.
    obs <- dplyr::filter(a_rows, .data$detected)
    if (nrow(obs) == 0L) obs <- a_rows
    obs <- dplyr::arrange(obs, .data$datetime)

    interp <- .interpolate_analyte(
      obs_dates     = obs$datetime,
      obs_values    = obs$value,
      obs_detected  = obs$detected,
      target_dates  = dates,
      interpolation = interpolation,
      leading_edge  = leading_edge,
      log_space     = a %in% tox_analytes
    )

    if (nrow(interp) == 0L) {
      return(NULL)
    }

    interp$analyte <- a
    ## Propagate ancillary columns from the last observed row.
    for (col in passthru_cols) {
      if (col %in% names(a_rows)) {
        interp[[col]] <- a_rows[[col]][nrow(a_rows)]
      }
    }
    interp
  })

  dplyr::bind_rows(Filter(Negate(is.null), rows_per_analyte))
}


#' Interpolate one analyte onto a target date vector
#'
#' Core per-analyte interpolation. Returns a tibble with `.date`, `value`,
#' `detected`, `.measured`.
#'
#' @keywords internal
.interpolate_analyte <- function(obs_dates, obs_values, obs_detected,
                                 target_dates, interpolation, leading_edge,
                                 log_space = FALSE) {
  if (length(obs_dates) == 0L) {
    return(tibble::tibble(
      .date     = as.Date(character()),
      value     = numeric(),
      detected  = logical(),
      .measured = logical()
    ))
  }

  ## De-duplicate by date: take the first occurrence per date.
  uniq_mask <- !duplicated(obs_dates)
  obs_d <- obs_dates[uniq_mask]
  obs_v <- obs_values[uniq_mask]
  obs_det <- obs_detected[uniq_mask]

  ## Sort ascending.
  ord <- order(obs_d)
  obs_d <- obs_d[ord]
  obs_v <- obs_v[ord]
  obs_det <- obs_det[ord]

  ## Work in log-space for toxicants (floor at machine epsilon to avoid log(0)).
  if (log_space) {
    obs_v_work <- log(pmax(obs_v, .Machine$double.eps))
  } else {
    obs_v_work <- obs_v
  }

  n <- length(target_dates)
  out_val <- numeric(n)
  out_det <- logical(n)
  out_meas <- logical(n)
  out_keep <- logical(n)

  ## For each target date find its bracketing observations.
  ## findInterval returns 0 when d < obs_d[1], else the index of the largest
  ## obs_d <= d.
  prev_idx <- findInterval(as.numeric(target_dates), as.numeric(obs_d))

  for (i in seq_len(n)) {
    d <- target_dates[i]
    pid <- prev_idx[i]

    exact_match <- pid > 0L && obs_d[pid] == d

    if (exact_match) {
      raw <- obs_v_work[pid]
      out_val[i] <- if (log_space) exp(raw) else raw
      out_det[i] <- obs_det[pid]
      out_meas[i] <- TRUE
      out_keep[i] <- TRUE
      next
    }

    if (pid == 0L) {
      ## Before the first observation.
      if (leading_edge == "backfill") {
        raw <- obs_v_work[1L]
        out_val[i] <- if (log_space) exp(raw) else raw
        out_det[i] <- obs_det[1L]
        out_meas[i] <- FALSE
        out_keep[i] <- TRUE
      }
      ## else: drop (out_keep stays FALSE)
      next
    }

    ## We have a previous observation but no exact match.
    nid <- pid + 1L ## index of the next observation (may be out-of-bounds)

    if (interpolation == "forward_fill" || nid > length(obs_d)) {
      ## Forward-fill, or past the last observation: carry forward.
      raw <- obs_v_work[pid]
      out_val[i] <- if (log_space) exp(raw) else raw
      out_det[i] <- obs_det[pid]
      out_meas[i] <- FALSE
      out_keep[i] <- TRUE
    } else {
      ## Linear / log-linear interpolation between pid and nid.
      frac <- as.numeric(d - obs_d[pid]) /
        as.numeric(obs_d[nid] - obs_d[pid])
      raw <- obs_v_work[pid] + frac * (obs_v_work[nid] - obs_v_work[pid])
      out_val[i] <- if (log_space) exp(raw) else raw
      ## A gap between detected/BDL passes detected only if both anchors are.
      out_det[i] <- obs_det[pid] && obs_det[nid]
      out_meas[i] <- FALSE
      out_keep[i] <- TRUE
    }
  }

  tibble::tibble(
    .date     = target_dates[out_keep],
    value     = out_val[out_keep],
    detected  = out_det[out_keep],
    .measured = out_meas[out_keep]
  )
}


#' Fill temperature from an external daily series on non-grab days
#'
#' Augments `daily_long` with temperature values from `temperature_df` for
#' dates where there is no directly measured grab temperature.  Dates with
#' a `.measured == TRUE` temperature row (actual grab measurement) are kept
#' unchanged.
#'
#' @param daily_long Output of `.build_daily_chem()`.
#' @param temperature_df Data frame with `datetime` and `value` columns.
#' @return `daily_long` with temperature rows augmented.
#' @keywords internal
.fill_external_temperature <- function(daily_long, temperature_df) {
  ext <- tibble::tibble(
    .date     = as.Date(temperature_df$datetime),
    value     = temperature_df$value,
    detected  = TRUE,
    .measured = FALSE,
    analyte   = "temperature"
  )

  ## Dates that already have a directly measured grab temperature.
  grab_temp_dates <- daily_long$.date[
    daily_long$analyte == "temperature" & daily_long$.measured
  ]

  ## Remove existing interpolated temperature rows that the external series
  ## will replace (do not disturb grab-measured rows).
  no_existing_grab <- !(ext$.date %in% grab_temp_dates)
  new_rows <- ext[no_existing_grab, , drop = FALSE]

  ## Drop any interpolated temp rows for those dates before adding new ones.
  daily_long <- dplyr::filter(
    daily_long,
    !(.data$analyte == "temperature" &
      .data$.date %in% new_rows$.date &
      !.data$.measured)
  )

  dplyr::bind_rows(daily_long, new_rows)
}


#' Compute per-day diagnostics: n_measured_analytes, days_since_last_sample
#'
#' Operates on SSD-eligible rows only (toxicants drive the msPAF; co-analyte
#' sampling frequency is generally higher and not the bottleneck).
#'
#' @param daily_long Output of `.build_daily_chem()` (possibly augmented).
#' @param tox_analytes SSD-eligible analyte names.
#' @param site Site identifier for the output `site_id` column.
#' @return Tibble `(date, site_id, n_measured_analytes, days_since_last_sample)`.
#' @keywords internal
.compute_daily_diag <- function(daily_long, tox_analytes, site) {
  tox_rows <- dplyr::filter(daily_long, .data$analyte %in% .env$tox_analytes)

  daily_dates <- sort(unique(tox_rows$.date))

  if (length(daily_dates) == 0L) {
    return(tibble::tibble(
      date = as.Date(character()),
      site_id = character(),
      n_measured_analytes = integer(),
      days_since_last_sample = integer()
    ))
  }

  ## All dates on which any SSD-eligible analyte had a grab measurement.
  grab_dates <- sort(unique(tox_rows$.date[tox_rows$.measured]))

  purrr::map_dfr(daily_dates, function(d) {
    n_meas <- sum(tox_rows$.date == d & tox_rows$.measured, na.rm = TRUE)
    prev <- grab_dates[grab_dates <= d]
    days_since <- if (length(prev) > 0L) {
      as.integer(d - max(prev))
    } else {
      NA_integer_
    }
    tibble::tibble(
      date = d,
      site_id = site,
      n_measured_analytes = as.integer(n_meas),
      days_since_last_sample = days_since
    )
  })
}


#' Build synthetic long-format daily samples from interpolated chemistry
#'
#' Assigns `sample_id = "daily_{YYYY-MM-DD}_{site}"` per day.  Keeps `.date`
#' as a column (the caller extracts it before passing to [add_mspaf()]).
#' No `focal_date` column is added -- this is deliberate so [add_mspaf()]
#' uses the instant (pointwise) ARA path, not the chronic integrated path.
#'
#' @param daily_long Output of `.build_daily_chem()` (after temperature fill).
#' @param site Site identifier string.
#' @return Long-format tibble ready for [add_mspaf()] (after removing `.date`).
#' @keywords internal
.build_synthetic_samples <- function(daily_long, site) {
  if (nrow(daily_long) == 0L) {
    return(tibble::tibble(
      .date = as.Date(character()), sample_id = character(),
      site_id = character(), datetime = as.Date(character()),
      analyte = character(), value = numeric(), detected = logical()
    ))
  }

  daily_long |>
    dplyr::mutate(
      sample_id = paste0(
        "daily_", format(.data$.date, "%Y-%m-%d"),
        "_", .env$site
      ),
      site_id = .env$site,
      datetime = .data$.date
    ) |>
    dplyr::select(-".measured")
}

#' Fit-once scaffold for the season-blind daily target model (issue #16)
#'
#' Calls [fit_target_model()] once, builds the static scaffolding needed for
#' per-draw prediction, and precomputes the OU bridge Cholesky factors for
#' every analyte over the full daily date grid.  The result ("fitted daily
#' model", `fdm`) is the input to [.predict_daily_tox()]; call
#' `.predict_daily_tox(fdm)` once for point mode or N times (with a perturbed
#' target model and OU-sampled `eps_paths`) for draw mode.
#'
#' @param site_rows Grab chemistry for the site.
#' @param reference_model,imputation_model,conc_units,meta,tox_analytes
#'   Forwarded to [fit_target_model()].
#' @param daily_long Forward-filled daily chemistry for the site; supplies the
#'   date grid (`qdates`) and co-analyte values.
#' @return A named list (`fdm`) or `NULL` on model failure.
#' @keywords internal
.fit_daily_target <- function(site_rows, reference_model, imputation_model,
                              conc_units, meta, tox_analytes, daily_long,
                              ou_scale = 1, grab_cv = NULL, kappa = 0.5,
                              method = "multi", guideline_dir = NULL,
                              transform = "pseudo_log") {
  qdates <- sort(unique(daily_long$.date))

  tm <- tryCatch(
    fit_target_model(
      target           = site_rows,
      reference_model  = reference_model,
      imputation_model = imputation_model,
      conc_units       = conc_units,
      analyte_metadata = meta,
      method           = method,
      guideline_dir    = guideline_dir,
      transform        = transform
    ),
    error = function(e) {
      cli::cli_warn(c(
        "Target model fit failed; falling back to forward-fill for toxicants.",
        "x" = conditionMessage(e)
      ))
      NULL
    }
  )
  if (is.null(tm) || length(tm$models) == 0L) {
    return(NULL)
  }

  modelled <- names(tm$models)

  ## Reference norm over the daily grid — static across draws (ARA cancels it),
  ## so resolve it ONCE here and reuse for the centre line and every draw,
  ## rather than re-evaluating the reference GAM per draw.
  ref_q <- tryCatch(
    .resolve_ref_norm_instant(
      tm$reference_model,
      tibble::tibble(sample_id = as.character(qdates), datetime = qdates)
    ) |>
      dplyr::mutate(date = as.Date(.data$sample_id)),
    error = function(e) NULL
  )

  ## Per-date co-analyte lookup (static across draws; Chunk D may override).
  co <- daily_long[!daily_long$analyte %in% tox_analytes &
    (is.na(daily_long$detected) | daily_long$detected), , drop = FALSE]
  co_split <- split(
    data.frame(analyte = co$analyte, value = co$value, stringsAsFactors = FALSE),
    as.character(co$.date)
  )

  ## WQ layer data for analytes with a fitted WQ→metal response.
  wq_long <- if (!is.null(tm$pca)) {
    tibble::tibble(
      sample_id = as.character(co$.date),
      analyte = co$analyte, value = co$value
    )
  } else {
    NULL
  }

  ## Normalisation formula lookup (one parsed formula per modelled analyte).
  meta_norm <- meta |>
    dplyr::select("analyte", "normalisation_formula", "coanalytes_required")
  fac_lookup <- stats::setNames(
    lapply(modelled, function(a) {
      row <- meta_norm[meta_norm$analyte == a, , drop = FALSE]
      list(
        parsed = if (nrow(row)) {
          .parse_normalisation_formula(
            row$normalisation_formula %||% ""
          )
        } else {
          NULL
        }
      )
    }), modelled
  )

  ## Measured (grab) dates — used to set .measured flags on synthetic rows.
  sr_mod <- site_rows[
    site_rows$analyte %in% modelled &
      (is.na(site_rows$detected) | site_rows$detected), ,
    drop = FALSE
  ]
  measured_key <- paste(sr_mod$analyte, as.Date(sr_mod$datetime))

  ## State-space residual smoother: build once per analyte over qdates.
  ## WQ-tier analytes smooth their residual d (d_anchors); others smooth the
  ## impact residual S (anchors). The CENTRE model pins to anchors (r -> 0) and
  ## its posterior mean is the deterministic centre line; the DRAW model carries
  ## the S6 grab measurement error as observation noise r_i. Hydrology modulates
  ## the process variance (q_mult) inside .analyte_residual_smoother().
  get_cv <- function(nm) {
    if (is.null(grab_cv)) {
      return(NA_real_)
    }
    if (length(grab_cv) == 1L) {
      return(as.numeric(grab_cv))
    }
    if (nm %in% names(grab_cv)) {
      return(as.numeric(grab_cv[[nm]]))
    }
    NA_real_
  }

  smoothers <- stats::setNames(
    lapply(modelled, function(nm) {
      m <- tm$models[[nm]]
      has_wq <- !is.null(m$wq_fit) && !is.null(m$d_anchors) &&
        nrow(m$d_anchors) >= 2L
      anch <- if (has_wq) m$d_anchors else m$anchors
      empty <- list(
        grid_dates = as.Date(character()), mean = numeric(0),
        var = numeric(0), anchor_dates = as.Date(character()),
        draw_model = NULL
      )
      if (is.null(anch) || nrow(anch) < 1L) {
        return(empty)
      }

      ## Centre smoother (r -> 0): posterior mean = deterministic centre line.
      sm_c <- .analyte_residual_smoother(m, tm, qdates,
        kappa = kappa,
        scale = ou_scale, r_vec = NULL
      )
      if (is.null(sm_c)) {
        return(empty)
      }

      ## Draw model: add S6 grab measurement error as observation noise r_i,
      ## scaled by C_norm_obs at the anchors (model: C_norm = ref + I;
      ## WQ-tier: C_norm = WQ_pred + d).
      draw_model <- sm_c$model
      cv <- get_cv(nm)
      if (!is.na(cv) && cv > 0 && !is.null(sm_c$model)) {
        c_norm_obs <- tryCatch(
          {
            if (!has_wq) {
              ref_q <- .resolve_ref_norm_instant(
                tm$reference_model,
                tibble::tibble(
                  sample_id = as.character(anch$date),
                  datetime = anch$date
                )
              )
              ref_lkp <- stats::setNames(
                ref_q$ref_norm[ref_q$analyte == nm],
                ref_q$sample_id[ref_q$analyte == nm]
              )
              rv <- as.numeric(ref_lkp[as.character(anch$date)])
              rv[is.na(rv)] <- 0
              pmax(anch$I + rv, 0)
            } else if (!is.null(tm$pca) && !is.null(wq_long)) {
              pc_anch <- .compute_pca_scores(wq_long, tm$pca)
              nd_anch <- dplyr::left_join(
                tibble::tibble(sample_id = as.character(anch$date)),
                pc_anch,
                by = "sample_id"
              )
              pmax(as.numeric(stats::predict(m$wq_fit, newdata = nd_anch)) +
                anch$S, 0)
            } else {
              NULL
            }
          },
          error = function(e) NULL
        )
        if (!is.null(c_norm_obs)) {
          ## S6 is a concentration-space variance; the smoother now works on the
          ## g = asinh(x/c) scale (issue #15), so map it via the delta method.
          ## The transformed level x is the impact I (impact tier) or the
          ## concentration (WQ tier); NA scale_c keeps the additive variance.
          var_c <- (cv * c_norm_obs)^2
          sc_a <- m$scale_c %||% NA_real_
          r_vec <- if (is.na(sc_a)) {
            var_c
          } else {
            level <- if (!has_wq) anch$I else c_norm_obs
            .s6_var_to_g(var_c, level, sc_a)
          }
          sm_d <- .analyte_residual_smoother(m, tm, qdates,
            kappa = kappa,
            scale = ou_scale, r_vec = r_vec
          )
          if (!is.null(sm_d) && !is.null(sm_d$model)) draw_model <- sm_d$model
        }
      }

      ## var + anchor_dates are carried so the #50 bracket can locate in-gap
      ## days (where the informative envelope freezes the residual at its mean).
      list(
        grid_dates = sm_c$grid_dates, mean = sm_c$mean, var = sm_c$var,
        anchor_dates = as.Date(anch$date), draw_model = draw_model
      )
    }),
    modelled
  )

  ## Co-source grab mapping for S7 (coherent per-grab co-analyte perturbation).
  ## co_grab_map[[analyte]][[date_str]] = source grab date string (or NA).
  co_analytes_nm <- unique(co$analyte)
  co_grab_map <- if (length(co_analytes_nm) > 0L) {
    co_grabs <- site_rows[
      site_rows$analyte %in% co_analytes_nm &
        (is.na(site_rows$detected) | site_rows$detected), ,
      drop = FALSE
    ]
    co_grab_by_a <- split(as.Date(co_grabs$datetime), co_grabs$analyte)
    stats::setNames(
      lapply(co_analytes_nm, function(a) {
        gd <- sort(unique(co_grab_by_a[[a]]))
        if (length(gd) == 0L) {
          return(stats::setNames(rep(NA_character_, length(qdates)), as.character(qdates)))
        }
        idx <- findInterval(qdates, gd) # 0 = before first grab
        src <- ifelse(idx == 0L, NA_character_, as.character(gd[pmax(idx, 1L)]))
        stats::setNames(src, as.character(qdates))
      }),
      co_analytes_nm
    )
  } else {
    NULL
  }

  ## Static per-analyte context reused across draws (the per-draw speed path):
  ##  - hydro features + ref vector (constant across draws),
  ##  - the impact-GAM design matrix (lpmatrix) so each draw's beta.f(hydro) is
  ##    `lp %*% coef` (identical to predict.gam, but no per-draw basis rebuild),
  ##  - the per-(analyte,date) normalisation factor (co_split is exact in the
  ##    C_raw reconstruction, so the factor is constant across draws).
  co_by_date <- stats::setNames(lapply(as.character(qdates), function(ds) {
    cd <- co_split[[ds]]
    if (is.null(cd)) numeric(0) else stats::setNames(cd$value, cd$analyte)
  }), as.character(qdates))

  static <- stats::setNames(lapply(modelled, function(nm) {
    m <- tm$models[[nm]]
    feats <- .compute_hydro_features(
      tm$hydro, qdates, m$tau_short,
      m$tau_long, tm$hydro_type
    )
    ref_lk <- stats::setNames(
      ref_q$ref_norm[ref_q$analyte == nm],
      as.character(ref_q$date[ref_q$analyte == nm])
    )
    ref_vec <- as.numeric(ref_lk[as.character(qdates)])
    ref_vec[is.na(ref_vec)] <- 0
    lp <- NULL
    if (m$tier == "model" && !is.null(m$impact_fit)) {
      nd <- dplyr::tibble(
        hydro_short = feats$hydro_short,
        hydro_long = feats$hydro_long
      )
      if (isTRUE(m$pooled)) nd$analyte <- factor(m$analyte, levels = m$pool_levels)
      lp <- tryCatch(stats::predict(m$impact_fit, newdata = nd, type = "lpmatrix"),
        error = function(e) NULL
      )
    }
    parsed <- fac_lookup[[nm]]$parsed
    fac <- if (is.null(parsed)) {
      rep(1, length(qdates))
    } else {
      vapply(seq_along(qdates), function(i) {
        f <- .apply_normalisation(parsed, 1, co_by_date[[i]])
        if (is.na(f) || f <= 0) NA_real_ else f
      }, numeric(1L))
    }
    names(fac) <- as.character(qdates)
    list(
      hydro_short = feats$hydro_short, hydro_long = feats$hydro_long,
      ref_vec = ref_vec, lp = lp, pooled = isTRUE(m$pooled),
      pool_center = m$pool_center, pool_scale = m$pool_scale, fac = fac
    )
  }), modelled)

  ## Flat (analyte date) -> factor lookup for vectorised C_raw reconstruction.
  fac_lut <- unlist(lapply(modelled, function(nm) {
    f <- static[[nm]]$fac
    stats::setNames(as.numeric(f), paste(nm, names(f)))
  }))

  list(
    tm           = tm,
    modelled     = modelled,
    qdates       = qdates,
    co_split     = co_split,
    wq_long      = wq_long,
    fac_lookup   = fac_lookup,
    measured_key = measured_key,
    smoothers    = smoothers,
    kappa        = kappa,
    ou_scale     = ou_scale,
    ref_q        = ref_q,
    static       = static,
    fac_lut      = fac_lut,
    co_grab_map  = co_grab_map
  )
}


#' Per-draw daily toxicant prediction from a pre-fitted scaffold (issue #16)
#'
#' Predicts normalised and raw concentrations for one draw (or for the
#' deterministic centre line when `residual_paths = NULL`).  Always uses the
#' precomputed static scaffolding from `fdm`; accepts an optionally-perturbed
#' target model (`tm_p`) and a per-analyte residual path (`residual_paths`) for
#' draw mode.  Pass overrides to `co_split` and `wq_long` (from
#' [.perturb_co_split()]) for S7 co-analyte measurement-error draws.
#'
#' @param fdm Output of [.fit_daily_target()].
#' @param tm_p Target model to predict with (default: `fdm$tm`).  For draw
#'   mode, pass a GAM-perturbed copy from [.perturb_target_model()].
#' @param residual_paths Named list of residual paths (impact `S` or WQ `d`),
#'   one numeric vector of length `length(fdm$qdates)` per analyte (`NA` outside
#'   the analyte's clipped grab span).  `NULL` (the default) uses the smoother
#'   posterior means — the deterministic centre line.
#' @param co_split Per-date co-analyte lookup (default: `fdm$co_split`).
#'   Chunk D supplies a perturbed version for S7 draws.
#' @param wq_long WQ layer data (default: `fdm$wq_long`).
#' @return Tibble of model rows (`.date`, `value`, `detected`, `.measured`,
#'   `analyte`, `units.analyte`) with `attr("impact_tiers")` attached, or
#'   `NULL` if prediction produced no finite rows.
#' @keywords internal
.predict_daily_tox <- function(fdm,
                               tm_p = fdm$tm,
                               residual_paths = NULL,
                               co_split = fdm$co_split,
                               wq_long = fdm$wq_long) {
  ## Residual paths (S for impact tier, d for WQ tier) on qdates, NA outside the
  ## per-analyte clipped grab span. NULL -> deterministic centre line (smoother
  ## posterior means from the fit-once scaffold).
  if (is.null(residual_paths)) {
    residual_paths <- stats::setNames(lapply(fdm$modelled, function(nm) {
      sm <- fdm$smoothers[[nm]]
      if (is.null(sm) || length(sm$grid_dates) == 0L) {
        return(rep(NA_real_, length(fdm$qdates)))
      }
      .residual_on_qdates(sm$grid_dates, sm$mean, fdm$qdates)
    }), fdm$modelled)
  }

  pred <- .resolve_target_impact(tm_p, tibble::tibble(date = fdm$qdates),
    fdm$modelled,
    wq = wq_long,
    residual_paths = residual_paths,
    kappa = fdm$kappa, scale = fdm$ou_scale,
    ref_q = fdm$ref_q, static = fdm$static
  )
  if (nrow(pred) == 0L) {
    return(NULL)
  }

  ## Reconstruct raw µg/L: C_raw = C_norm / normalisation_factor(co-analytes).
  ## co_split is exact in this reconstruction (S7 enters via the synthetic
  ## co-rows), so the factor is constant across draws: use the precomputed flat
  ## (analyte date) -> factor lookup. Fall back to per-row evaluation only when
  ## the static context is unavailable (standalone callers).
  if (!is.null(fdm$fac_lut)) {
    key <- paste(pred$analyte, as.character(pred$date))
    pred$C_raw <- pred$C_norm / as.numeric(fdm$fac_lut[key])
  } else {
    co_vec_for <- function(d) {
      cd <- co_split[[as.character(d)]]
      if (is.null(cd)) {
        return(numeric(0))
      }
      stats::setNames(cd$value, cd$analyte)
    }
    pred$C_raw <- vapply(seq_len(nrow(pred)), function(i) {
      a <- pred$analyte[i]
      d <- pred$date[i]
      parsed <- fdm$fac_lookup[[a]]$parsed
      if (is.null(parsed)) {
        return(pred$C_norm[i])
      }
      factor <- .apply_normalisation(parsed, 1, co_vec_for(d))
      if (is.na(factor) || factor <= 0) {
        return(NA_real_)
      }
      pred$C_norm[i] / factor
    }, numeric(1L))
  }

  pred_ok <- pred[is.finite(pred$C_raw), , drop = FALSE]
  if (nrow(pred_ok) == 0L) {
    return(NULL)
  }

  model_rows <- tibble::tibble(
    .date         = pred_ok$date,
    value         = pred_ok$C_raw,
    detected      = TRUE,
    .measured     = paste(pred_ok$analyte, pred_ok$date) %in% fdm$measured_key,
    analyte       = pred_ok$analyte,
    units.analyte = "ug/L"
  )
  attr(model_rows, "impact_tiers") <- dplyr::distinct(
    pred[, c("analyte", "impact_tier")]
  )
  model_rows
}


## ── Chunk D: co-analyte measurement-error perturbation helpers (S7) ──────────
## (S6 grab measurement error now enters as the Kalman observation noise r_i in
## .fit_daily_target(); there is no separate anchor-perturbation step.)


#' Draw coherent co-analyte perturbations for one draw (S7)
#'
#' For each co-analyte grab in `fdm`, draws one lognormal multiplier
#' (geo-mean = 1, per-analyte or global CV from `grab_cv_co`).  Each day in
#' `co_split` inherits the multiplier of its source grab (from
#' `fdm$co_grab_map`), preserving temporal coherence: consecutive forward-
#' filled days that came from the same grab receive the same perturbation.
#'
#' @param fdm Fitted daily scaffold from [.fit_daily_target()].
#' @param grab_cv_co Named numeric vector of CVs per co-analyte, or a single
#'   scalar.  `NULL` or `NA` → return `fdm$co_split` and `fdm$wq_long` unchanged.
#' @return Named list `list(co_split, wq_long)` with perturbed values.
#' @keywords internal
.perturb_co_split <- function(fdm, grab_cv_co) {
  if (is.null(fdm$co_grab_map) || is.null(grab_cv_co) ||
    (length(grab_cv_co) == 1L && is.na(grab_cv_co))) {
    return(list(co_split = fdm$co_split, wq_long = fdm$wq_long))
  }

  get_cv_co <- function(a) {
    if (length(grab_cv_co) == 1L) {
      return(as.numeric(grab_cv_co))
    }
    if (a %in% names(grab_cv_co)) {
      return(as.numeric(grab_cv_co[[a]]))
    }
    return(NA_real_)
  }

  co_analytes_nm <- names(fdm$co_grab_map)

  ## One multiplier per unique (analyte, source-grab-date) pair.
  all_pairs <- unique(unlist(lapply(co_analytes_nm, function(a) {
    src <- fdm$co_grab_map[[a]]
    src_ok <- src[!is.na(src)]
    paste(a, src_ok, sep = "::")
  })))

  if (length(all_pairs) == 0L) {
    return(list(co_split = fdm$co_split, wq_long = fdm$wq_long))
  }

  analyte_of_pair <- sub("::.*", "", all_pairs)
  sigma_lns <- vapply(analyte_of_pair, function(a) {
    cv <- get_cv_co(a)
    if (is.na(cv) || cv <= 0) 0 else sqrt(log(1 + cv^2))
  }, numeric(1L))

  mults <- stats::setNames(
    exp(stats::rnorm(length(all_pairs), -sigma_lns^2 / 2, sigma_lns)),
    all_pairs
  )

  ## Apply multipliers to co_split (one date at a time).
  co_split_d <- lapply(names(fdm$co_split), function(d) {
    co_day <- fdm$co_split[[d]]
    for (i in seq_len(nrow(co_day))) {
      a <- co_day$analyte[i]
      cgm <- fdm$co_grab_map[[a]]
      if (is.null(cgm)) next
      src_d <- cgm[d]
      if (is.na(src_d)) next
      key <- paste(a, src_d, sep = "::")
      m <- mults[[key]]
      if (!is.null(m) && is.finite(m) && m > 0) co_day$value[i] <- co_day$value[i] * m
    }
    co_day
  })
  names(co_split_d) <- names(fdm$co_split)

  ## Rebuild wq_long (needed for WQ-layer PC-score computation) if present.
  wq_long_d <- if (!is.null(fdm$wq_long)) {
    rows <- lapply(names(co_split_d), function(d) {
      co_d <- co_split_d[[d]]
      if (nrow(co_d) == 0L) {
        return(NULL)
      }
      data.frame(
        sample_id = d, analyte = co_d$analyte, value = co_d$value,
        stringsAsFactors = FALSE
      )
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) > 0L) tibble::as_tibble(do.call(rbind, rows)) else NULL
  } else {
    NULL
  }

  list(co_split = co_split_d, wq_long = wq_long_d)
}


#' Build per-draw co-analyte rows for the synthetic frame (S7, G2)
#'
#' Creates a copy of the exact co-analyte rows from `daily_long_exact` with
#' values replaced by the corresponding perturbed values from `co_pert_split`
#' (output of [.perturb_co_split()]), and tags the result with `draw_id`.
#' Analytes absent from `co_pert_split` (e.g. non-modelled toxicants) are left
#' unchanged.  The resulting rows go into the synthetic frame so that
#' [add_mspaf()]'s normalisation sees the per-draw co-analyte perturbation.
#'
#' @param daily_long_exact Rows of the forward-filled daily chemistry that
#'   belong to the non-modelled subset (co-analytes + unmodelled toxicants).
#' @param co_pert_split Named list (date-string → data.frame(analyte, value))
#'   from [.perturb_co_split()].
#' @param draw_id Integer draw index to assign to the output rows.
#' @return A copy of `daily_long_exact` with perturbed values and `draw_id`
#'   set to `draw_id`.
#' @keywords internal
.co_draw_rows <- function(daily_long_exact, co_pert_split, draw_id) {
  ## Flatten the per-date perturbed co-analyte values into a single lookup
  ## table, then join once — avoids the row-by-row loop over daily_long_exact.
  pert_tbl <- dplyr::bind_rows(
    lapply(names(co_pert_split), function(d) {
      cd <- co_pert_split[[d]]
      if (is.null(cd) || nrow(cd) == 0L) {
        return(NULL)
      }
      tibble::tibble(.date_chr = d, analyte = cd$analyte, .new_val = cd$value)
    })
  )

  result <- daily_long_exact
  if (nrow(pert_tbl) > 0L) {
    result$.date_chr <- as.character(result$.date)
    result <- dplyr::left_join(result, pert_tbl, by = c(".date_chr", "analyte"))
    idx_replace <- !is.na(result$.new_val)
    result$value[idx_replace] <- result$.new_val[idx_replace]
    result$.date_chr <- NULL
    result$.new_val <- NULL
  }

  result$draw_id <- as.integer(draw_id)
  result
}


#' Model-based daily toxicant interpolation (season-blind impact model)
#'
#' Thin wrapper: calls [.fit_daily_target()] once then [.predict_daily_tox()]
#' once (point mode, no ε).  Modelled-toxicant rows in `daily_long` are
#' replaced; co-analytes and non-modelled toxicants are left untouched.  On
#' any failure the input is returned unchanged.
#'
#' For draw-mode orchestration (Chunk E), call the two underlying helpers
#' directly: fit once with `.fit_daily_target()`, then loop N times over
#' `.predict_daily_tox(fdm, perturbed_tm, eps_paths)`.
#'
#' @param daily_long Output of `.build_daily_chem()` (+ temperature fill).
#' @param site_rows Grab chemistry for the site (passed to the target model).
#' @param reference_model A `reference_model`.
#' @param imputation_model Optional `imputation_model` (tier-2 enrichment).
#' @param conc_units Units for `site_rows` toxicants when no `units.analyte`.
#' @param meta Analyte metadata (normalisation formulas).
#' @param tox_analytes SSD-eligible analyte names.
#' @return `daily_long` with modelled-toxicant rows replaced.
#' @keywords internal
.daily_tox_from_model <- function(daily_long, site_rows, reference_model,
                                  imputation_model, conc_units, meta,
                                  tox_analytes, method = "multi",
                                  guideline_dir = NULL, transform = "pseudo_log") {
  fdm <- .fit_daily_target(
    site_rows        = site_rows,
    reference_model  = reference_model,
    imputation_model = imputation_model,
    conc_units       = conc_units,
    meta             = meta,
    tox_analytes     = tox_analytes,
    daily_long       = daily_long,
    method           = method,
    guideline_dir    = guideline_dir,
    transform        = transform
  )
  if (is.null(fdm)) {
    return(daily_long)
  }

  model_rows <- .predict_daily_tox(fdm)
  if (is.null(model_rows)) {
    return(daily_long)
  }

  modelled <- fdm$modelled
  keep <- dplyr::filter(daily_long, !.data$analyte %in% .env$modelled)
  if (!"units.analyte" %in% names(keep)) keep$units.analyte <- NA_character_
  out <- dplyr::bind_rows(keep, model_rows)
  attr(out, "impact_tiers") <- attr(model_rows, "impact_tiers")
  out
}


#' Empty tibble matching the mspaf_daily() return schema
#' @param mode One of `"point"`, `"summary"`, or `"draws"` — governs which
#'   extra columns are included.
#' @param gap_uncertainty Bracket mode (`"bracket"`/`"ignorable"`/
#'   `"informative"`); governs the envelope columns in draws/summary mode.
#'   Ignored for point mode.
#' @keywords internal
.empty_daily_result <- function(mode = c("point", "summary", "draws"),
                                gap_uncertainty = "bracket") {
  mode <- match.arg(mode)
  diag <- tibble::tibble(
    n_analytes_used        = integer(),
    dominant_analyte       = character(),
    max_paf                = numeric(),
    n_measured_analytes    = integer(),
    days_since_last_sample = integer(),
    analyte_pafs           = list()
  )
  if (mode == "point") {
    return(dplyr::bind_cols(
      tibble::tibble(
        date = as.Date(character()), site_id = character(),
        mspaf = numeric()
      ),
      diag
    ))
  }
  if (mode == "draws") {
    env <- stats::setNames(
      rep(list(numeric()), length(.bracket_draw_cols(gap_uncertainty))),
      .bracket_draw_cols(gap_uncertainty)
    )
    return(dplyr::bind_cols(
      tibble::tibble(
        date = as.Date(character()), site_id = character(),
        draw_id = integer()
      ),
      tibble::as_tibble(env),
      diag
    ))
  }
  ## summary
  env <- stats::setNames(
    rep(list(numeric()), length(.bracket_summary_cols(gap_uncertainty))),
    .bracket_summary_cols(gap_uncertainty)
  )
  dplyr::bind_cols(
    tibble::tibble(
      date = as.Date(character()), site_id = character(),
      deterministic = numeric()
    ),
    tibble::as_tibble(env),
    diag
  )
}
