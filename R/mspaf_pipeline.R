# ── mspaf_pipeline() ──────────────────────────────────────────────────────────

#' Run the full chronic msPAF pipeline end to end
#'
#' Orchestrates the four steps a chronic msPAF analysis usually chains by hand:
#' optionally fit a Bayesian imputation model ([fit_imputation_model()]), fit
#' the reference/impact model ([fit_reference_model()]), compute the daily
#' multi-substance PAF ([mspaf_daily()]) for the target site, then aggregate the
#' daily series into a **chronic** time-weighted msPAF
#' ([time_weighted_aggregate()]). The fitted imputation model is threaded into
#' **both** downstream chemistry steps so that below-detection and
#' entirely-absent analytes are filled before the toxicity calculation.
#'
#' The chronic step "smears" the daily msPAF line with an exponential-decay
#' memory kernel: because the daily series is equally spaced, forward-step
#' duration weighting is uniform (Δt = 1 day) and the chronic value is shaped
#' only by the `tau` decay. Set `chronic = FALSE` to stop at the daily series
#' and reproduce the previous (daily) return.
#'
#' **Imputation is on by default here.** The individual functions treat
#' imputation as opt-in (you must fit a model and pass it). This orchestrator
#' flips that default: with `impute = TRUE` it fits an imputation model and
#' applies it throughout. This is convenient but **not free of consequences** —
#' imputing below-detection cells (replacing detection-limit values with
#' modelled sub-limit estimates) and fabricating entirely-absent analytes
#' generally *moves* the msPAF result. Whether that is appropriate is a
#' missingness (MAR/MNAR) modelling judgement, not a neutral default. Set
#' `impute = FALSE` to reproduce the non-imputed behaviour of calling
#' [fit_reference_model()] and [mspaf_daily()] directly.
#'
#' With `reference = NULL` (no reference site of the right water type) no
#' reference model is fitted. `daily_args$interpolation = "model"` is then an
#' error (class `hydroSense_model_needs_reference`); the default is
#' `"linear"`. All such argument errors are raised before any model is fitted.
#'
#' @param target Long-format target-site chemistry (same schema as
#'   [mspaf_daily()]'s `df`): `sample_id`, `site_id`, `datetime`, `analyte`,
#'   `value`, `detected`, and a unit column (`units.analyte`, or sampleTidy's
#'   `units`). Analyte names may be hydroSense's or sampleTidy's (see
#'   [analyte_aliases()]).
#' @param reference Long-format reference-site chemistry passed to
#'   [fit_reference_model()]. `NULL` (default) runs the pipeline without a
#'   reference site: no reference model is fitted (the `reference_model`
#'   attribute is `NULL`), `ara` defaults to `"none"`, daily interpolation
#'   defaults to `"linear"` and `impute_on` defaults to `"target"`. See
#'   Details.
#' @param hydro Optional hydrology series forwarded to [fit_reference_model()].
#' @param impute Logical. Fit and apply a Bayesian imputation model? Default
#'   `TRUE`. When `FALSE`, no imputation model is fitted and `NULL` is passed
#'   downstream (identical to chaining the functions by hand without one).
#' @param impute_on Which chemistry to fit the imputation model on:
#'   `"reference"` (default, matching the documented design) or `"target"`.
#'   With `reference = NULL` the default becomes `"target"`; an explicit
#'   `"reference"` is an error of class `hydroSense_impute_needs_reference`
#'   when `impute = TRUE`.
#' @param required_vars Passed to [fit_imputation_model()]. Default
#'   `c("pH", "EC")`.
#' @param impute_groups Optional `groups` for [fit_imputation_model()]; `NULL`
#'   (default) uses [leachate_impute_groups()].
#' @param impute_seed Optional integer seed for the imputation fit (forwarded to
#'   [fit_imputation_model()]) so a stochastic fit is reproducible. The default
#'   `"marginal"` method's GAM fit is deterministic; the seed matters for the
#'   Stan-based methods (`"factor"` and the brms methods).
#' @param reference_args,daily_args Named lists of additional arguments
#'   forwarded to [fit_reference_model()] and [mspaf_daily()] respectively. Do
#'   not include arguments the orchestrator sets itself (`reference`/`hydro`/
#'   `imputation_model` for the reference fit; `df`/`reference_model`/
#'   `imputation_model` for the daily call). `daily_args` may override the
#'   `interpolation = "model"` default. When `chronic = TRUE` and `daily_args`
#'   requests draws (`ndraws`), `return = "draws"` is forced so per-draw rows
#'   can be propagated into the chronic summary; the `gap_uncertainty`,
#'   `interval` and `central` values from `daily_args` (or [mspaf_daily()]'s
#'   defaults `"bracket"`, `0.9`, `"median"`) shape the chronic envelope.
#' @param chronic Logical. Aggregate the daily series into chronic time-weighted
#'   msPAF as a final step? Default `TRUE`. When `FALSE`, the daily
#'   [mspaf_daily()] result is returned unchanged (with model attributes only).
#' @param focal_dates Focal dates for the chronic aggregation, passed to
#'   [time_weighted_aggregate()]. `NULL` (default) derives a sequence spanning
#'   the daily date range via [expand_focal_dates()] (spacing `focal_by`). May
#'   also be a `Date` vector or a `focal_date`/`site_id` data frame.
#' @param focal_by Spacing for the derived `focal_dates` when `focal_dates` is
#'   `NULL`; passed to [expand_focal_dates()]'s `by` ([base::seq.Date()]
#'   semantics): `"day"` (default), `7` (weekly), `"week"`, `"month"`, etc.
#'   Ignored when `focal_dates` is supplied.
#' @param tau Exponential-decay parameter in days for the chronic kernel,
#'   forwarded to [time_weighted_aggregate()]. `NULL` (default) uses its default
#'   of 90 days.
#' @param window Look-back window in days for the chronic aggregation, forwarded
#'   to [time_weighted_aggregate()]. `NULL` (default) uses its default of 365
#'   days.
#' @param chronic_summary Aggregation method for the chronic step, passed to
#'   [time_weighted_aggregate()]'s `summary`. Default `"arith_mean"` (the
#'   recommended aggregation for bounded msPAF percentages).
#' @param also_total Logical (default `FALSE`). Point mode only. Passed to
#'   [mspaf_daily()], which adds a total-concentration (no ARA) `mspaf_total`
#'   daily series alongside the ARA one; the chronic step aggregates it exactly
#'   as it aggregates `mspaf` and returns it as `value_total`. The static ARA
#'   reference is `daily_args$reference` (the `reference` argument here is the
#'   reference-model fit data); without it `value_total` equals `value`. An
#'   error in draws mode (`daily_args$ndraws` set). Composes with every `ara`
#'   mode: `value_total` is always the no-ARA value.
#' @param ara Area-reference adjustment (ARA) applied by [mspaf_daily()]:
#'   - `"static"`: static ARA against the per-analyte summary of
#'     `daily_args$reference` (a raw reference chemistry frame, see
#'     [prepare_reference()]), which is then required.
#'   - `"contemporaneous"`: the reference model this pipeline fits is passed
#'     to [mspaf_daily()] as `reference`, so each day is adjusted by the
#'     reference model's resolution for that day (tier-1 match, GAM
#'     prediction or static fallback; see [fit_reference_model()]). The model
#'     is fitted once and reused. `daily_args$reference` must not be given.
#'   - `"none"`: no ARA. `daily_args$reference` must not be given.
#'
#'   With `reference = NULL`, `"contemporaneous"` (and `"static"` without
#'   `daily_args$reference`) is an error of class
#'   `hydroSense_ara_needs_reference`.
#'
#'   When `ara` is not supplied the behaviour is unchanged from earlier
#'   versions: static ARA if `daily_args$reference` is given, otherwise none.
#'   Conflicting combinations are errors raised before any model is fitted.
#'
#' @return When `chronic = TRUE` (default), the chronic msPAF frame, one row per
#'   (focal date × site):
#'   - **point mode** (`ndraws` absent): columns `focal_date`, `site_id`,
#'     `sample_id`, `analyte`, `value`, `detected`, `n_samples_in_window`,
#'     `n_imputed_in_window` (a single chronic `value`, no interval columns),
#'     plus `value_total` when `also_total = TRUE`.
#'   - **draws mode** (`ndraws` set): columns `focal_date`, `site_id` plus the
#'     bracket envelope summary (`median_*`/`lo_*`/`hi_*` per envelope, and
#'     `precautionary_lo`/`precautionary_hi` in `"bracket"` mode).
#'   The chronic frame carries the daily frame it was built from as attribute
#'   `"daily"`. When `chronic = FALSE`, the [mspaf_daily()] result is returned
#'   unchanged. In all cases the fitted models are attached as attributes
#'   `"reference_model"` and `"imputation_model"` (the latter `NULL` when
#'   `impute = FALSE` or the fit was skipped). Note that most \pkg{dplyr} verbs
#'   drop attributes, so read them off the returned frame directly.
#'
#' @seealso [fit_imputation_model()], [fit_reference_model()], [mspaf_daily()],
#'   [time_weighted_aggregate()], [expand_focal_dates()]
#' @examples
#' \dontrun{
#' demo <- leachate_demo()
#' out <- mspaf_pipeline(
#'   target = subset(demo, site_id == "downstream"),
#'   reference = subset(demo, site_id == "reference"),
#'   daily_args = list(require_temperature = FALSE, ndraws = 50L, seed = 1L)
#' )
#' }
#' @export
mspaf_pipeline <- function(target,
                           reference = NULL,
                           hydro = NULL,
                           impute = TRUE,
                           impute_on = c("reference", "target"),
                           required_vars = c("pH", "EC"),
                           impute_groups = NULL,
                           impute_seed = NULL,
                           reference_args = list(),
                           daily_args = list(),
                           chronic = TRUE,
                           focal_dates = NULL,
                           focal_by = "day",
                           tau = NULL,
                           window = NULL,
                           chronic_summary = "arith_mean",
                           also_total = FALSE,
                           ara = c("static", "contemporaneous", "none")) {
  checkmate::assert_data_frame(target)
  checkmate::assert_data_frame(reference, null.ok = TRUE)
  no_ref <- is.null(reference)
  checkmate::assert_flag(impute)
  impute_on_default <- missing(impute_on)
  impute_on <- match.arg(impute_on)
  if (no_ref && impute_on_default) impute_on <- "target"
  checkmate::assert_list(reference_args)
  checkmate::assert_list(daily_args)
  checkmate::assert_flag(chronic)
  checkmate::assert_string(chronic_summary)
  checkmate::assert_flag(also_total)
  # Fail before any model is fitted rather than inside mspaf_daily().
  if (also_total && !is.null(daily_args$ndraws)) {
    cli::cli_abort(c(
      "{.code also_total = TRUE} is supported in point mode only.",
      "i" = "Drop {.code daily_args$ndraws}."
    ), class = "hydroSense_also_total_draws")
  }
  # Resolve the ARA mode, also before any fitting. Without `ara`, the presence
  # of daily_args$reference decides (static or none), exactly as before the
  # argument existed.
  has_static_ref <- !is.null(daily_args$reference)
  ara <- if (missing(ara)) {
    if (has_static_ref) "static" else "none"
  } else {
    match.arg(ara)
  }
  if (no_ref) {
    if (ara == "contemporaneous" || (ara == "static" && !has_static_ref)) {
      cli::cli_abort(c(
        "{.code ara = \"{ara}\"} needs a reference site.",
        "i" = "With {.code reference = NULL} use {.code ara = \"none\"}, or \
               {.code ara = \"static\"} with {.code daily_args$reference}."
      ), class = "hydroSense_ara_needs_reference")
    }
    if (identical(daily_args$interpolation, "model")) {
      cli::cli_abort(c(
        "{.code daily_args$interpolation = \"model\"} needs a reference \
         model.",
        "i" = "With {.code reference = NULL} use {.code \"linear\"} \
               (the default) or another non-model interpolation."
      ), class = "hydroSense_model_needs_reference")
    }
    if (impute && impute_on == "reference") {
      cli::cli_abort(c(
        "{.code impute_on = \"reference\"} needs a reference site.",
        "i" = "Use {.code impute_on = \"target\"} or {.code impute = FALSE}."
      ), class = "hydroSense_impute_needs_reference")
    }
  }
  if (ara == "static" && !has_static_ref) {
    cli::cli_abort(c(
      "{.code ara = \"static\"} needs a raw reference frame in \\
       {.code daily_args$reference}.",
      "i" = "Use {.code ara = \"contemporaneous\"} to adjust against the \\
             fitted reference model instead."
    ), class = "hydroSense_ara_static_needs_reference")
  }
  if (ara != "static" && has_static_ref) {
    cli::cli_abort(c(
      "{.code ara = \"{ara}\"} conflicts with {.code daily_args$reference}.",
      "i" = "Drop {.code daily_args$reference}, or use \\
             {.code ara = \"static\"}."
    ), class = "hydroSense_ara_conflict")
  }

  # Accept sampleTidy's analyte names and `units` column.
  target    <- .canonicalise_mspaf_input(target)
  if (!no_ref) reference <- .canonicalise_mspaf_input(reference)

  # ── 1. (Optional) imputation model — fit once, thread into both steps ───────
  imodel <- NULL
  if (impute) {
    fit_chem <- if (impute_on == "reference") reference else target
    imp_args <- list(fit_chem,
      required_vars = required_vars, groups = impute_groups
    )
    if (!is.null(impute_seed)) imp_args$seed <- impute_seed
    imodel <- tryCatch(
      do.call(fit_imputation_model, imp_args),
      error = function(e) {
        cli::cli_warn(c(
          "Imputation skipped: {conditionMessage(e)}",
          "i" = "Proceeding without imputation."
        ))
        NULL
      }
    )
  }

  # ── 2. Reference / impact model ─────────────────────────────────────────────
  ref_model <- if (no_ref) {
    NULL
  } else {
    do.call(
      fit_reference_model,
      c(list(reference, hydro = hydro, imputation_model = imodel),
        reference_args)
    )
  }

  # ── 3. Daily msPAF on the target ────────────────────────────────────────────
  daily_call <- utils::modifyList(
    list(interpolation = if (no_ref) "linear" else "model"), daily_args
  )
  draws_mode <- !is.null(daily_call$ndraws)
  # When chronic && draws_mode, force per-draw rows so they can be propagated
  # into the chronic summary; otherwise honour daily_call$return.
  if (chronic && draws_mode) daily_call$return <- "draws"
  if (also_total) daily_call$also_total <- TRUE
  # Contemporaneous ARA: the reference model fitted above is the reference.
  if (ara == "contemporaneous") daily_call$reference <- ref_model
  daily <- do.call(
    mspaf_daily,
    c(
      list(target, reference_model = ref_model, imputation_model = imodel),
      daily_call
    )
  )

  if (!chronic) {
    attr(daily, "reference_model") <- ref_model
    attr(daily, "imputation_model") <- imodel
    return(daily)
  }

  # ── 4. Chronic time-weighted aggregation of the daily series ────────────────
  focal <- focal_dates %||%
    expand_focal_dates(min(daily$date), max(daily$date), by = focal_by)

  gap_uncertainty <- daily_call$gap_uncertainty %||% "bracket"
  interval        <- daily_call$interval %||% 0.9
  central         <- daily_call$central %||% "median"

  out <- .chronic_from_daily(
    daily           = daily,
    draws_mode      = draws_mode,
    gap_uncertainty = gap_uncertainty,
    focal           = focal,
    tau             = tau,
    window          = window,
    chronic_summary = chronic_summary,
    interval        = interval,
    central         = central,
    also_total      = also_total
  )

  attr(out, "daily")            <- daily
  attr(out, "reference_model")  <- ref_model
  attr(out, "imputation_model") <- imodel
  out
}

# ── .chronic_from_daily() ─────────────────────────────────────────────────────

#' Aggregate a daily msPAF frame into chronic time-weighted msPAF
#'
#' Reshapes the [mspaf_daily()] output into the long `analyte = "msPAF"` schema
#' [time_weighted_aggregate()] requires, then aggregates. In point mode the daily
#' `mspaf` line is smeared into a single chronic value per (focal date × site);
#' in draws mode each envelope column (`mspaf_ignorable` / `mspaf_informative`)
#' is propagated per draw and collapsed with [.summarise_bracket()]. With
#' `also_total` (point mode) the daily `mspaf_total` line is aggregated the
#' same way and joined as `value_total`.
#'
#' @keywords internal
.chronic_from_daily <- function(daily, draws_mode, gap_uncertainty, focal,
                                tau, window, chronic_summary, interval,
                                central, also_total = FALSE) {
  # Reshape a daily frame to the long msPAF schema, taking the value from
  # `value_col` and carrying `draw_id` when present.
  to_long <- function(value_col) {
    long <- tibble::tibble(
      sample_id = paste0("daily_", daily$site_id, "_", daily$date),
      site_id   = daily$site_id,
      datetime  = daily$date,
      analyte   = "msPAF",
      detected  = TRUE,
      value     = daily[[value_col]]
    )
    if (draws_mode) long$draw_id <- daily$draw_id
    long
  }

  if (!draws_mode) {
    point_agg <- function(value_col) {
      time_weighted_aggregate(
        to_long(value_col),
        focal_dates  = focal,
        tau          = tau,
        tau_units    = "d",
        window       = window,
        window_units = "d",
        summary      = chronic_summary,
        return       = "summary"
      )
    }
    out <- point_agg("mspaf")
    if (also_total) {
      tot <- dplyr::select(point_agg("mspaf_total"),
        "focal_date", "site_id", value_total = "value"
      )
      out <- dplyr::left_join(out, tot, by = c("focal_date", "site_id"))
      out <- dplyr::relocate(out, "value_total", .after = "value")
    }
    return(out)
  }

  # Draws mode: aggregate each present envelope's draws, then bracket-summarise.
  env_cols <- intersect(.bracket_draw_cols(gap_uncertainty), names(daily))

  per_env <- lapply(env_cols, function(col) {
    long <- to_long(col)
    e <- time_weighted_aggregate(
      long,
      focal_dates  = focal,
      tau          = tau,
      tau_units    = "d",
      window       = window,
      window_units = "d",
      summary      = chronic_summary,
      return       = "draws"
    )
    e <- dplyr::select(
      e, "focal_date", "site_id", "draw_id", "value"
    )
    dplyr::rename(e, !!col := "value")
  })

  draws_wide <- Reduce(
    function(a, b) {
      dplyr::full_join(a, b, by = c("focal_date", "site_id", "draw_id"))
    },
    per_env
  )

  draws_wide <- dplyr::rename(draws_wide, date = "focal_date")
  summ <- .summarise_bracket(
    draws_wide,
    interval        = interval,
    central         = central,
    gap_uncertainty = gap_uncertainty
  )
  dplyr::rename(summ, focal_date = "date")
}
