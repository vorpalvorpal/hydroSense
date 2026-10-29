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
#' @param contributions Logical (default `FALSE`). Point mode only. When
#'   `TRUE`, the result carries attribute `"contributions"`: the chronic msPAF
#'   apportioned to analytes and mode-of-action groups, with columns
#'   `site_id`, `focal_date`, `analyte`, `group` and `contribution` (%).
#'   Each analyte's daily contribution series (from [mspaf_daily()]; zero on
#'   days it did not enter the msPAF) is aggregated with the same time weights,
#'   window and rows as the chronic `value`, so the contributions sum to
#'   `value` for every (site, focal date). With `also_total = TRUE` a
#'   `variant` column (`"ara"` / `"total"`) distinguishes the contributions to
#'   `value` and `value_total`. Requires `chronic_summary = "arith_mean"` (the
#'   only additive aggregation) when `chronic = TRUE`. An error of class
#'   `hydroSense_contributions_draws` in draws mode (`daily_args$ndraws` set).
#'   With `chronic = FALSE` the daily breakdown from [mspaf_daily()] is
#'   attached to the daily result instead. For `ara_form = "added_hazard"`
#'   the definition differs; see the section *Added-hazard form*.
#' @param ara_form The form of the area-reference adjustment:
#'   - `"excess"` (default): the existing behaviour described above, where the
#'     reference is subtracted from each analyte's concentration before the
#'     msPAF is formed.
#'   - `"added_hazard"`: the reference is compared at the msPAF level, as the
#'     added hazard of the target over the reference (see the section
#'     *Added-hazard form*). It has its own engine: no model is fitted and
#'     nothing is imputed.
#' @param min_coverage Used only by `ara_form = "added_hazard"`: the smallest
#'   share (0 to 1) of the time weights that covered days must carry for a
#'   chronic value to be reported. Default `0.2`.
#'
#' @section Added-hazard form:
#' With `ara_form = "added_hazard"` each point (every target site, and the
#' reference) is built the same way from its own samples, on a daily grid.
#'
#' **Input.** `target` and `reference` are raw long frames in the same format
#' (`sample_id`, `site_id`, `datetime`, `analyte`, `value`, `detected`, a unit
#' column). A below-detection row has `detected = FALSE` and its detection
#' limit (LOR) as `value`. Analyte names go through [analyte_aliases()].
#' Toxicants are the analytes with SSD parameters for the given
#' `daily_args$analyte_metadata` / `method` / `guideline_dir`; the rest are
#' co-analytes. Rows are reduced to one value per (site, analyte, day): the
#' mean of the detected results when there is one, otherwise the smallest LOR
#' (below detection). `reference` must hold exactly one `site_id` (error
#' class `hydroSense_reference_sites`).
#'
#' **Grid.** Daily, from the earliest sample date of the targets and the
#' reference to the latest target sample date plus the span cap (the latest
#' date itself when the cap is infinite), unless `daily_args$start` / `end`
#' are given.
#'
#' **Point build.** Toxicants are interpolated with [interpolate_loglinear()]
#' and the span cap `daily_args$span_cap` (default 120 days): days in a gap
#' between samples longer than the cap, more than the cap after the last
#' sample, or before the first sample are uncovered. Below-detection results
#' are anchored at the LOR for a target and at LOR/2 for the reference. The
#' co-analytes are interpolated without a cap and held constant beyond the
#' first and last sample (detected positive results only; one result is a
#' constant): pH and temperature linearly, DOC, hardness, Ca and Mg linearly
#' in log. Hardness is the measured value where it exists, else
#' `2.497 Ca + 4.118 Mg` on days with both. `daily_args$temperature` (a
#' `datetime`/`value` series, per site only when it has a `site_id` column)
#' replaces every point's own temperature results when given; it is
#' interpolated linearly between its days and held beyond its ends. Each
#' toxicant is normalised with the point's own co-analytes using the metadata
#' normalisation formulas, as [add_mspaf()] does; a day lacking a required
#' co-analyte leaves that analyte out of that point on that day.
#'
#' **Detection-limit rules and the gap-filler**, per target, analyte and day:
#' where the target is covered and the reference is not, the reference takes
#' its static summary `s = max(geometric mean of its results with
#' below-detection results at LOR/2, median(LOR of its below-detection
#' results) / 2)` (the geometric mean alone when it has none), normalised with
#' the reference's co-analytes that day. An analyte the reference never
#' measured is dropped. When the target is below detection that day (both
#' bracketing samples below detection) its value is set to the reference's,
#' so the analyte adds no hazard. An analyte enters a day only when both
#' sides have a value.
#'
#' **Daily added hazard.** `msPAF` combines the normalised concentrations by
#' concentration addition within mode-of-action group (`TU = sum(C / HC50)`,
#' mean slope) and independent action across groups. With `msT` and `msR` the
#' target's and reference's msPAF on the same analytes,
#' `D = max((msT - msR) / (1 - msR), 0)` (0 when `msR` is 1). Without a
#' reference (`ara = "none"`) below-detection target days enter at 0 and
#' `D = msT`. On this path `ara` defaults to `"contemporaneous"` with a
#' reference and `"none"` without; `"static"`, `"contemporaneous"` without a
#' reference, `impute = TRUE`, a `daily_args$interpolation` other than
#' `"loglinear"` and a `chronic_summary` other than `"arith_mean"` are errors
#' of class `hydroSense_ara_conflict`, raised before anything is computed
#' (`daily_args$interpolation = "loglinear"` with `ara_form = "excess"` is the
#' same error). An unspecified `impute` means `FALSE` here. `hydro`,
#' `reference_args` and the `daily_args` entries other than `interpolation`,
#' `span_cap`, `temperature`, `start`, `end`, `analyte_metadata`, `method`,
#' `guideline_dir` and `conc_units` are ignored.
#'
#' **Time-weighting and coverage.** For a focal date the days of age
#' `a = 0..window` carry weight `exp(-a / tau)` (`tau` 90 and `window` 365
#' days by default). A day is covered when `D` exists. The chronic `value` is
#' `100 * sum(w D) / sum(w)` over the covered days; `coverage` is the covered
#' share of all the weights. `value` is `NA` when `coverage < min_coverage` or
#' nothing is covered. `value_total` (with `also_total = TRUE`) is the same
#' time-weighting (same weights, covered days and `NA` gate) of `msT0`, the
#' target's own msPAF with no reference: its normalised values with
#' below-detection analyte-days at concentration 0, on the same analytes as
#' that day's `D`. Non-detects are therefore not counted at the reference's
#' level in the total.
#'
#' **Contributions.** Per analyte and day the added single-substance PAF is
#' `max(PAF_i(T) - PAF_i(R), 0)`, with `PAF_i(C) = pnorm(log10(C / HC50_i) /
#' sigma_i)`. The day's `D` is apportioned in proportion to it (all zero when
#' no analyte adds any), and the shares are time-weighted with exactly the
#' weights and denominator of `value`, so the contributions of a focal date sum
#' to `value` and are `NA` where it is. A target below detection gives that
#' analyte exactly zero. The `"total"` variant apportions `msT0` by
#' `PAF_i(T0)` (`T0` the target's values with below-detection days at 0, so a
#' below-detection analyte-day contributes 0) the same way and sums to
#' `value_total`. A group's contribution is the sum over its analytes.
#'
#' **Output.** One row per (target site, focal date): `focal_date`, `site_id`,
#' `sample_id` (`"chronic_<focal_date>_<site_id>"`), `analyte` (`"msPAF"`),
#' `value`, [`value_total`], `detected`, `n_samples_in_window` (target sample
#' dates with a toxicant result in the window), `n_imputed_in_window`
#' (covered days in the window that are not sample dates) and `coverage`.
#' Attribute `"daily"` holds `date`, `site_id`, `mspaf` (`100 D`),
#' `mspaf_total` (`100 msT0`), `mspaf_reference` and `n_analytes`; attribute
#' `"sample_days"` has one row per target sample date (`site_id`, `date`,
#' `paired`, `TRUE` when the reference has a toxicant result that day). With
#' `chronic = FALSE` the daily frame is returned, carrying `"sample_days"`
#' (and the daily `"contributions"`, by `date`, when requested).
#'
#' @return When `chronic = TRUE` (default), the chronic msPAF frame, one row per
#'   (focal date × site):
#'   - **point mode** (`ndraws` absent): columns `focal_date`, `site_id`,
#'     `sample_id`, `analyte`, `value`, `detected`, `n_samples_in_window`,
#'     `n_imputed_in_window` (the number of daily rows in the window with no grab
#'     sample, i.e. model-predicted or forward-filled; a single chronic `value`,
#'     no interval columns),
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
#'   For `ara_form = "added_hazard"` see the section *Added-hazard form*.
#'
#' @seealso [fit_imputation_model()], [fit_reference_model()], [mspaf_daily()],
#'   [time_weighted_aggregate()], [expand_focal_dates()],
#'   [interpolate_loglinear()]
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
                           ara = c("static", "contemporaneous", "none"),
                           contributions = FALSE,
                           ara_form = c("excess", "added_hazard"),
                           min_coverage = 0.2) {
  checkmate::assert_data_frame(target)
  checkmate::assert_data_frame(reference, null.ok = TRUE)
  checkmate::assert_list(daily_args)
  ara_form <- match.arg(ara_form)

  # ── Chronic added-hazard form: its own engine, validated before anything ──
  if (ara_form == "added_hazard") {
    checkmate::assert_flag(impute)
    checkmate::assert_flag(chronic)
    checkmate::assert_string(chronic_summary)
    checkmate::assert_flag(also_total)
    checkmate::assert_flag(contributions)
    v <- .ah_validate(
      reference = reference, impute = impute, impute_missing = missing(impute),
      daily_args = daily_args, chronic_summary = chronic_summary,
      ara = ara, ara_missing = missing(ara), min_coverage = min_coverage,
      tau = tau, window = window
    )
    return(.added_hazard_pipeline(
      target = target, reference = v$reference, daily_args = daily_args,
      chronic = chronic, focal_dates = focal_dates, focal_by = focal_by,
      tau = tau, window = window, also_total = also_total,
      contributions = contributions, min_coverage = min_coverage
    ))
  }
  if (identical(daily_args$interpolation, "loglinear")) {
    cli::cli_abort(c(
      "{.code daily_args$interpolation = \"loglinear\"} belongs to \
       {.code ara_form = \"added_hazard\"}.",
      "i" = "Set {.code ara_form = \"added_hazard\"}, or choose another \
             interpolation."
    ), class = "hydroSense_ara_conflict")
  }
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
  checkmate::assert_flag(contributions)
  if (contributions && !is.null(daily_args$ndraws)) {
    cli::cli_abort(c(
      "{.code contributions = TRUE} is supported in point mode only.",
      "i" = "Drop {.code daily_args$ndraws}."
    ), class = "hydroSense_contributions_draws")
  }
  if (contributions && chronic && chronic_summary != "arith_mean") {
    cli::cli_abort(c(
      "{.code contributions = TRUE} needs {.code chronic_summary = \
       \"arith_mean\"}.",
      "i" = "Only the arithmetic time-weighted mean is additive across \
             analytes."
    ), class = "hydroSense_contributions_summary")
  }
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
  if (contributions) daily_call$contributions <- TRUE
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

  if (contributions) {
    attr(out, "contributions") <- .chronic_contributions(
      daily, focal = focal, tau = tau, window = window,
      chronic_summary = chronic_summary, also_total = also_total
    )
    # The daily breakdown is folded into the chronic one; keep the attached
    # daily frame identical to a run without contributions.
    attr(daily, "contributions") <- NULL
  }

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
      value     = daily[[value_col]],
      imputed   = if ("imputed" %in% names(daily)) daily$imputed else FALSE
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


# ── .chronic_contributions() ──────────────────────────────────────────────────

#' Aggregate daily msPAF contributions into chronic contributions
#'
#' Each analyte's daily contribution series (zero on days it did not enter that
#' day's msPAF) is passed through [time_weighted_aggregate()] with the same
#' arguments, focal grid and daily schedule as the chronic msPAF itself, so the
#' chronic contributions sum to the chronic value.
#'
#' @param daily [mspaf_daily()] result carrying a `"contributions"` attribute.
#' @inheritParams .chronic_from_daily
#' @return Tibble `(site_id, focal_date, analyte, group, contribution)`, plus a
#'   leading `variant` column when `also_total`.
#' @keywords internal
.chronic_contributions <- function(daily, focal, tau, window, chronic_summary,
                                   also_total = FALSE) {
  ctb <- attr(daily, "contributions")
  one_variant <- function(x) {
    groups <- dplyr::distinct(x, .data$analyte, .data$group)
    sched <- dplyr::distinct(daily[!is.na(daily$mspaf), , drop = FALSE],
                             .data$site_id, .data$date,
                             .data$imputed)
    sites <- unique(x$site_id)
    long <- dplyr::bind_rows(lapply(sites, function(sid) {
      an <- unique(x$analyte[x$site_id == sid])
      sc <- sched[sched$site_id == sid, , drop = FALSE]
      do.call(rbind, lapply(an, function(a) {
        xa <- x[x$site_id == sid & x$analyte == a, , drop = FALSE]
        v  <- xa$contribution[match(sc$date, xa$date)]
        v[is.na(v)] <- 0
        data.frame(
          sample_id = paste0("daily_", sid, "_", sc$date),
          site_id = sid, datetime = sc$date, analyte = a, detected = TRUE,
          value = v, imputed = sc$imputed, stringsAsFactors = FALSE
        )
      }))
    }))
    agg <- time_weighted_aggregate(
      long, focal_dates = focal, tau = tau, tau_units = "d",
      window = window, window_units = "d", summary = chronic_summary,
      return = "summary"
    )
    agg |>
      dplyr::left_join(groups, by = "analyte") |>
      dplyr::select("site_id", "focal_date", "analyte", "group",
                    contribution = "value") |>
      dplyr::arrange(.data$site_id, .data$focal_date, .data$analyte)
  }
  if (!also_total) return(one_variant(ctb))
  dplyr::bind_rows(
    dplyr::mutate(one_variant(ctb[ctb$variant == "ara", , drop = FALSE]),
                  variant = "ara", .before = 1L),
    dplyr::mutate(one_variant(ctb[ctb$variant == "total", , drop = FALSE]),
                  variant = "total", .before = 1L)
  )
}
