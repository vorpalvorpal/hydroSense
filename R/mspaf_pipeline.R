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
#' @param target Long-format target-site chemistry (same schema as
#'   [mspaf_daily()]'s `df`): `sample_id`, `site_id`, `datetime`, `analyte`,
#'   `value`, `detected`.
#' @param reference Long-format reference-site chemistry passed to
#'   [fit_reference_model()].
#' @param hydro Optional hydrology series forwarded to [fit_reference_model()].
#' @param impute Logical. Fit and apply a Bayesian imputation model? Default
#'   `TRUE`. When `FALSE`, no imputation model is fitted and `NULL` is passed
#'   downstream (identical to chaining the functions by hand without one).
#' @param impute_on Which chemistry to fit the imputation model on:
#'   `"reference"` (default, matching the documented design) or `"target"`.
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
#'
#' @return When `chronic = TRUE` (default), the chronic msPAF frame, one row per
#'   (focal date × site):
#'   - **point mode** (`ndraws` absent): columns `focal_date`, `site_id`,
#'     `sample_id`, `analyte`, `value`, `detected`, `n_samples_in_window`,
#'     `n_imputed_in_window` (a single chronic `value`, no interval columns).
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
                           reference,
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
                           chronic_summary = "arith_mean") {
  checkmate::assert_data_frame(target)
  checkmate::assert_data_frame(reference)
  checkmate::assert_flag(impute)
  impute_on <- match.arg(impute_on)
  checkmate::assert_list(reference_args)
  checkmate::assert_list(daily_args)
  checkmate::assert_flag(chronic)
  checkmate::assert_string(chronic_summary)

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
  ref_model <- do.call(
    fit_reference_model,
    c(list(reference, hydro = hydro, imputation_model = imodel), reference_args)
  )

  # ── 3. Daily msPAF on the target ────────────────────────────────────────────
  daily_call <- utils::modifyList(list(interpolation = "model"), daily_args)
  draws_mode <- !is.null(daily_call$ndraws)
  # When chronic && draws_mode, force per-draw rows so they can be propagated
  # into the chronic summary; otherwise honour daily_call$return.
  if (chronic && draws_mode) daily_call$return <- "draws"
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
    central         = central
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
#' is propagated per draw and collapsed with [.summarise_bracket()].
#'
#' @keywords internal
.chronic_from_daily <- function(daily, draws_mode, gap_uncertainty, focal,
                                tau, window, chronic_summary, interval,
                                central) {
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
    long <- to_long("mspaf")
    return(time_weighted_aggregate(
      long,
      focal_dates  = focal,
      tau          = tau,
      tau_units    = "d",
      window       = window,
      window_units = "d",
      summary      = chronic_summary,
      return       = "summary"
    ))
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
