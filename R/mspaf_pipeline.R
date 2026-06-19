# ── mspaf_pipeline() ──────────────────────────────────────────────────────────

#' Run the full daily msPAF pipeline end to end
#'
#' Orchestrates the three steps a daily msPAF analysis usually chains by hand:
#' optionally fit a Bayesian imputation model ([fit_imputation_model()]), fit
#' the reference/impact model ([fit_reference_model()]), then compute the daily
#' multi-substance PAF ([mspaf_daily()]) for the target site. The fitted
#' imputation model is threaded into **both** downstream steps so that
#' below-detection and entirely-absent analytes are filled before the toxicity
#' calculation.
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
#'   `brms::brm()` via [fit_imputation_model()]) so the default-on path is
#'   reproducible.
#' @param reference_args,daily_args Named lists of additional arguments
#'   forwarded to [fit_reference_model()] and [mspaf_daily()] respectively. Do
#'   not include arguments the orchestrator sets itself (`reference`/`hydro`/
#'   `imputation_model` for the reference fit; `df`/`reference_model`/
#'   `imputation_model` for the daily call). `daily_args` may override the
#'   `interpolation = "model"` default.
#'
#' @return The [mspaf_daily()] result, with the fitted models attached as
#'   attributes `"reference_model"` and `"imputation_model"` (the latter `NULL`
#'   when `impute = FALSE` or the fit was skipped). Note that most \pkg{dplyr}
#'   verbs drop attributes, so read them off the returned frame directly.
#'
#' @seealso [fit_imputation_model()], [fit_reference_model()], [mspaf_daily()]
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
                           daily_args = list()) {
  checkmate::assert_data_frame(target)
  checkmate::assert_data_frame(reference)
  checkmate::assert_flag(impute)
  impute_on <- match.arg(impute_on)
  checkmate::assert_list(reference_args)
  checkmate::assert_list(daily_args)

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
  out <- do.call(
    mspaf_daily,
    c(
      list(target, reference_model = ref_model, imputation_model = imodel),
      daily_call
    )
  )

  attr(out, "reference_model") <- ref_model
  attr(out, "imputation_model") <- imodel
  out
}
