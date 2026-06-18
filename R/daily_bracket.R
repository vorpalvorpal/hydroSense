## ============================================================================
## daily_bracket.R -- summarise the informative/ignorable gap-uncertainty
## bracket into a tidy frame (issue #50)
## ============================================================================
##
## amspaf_daily() in draws mode produces, per (site, day, draw), the AmsPAF under
## two missingness assumptions:
##   * ignorable  (upper) -- the residual is drawn from the simulation smoother;
##                           its variance balloons across gaps (MAR/today).
##   * informative (lower) -- the residual is frozen at its posterior mean on
##                           in-gap days (the fully-informative extreme: gaps
##                           perfectly predictable). Rubin (1976).
## The two are nested and coincide at observation days. .summarise_bracket()
## collapses the raw draws to per-day central tendency + credible bounds for the
## requested envelope(s), plus -- in "bracket" mode -- the precautionary
## composite [lo_informative, hi_ignorable] (a *decision* bound, NOT a calibrated
## interval; its coverage exceeds nominal and is undefined).

#' Summarise per-draw bracket AmsPAF into a tidy per-day frame
#'
#' @param draws_df Long per-draw frame with `date`, `site_id`, `draw_id` and the
#'   envelope value columns needed for `gap_uncertainty`: `amspaf_ignorable`
#'   and/or `amspaf_informative`.
#' @param interval Credible-interval width (default 0.9). The lower bound is the
#'   `(1 - interval)/2` quantile, the upper the `1 - (1 - interval)/2` quantile
#'   (type-7, matching the package's other draw summaries).
#' @param central Per-day central tendency: `"median"` (default) or `"mean"` of
#'   the draws (issue #42 -- the centre is a summary of the draws themselves).
#' @param gap_uncertainty `"bracket"` (both envelopes + precautionary),
#'   `"ignorable"` (ignorable columns only), or `"informative"` (informative
#'   columns only).
#' @return A tibble keyed by (`date`, `site_id`) with the envelope columns for
#'   the requested mode: `median_*`, `lo_*`, `hi_*` per envelope, and
#'   `precautionary_lo`/`precautionary_hi` in `"bracket"` mode.
#' @keywords internal
.summarise_bracket <- function(draws_df, interval = 0.9, central = "median",
                               gap_uncertainty = "bracket") {
  lo_p <- (1 - interval) / 2
  hi_p <- 1 - lo_p
  centre_fun <- if (central == "mean") base::mean else stats::median

  want_ig <- gap_uncertainty %in% c("bracket", "ignorable")
  want_inf <- gap_uncertainty %in% c("bracket", "informative")

  ## Summarise one envelope's value column into median_/lo_/hi_<suffix>.
  summarise_one <- function(value_col, suffix) {
    draws_df |>
      dplyr::group_by(.data$date, .data$site_id) |>
      dplyr::summarise(
        "median_{suffix}" := centre_fun(.data[[value_col]]),
        "lo_{suffix}" := stats::quantile(.data[[value_col]], lo_p, names = FALSE),
        "hi_{suffix}" := stats::quantile(.data[[value_col]], hi_p, names = FALSE),
        .groups = "drop"
      )
  }

  parts <- list()
  if (want_ig) parts$ig <- summarise_one("amspaf_ignorable", "ignorable")
  if (want_inf) parts$inf <- summarise_one("amspaf_informative", "informative")

  out <- Reduce(
    function(a, b) dplyr::full_join(a, b, by = c("date", "site_id")),
    parts
  )

  ## Precautionary composite: only defined when both envelopes are present.
  if (gap_uncertainty == "bracket") {
    out$precautionary_lo <- out$lo_informative
    out$precautionary_hi <- out$hi_ignorable
  }
  out
}

#' Column names of the bracket envelope summary for a given mode
#'
#' Single source of truth for the per-mode summary schema, shared by
#' [.summarise_bracket()]'s consumers and [.empty_daily_result()].
#' @keywords internal
.bracket_summary_cols <- function(gap_uncertainty = "bracket") {
  ig <- c("median_ignorable", "lo_ignorable", "hi_ignorable")
  inf <- c("median_informative", "lo_informative", "hi_informative")
  switch(gap_uncertainty,
    ignorable   = ig,
    informative = inf,
    bracket     = c(inf, ig, "precautionary_lo", "precautionary_hi")
  )
}

#' Column name of the bracket per-draw value(s) for a given mode
#' @keywords internal
.bracket_draw_cols <- function(gap_uncertainty = "bracket") {
  switch(gap_uncertainty,
    ignorable   = "amspaf_ignorable",
    informative = "amspaf_informative",
    bracket     = c("amspaf_ignorable", "amspaf_informative")
  )
}
