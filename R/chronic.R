#' Time-weighted chronic aggregation for any long-format value column
#'
#' For each (focal date × monitoring feature × analyte) combination, aggregates
#' values from the preceding `window` (in days) using exponential-decay temporal
#' weighting and forward-step duration weighting.  This is the chronic
#' exposure / response predictor used for downstream calibration against
#' biological community state.
#'
#' This function is value-agnostic: pass any long-format data frame with
#' `analyte` and `value` columns and it will compute one time-weighted value
#' per (focal_date × site_id × analyte).  Use it to aggregate raw chemistry
#' (`analyte` = chemical species, `summary = "geom_mean"`) or per-sample
#' msPAF values (`analyte = "msPAF"`, `summary = "arith_mean"`) into a
#' chronic predictor.
#'
#' **Forward-step duration weighting** treats each sample as representing the
#' period from its collection date to the next sample's date (or to
#' `focal_date` for the most recent sample).  This corrects for pulse-biased
#' sampling where storm events are sampled more frequently than base-flow
#' periods, which would otherwise over-weight episodic concentrations.
#'
#' **Exponential-decay temporal weighting** assigns higher weight to recent
#' samples, with a half-life of approximately `tau * log(2)` days (at the
#' default of 90 d).
#'
#' These two components are the only weighting scheme offered, by design:
#' forward-step duration weighting is the minimal correction for irregular /
#' pulse-biased sampling, and exponential decay is the standard memory kernel
#' for a community integrating a fluctuating exposure (one interpretable
#' parameter, `tau`, tied to the target biology's response time).  Richer
#' kernels would add parameters without a defensible way to fit them from
#' routine monitoring data, so they are left out — tune `tau` rather than
#' swapping the kernel.
#'
#' **Combined weight** for sample *i* is
#' `w_i = Δt_i × exp(-(focal_date - midpoint_i) / tau_d)`, where
#' `midpoint_i` is the midpoint of sample *i*'s representative interval.
#'
#' **Choosing `summary`:**
#'
#' - `"geom_mean"`: max-likelihood central tendency for log-normal data.
#'   Use for chemistry concentrations.  Caveat: with highly pulsed exposure
#'   and a non-linear SSD response, this underestimates the time-averaged
#'   risk metric (Jensen's inequality on the upper tail of the SSD).
#' - `"arith_mean"`: weighted arithmetic mean.  Use for bounded indices
#'   like msPAF percentages, or for chemistry when comparison against an
#'   arithmetic-mean compliance trigger is wanted.
#' - `"p90"`: 90th percentile (duration-weighted empirical CDF).  Diagnostic
#'   for "what's the upper-end exposure the community sees most of the
#'   time."
#'
#' For chronic msPAF specifically, the recommended pipeline is
#' `add_mspaf()` on per-sample chemistry, then `time_weighted_aggregate()`
#' on the resulting msPAF rows with `summary = "arith_mean"`.  This
#' computes the time-averaged msPAF, which integrates the toxic response
#' over time (consistent with how biological communities respond to
#' fluctuating exposure) — rather than the msPAF computed at a single
#' time-averaged chemistry, which can substantially under-state risk for
#' pulsed exposures.
#'
#' @param df Long-format data frame.  Required columns:
#'   - `sample_id` (character) — unique sample identifier
#'   - `site_id` (character) — monitoring feature identifier
#'   - `datetime` (Date or POSIXct) — sample collection date/time
#'   - `analyte` (character) — analyte / index name
#'   - `value` (numeric) — concentration or index value
#'   - `detected` (logical) — whether value is a detected observation
#'   Optional: `imputed` (logical) — propagated to `n_imputed_in_window`.
#' @param focal_dates Either a `Date` vector (applied to all features in `df`)
#'   or a data frame with columns `focal_date` (Date) and `site_id`
#'   (character) specifying per-feature focal dates.
#' @param tau Exponential-decay rate parameter.  Numeric or `units` object;
#'   bare numeric requires `tau_units`.  Default `NULL` → 90 days.  The
#'   effective half-life is `tau * log(2)` ≈ 62 days at the default.  Choose
#'   to match the response timescale of the downstream biology (algae:
#'   days–weeks; macroinvertebrates: weeks–months; fish: months–years).
#' @param tau_units Character unit string for `tau` when it is bare numeric,
#'   e.g. `"d"`.  Ignored when `tau` is a `units` object or `NULL`.
#' @param window Look-back window.  Numeric or `units` object; bare numeric
#'   requires `window_units`.  Default `NULL` → 365 days.
#' @param window_units Character unit string for `window` when it is bare
#'   numeric, e.g. `"d"`.  Ignored when `window` is a `units` object or
#'   `NULL`.
#' @param summary Aggregation method: `"geom_mean"` (default), `"arith_mean"`,
#'   `"p90"`.  See *Choosing `summary`* above.
#' @param anchor_outside_window Logical (default `TRUE`).  If `TRUE`, the most
#'   recent sample *before* the window start is included as an anchor to
#'   provide duration weighting at the leading edge of the window. Its
#'   duration is clipped to start at `window_start` so it only contributes
#'   the in-window portion of its interval.
#' @param eps Small positive guard added inside the log for geometric mean
#'   to avoid `log(0)`. Default `1e-9`.
#' @param return Output mode for draw-carrier input (see [summarise_draws()]).
#'   `"summary"` (default) collapses posterior draws to a central estimate plus
#'   a credible interval (`value`, `value_lower`, `value_upper`, `n_draws`);
#'   `"draws"` returns the raw per-draw chronic rows (`draw_id 1..N`).  For
#'   point (non-draw) input both modes return byte-identical output with no
#'   interval columns.  Draws are paired across time by `draw_id` (an
#'   index-pairing approximation that assumes temporal independence; see the
#'   OU/Kalman smoothing roadmap item).
#' @param interval Width of the credible interval when `return = "summary"`.
#'   Default `0.90` (5th/95th percentile bounds).
#' @param central Central-tendency statistic when `return = "summary"`:
#'   `"median"` (default) or `"mean"`.
#'
#' @return A long-format tibble with columns:
#'   - `focal_date` (Date)
#'   - `site_id` (character)
#'   - `sample_id` (character) — synthetic key `"chronic_<focal_date>_<site>"`
#'   - `analyte` (character)
#'   - `value` (numeric) — chronic time-weighted value
#'   - `detected` (logical) — always `TRUE`
#'   - `n_samples_in_window` (integer)
#'   - `n_imputed_in_window` (integer) — count of `imputed == TRUE` samples
#'     contributing; 0 if `imputed` column absent from `df`
#'
#' @examples
#' # Chronic (time-weighted geometric-mean) chemistry for one downstream
#' # analyte at two focal dates, from the bundled demo data.
#' cu <- subset(leachate_demo(), site_id == "downstream" & analyte == "Cu")
#' time_weighted_aggregate(
#'   cu,
#'   focal_dates = as.Date(c("2024-06-01", "2024-12-01")),
#'   tau = 90, tau_units = "d", window = 365, window_units = "d",
#'   summary = "geom_mean"
#' )
#'
#' @export
time_weighted_aggregate <- function(
    df,
    focal_dates,
    tau         = NULL,
    tau_units   = NULL,
    window      = NULL,
    window_units = NULL,
    summary     = c("geom_mean", "arith_mean", "p90"),
    anchor_outside_window = TRUE,
    eps         = 1e-9,
    return      = c("summary", "draws"),
    interval    = 0.90,
    central     = c("median", "mean")
) {
  summary <- match.arg(summary)
  return  <- match.arg(return)
  central <- match.arg(central)

  # ── Input validation ───────────────────────────────────────────────────────
  checkmate::assert_data_frame(df)
  checkmate::assert_names(
    names(df),
    must.include = c("sample_id", "site_id", "datetime",
                     "analyte", "value", "detected")
  )
  tau_days    <- if (is.null(tau))    90  else .resolve_to(tau,    "d", tau_units,    "tau")
  window_days <- if (is.null(window)) 365 else .resolve_to(window, "d", window_units, "window")
  checkmate::assert_number(tau_days,    lower = 0.001)
  checkmate::assert_number(window_days, lower = 1)
  checkmate::assert_flag(anchor_outside_window)
  checkmate::assert_number(eps, lower = 0)

  has_imputed   <- "imputed" %in% names(df)
  is_draws_mode <- "draw_id" %in% names(df) && !all(is.na(df[["draw_id"]]))
  draws         <- .draw_domain(df)

  # Normalise datetime to Date for day-level arithmetic
  df <- dplyr::mutate(df,
    .date = as.Date(.data$datetime)
  )

  # ── Build the focal_dates × site_id grid ──────────────────────────────────
  if (is.data.frame(focal_dates)) {
    checkmate::assert_names(names(focal_dates),
      must.include = c("focal_date", "site_id"))
    focal_grid <- dplyr::mutate(focal_dates,
      focal_date = as.Date(.data$focal_date)
    )
  } else {
    focal_dates <- as.Date(focal_dates)
    features    <- unique(df$site_id)
    focal_grid  <- tidyr::expand_grid(
      focal_date = focal_dates,
      site_id    = features
    )
  }

  # Split the data by site ONCE, up front, rather than re-filtering the full
  # df inside every (focal_date × site_id) iteration below.  `split()` keys
  # the resulting list by site_id, so each iteration is an O(1) list lookup
  # instead of an O(nrow(df)) scan — a large saving when there are many focal
  # dates.  (dplyr::group_split() would also work but returns an unnamed list,
  # forcing a separate key lookup; base split() gives us the named list
  # directly.)
  df_by_site <- split(df, df$site_id)

  # ── Per (focal_date, site_id): compute chronic values ─────────────────────
  results <- purrr::pmap(focal_grid, function(focal_date, site_id) {
    window_start <- focal_date - window_days

    # Samples for this feature (pre-split list lookup)
    feat_df <- df_by_site[[site_id]]

    if (is.null(feat_df) || nrow(feat_df) == 0L) {
      return(NULL)
    }

    # Unique samples (one row per sample_id × date)
    sample_dates <- feat_df |>
      dplyr::distinct(.data$sample_id, .data$.date) |>
      dplyr::arrange(.data$.date)

    in_window <- sample_dates$.date >= window_start & sample_dates$.date <= focal_date

    if (!any(in_window)) {
      return(NULL)
    }

    # Identify anchor: most recent sample strictly before window_start
    anchor_idx <- which(sample_dates$.date < window_start)
    anchor_uid <- if (anchor_outside_window && length(anchor_idx) > 0L) {
      sample_dates$sample_id[max(anchor_idx)]
    } else {
      NA_character_
    }

    # Samples to use: anchor (if any) + in-window samples
    use_uids <- c(
      if (!is.na(anchor_uid)) anchor_uid,
      sample_dates$sample_id[in_window]
    )
    use_dates <- c(
      if (!is.na(anchor_uid)) sample_dates$.date[max(anchor_idx)],
      sample_dates$.date[in_window]
    )

    n_use <- length(use_uids)

    # ── Forward-step interval boundaries ────────────────────────────────────
    # Interval i: [use_dates[i], use_dates[i+1]) for i < n_use
    # Last sample: [use_dates[n_use], focal_date]
    interval_end   <- c(use_dates[-1L], focal_date)
    interval_start <- use_dates

    # Clip anchor interval to start at window_start
    if (!is.na(anchor_uid) && n_use >= 1L) {
      interval_start[1L] <- max(interval_start[1L], window_start)
    }

    delta_t   <- as.numeric(interval_end - interval_start)   # days
    midpoints <- interval_start + delta_t / 2

    # Exponential decay weights at midpoint
    w_time <- exp(-as.numeric(focal_date - midpoints) / tau_days)
    w      <- w_time * pmax(delta_t, 0)

    # Map sample_id → (delta_t, midpoint, weight)
    sample_wt <- tibble::tibble(
      sample_id  = use_uids,
      .weight    = w,
      .n_in_win  = as.integer(in_window[match(use_uids, sample_dates$sample_id)])
    )
    # Anchor is not counted in n_samples_in_window (it's outside the window)
    if (!is.na(anchor_uid)) {
      sample_wt$.n_in_win[sample_wt$sample_id == anchor_uid] <- 0L
    }

    # Filter chemistry to used samples
    chem_use <- feat_df |>
      dplyr::filter(.data$sample_id %in% use_uids) |>
      dplyr::left_join(sample_wt, by = "sample_id")

    # ── Broadcast exact cells and aggregate per (analyte, draw) ──────────
    # Weights and in-window flags are date-derived (draw-agnostic), so
    # broadcasting copies them identically across draws — correct behaviour.
    # In the point case draws=integer(0), broadcast adds draw_id=1L and
    # the group key is effectively just analyte (one group per analyte).
    chem_use <- .broadcast_draws(chem_use, draws)

    synth_uid <- paste0("chronic_", focal_date, "_", site_id)

    analyte_out <- chem_use |>
      dplyr::group_by(.data$analyte, .data$draw_id) |>
      dplyr::summarise(
        value             = .aggregate_weighted(
          .data$value, .data$.weight, summary, eps
        ),
        n_samples_in_window = sum(.data$.n_in_win, na.rm = TRUE),
        n_imputed_in_window = if (has_imputed) {
          sum(.data$imputed & .data$.n_in_win > 0L, na.rm = TRUE)
        } else {
          0L
        },
        .groups = "drop"
      ) |>
      dplyr::mutate(
        focal_date = .env$focal_date,
        site_id    = .env$site_id,
        sample_id  = synth_uid,
        detected   = TRUE,
        .before    = 1L
      )

    analyte_out
  })

  results <- purrr::compact(results)

  if (length(results) == 0L) {
    cli::cli_abort(c(
      "No values could be aggregated for any (focal_date, site_id).",
      "i" = "Check that {.arg focal_dates} and {.arg df$datetime} overlap \\
             within window = {window_days} days, and that \\
             {.arg df$site_id} matches the sites in {.arg focal_dates}."
    ))
  }

  out <- dplyr::bind_rows(results) |>
    dplyr::select(
      "focal_date", "site_id", "sample_id",
      "analyte", "value", "detected",
      "n_samples_in_window", "n_imputed_in_window",
      dplyr::any_of("draw_id")
    )

  if (!is_draws_mode) out <- dplyr::select(out, -dplyr::any_of("draw_id"))
  if (return == "summary") out <- summarise_draws(out, interval, central)
  out
}

# ── expand_focal_dates ────────────────────────────────────────────────────────

#' Generate a sequence of focal dates for chronic msPAF computation
#'
#' A thin convenience wrapper around [base::seq.Date()] for generating the
#' `focal_dates` vector passed to [time_weighted_aggregate()].  The most
#' common use is a daily sequence for time-series analysis.
#'
#' @param start Start date (character `"YYYY-MM-DD"` or `Date`).
#' @param end   End date (character `"YYYY-MM-DD"` or `Date`).
#' @param by    Increment. Passed to [base::seq.Date()]. Common values:
#'   `"day"`, `"week"`, `"month"`. Default `"day"`.
#'
#' @return A `Date` vector from `start` to `end` at the specified increment.
#'
#' @examples
#' # Monthly focal dates across 2024, e.g. to feed time_weighted_aggregate().
#' expand_focal_dates("2024-01-01", "2024-12-31", by = "month")
#'
#' @export
expand_focal_dates <- function(start, end, by = "day") {
  start <- as.Date(start)
  end   <- as.Date(end)
  if (is.na(start)) cli::cli_abort("{.arg start} could not be parsed as a date.")
  if (is.na(end))   cli::cli_abort("{.arg end} could not be parsed as a date.")
  if (end < start)  cli::cli_abort("{.arg end} must be on or after {.arg start}.")
  seq.Date(start, end, by = by)
}


# ── Internal aggregation helper ───────────────────────────────────────────────

.aggregate_weighted <- function(values, weights, method, eps) {
  if (all(is.na(values))) return(NA_real_)
  ok  <- !is.na(values) & !is.na(weights) & weights > 0
  v   <- values[ok]
  w   <- weights[ok]
  if (length(v) == 0L) return(NA_real_)

  switch(method,
    geom_mean = exp(sum(w * log(v + eps)) / sum(w)),
    arith_mean = sum(w * v) / sum(w),
    p90 = {
      # Duration-weighted pseudo-quantile: weight each value by its w, then
      # find the 90th percentile of the weighted empirical distribution.
      ord   <- order(v)
      v_ord <- v[ord]
      w_cum <- cumsum(w[ord]) / sum(w)
      v_ord[which(w_cum >= 0.9)[1L]]
    }
  )
}
