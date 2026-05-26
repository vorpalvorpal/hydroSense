#' Compute time-weighted chronic chemistry for given focal dates
#'
#' For each (focal date × monitoring feature × analyte) combination, aggregates
#' chemistry samples from the preceding `window_days` using exponential-decay
#' temporal weighting and forward-step duration weighting. The result is the
#' chronic exposure predictor suitable for regression against macroinvertebrate
#' community data (e.g. Baselga nestedness component of beta-diversity).
#'
#' **Forward-step duration weighting** treats each sample as representing the
#' period from its collection date to the next sample's date (or to
#' `focal_date` for the most recent sample). This corrects for pulse-biased
#' sampling where storm events are sampled more frequently than base-flow
#' periods, which would otherwise over-weight episodic concentrations.
#'
#' **Exponential-decay temporal weighting** assigns higher weight to recent
#' samples, with a half-life of approximately `tau_days * log(2)` days.
#'
#' **Combined weight** for sample *i* is
#' `w_i = Δt_i × exp(-(focal_date - midpoint_i) / tau_days)`, where
#' `midpoint_i` is the midpoint of sample *i*'s representative interval.
#'
#' @param df Long-format chemistry data frame, typically the output of
#'   `impute_chemistry()`. Required columns:
#'   - `sample_id` (character) — unique sample identifier
#'   - `site_id` (character) — monitoring feature identifier
#'   - `datetime` (Date or POSIXct) — sample collection date/time
#'   - `analyte` (character) — analyte name
#'   - `value` (numeric) — concentration
#'   - `detected` (logical) — whether value is a detected concentration
#'   Optional: `imputed` (logical) — propagated to `n_imputed_in_window`.
#' @param focal_dates Either a `Date` vector (applied to all features in `df`)
#'   or a data frame with columns `focal_date` (Date) and `site_id`
#'   (character) specifying per-feature focal dates.
#' @param tau_days Exponential decay half-life parameter in days. Default 90.
#' @param window_days Look-back window in days. Default 365.
#' @param summary Aggregation method. `"geom_mean"` (default) computes the
#'   time-weighted geometric mean; `"arith_mean"` the arithmetic mean;
#'   `"p90"` the 90th percentile (ignores temporal weights, uses duration
#'   weights only as selection, not for p90 calculation — mostly for
#'   diagnostic use).
#' @param anchor_outside_window Logical (default `TRUE`). If `TRUE`, the most
#'   recent sample *before* the window start is included as an anchor to
#'   provide duration weighting at the leading edge of the window. Its
#'   duration is clipped to start at `window_start` so it only contributes
#'   the in-window portion of its interval.
#' @param eps Small positive guard added inside the log for geometric mean
#'   to avoid `log(0)`. Default `1e-9`.
#'
#' @return A long-format tibble with columns:
#'   - `focal_date` (Date)
#'   - `site_id` (character)
#'   - `sample_id` (character) — synthetic key `"chronic_<focal_date>_<site>"`
#'   - `analyte` (character)
#'   - `value` (numeric) — chronic aggregated concentration
#'   - `detected` (logical) — always `TRUE`
#'   - `n_samples_in_window` (integer)
#'   - `n_imputed_in_window` (integer) — count of `imputed == TRUE` samples
#'     contributing; 0 if `imputed` column absent from `df`
#'
#' @examples
#' \dontrun{
#' bio_dates <- as.Date(c("2024-04-01", "2025-04-01", "2026-04-01"))
#' chr_chem <- compute_chronic_chemistry(
#'   imp,
#'   focal_dates = bio_dates,
#'   tau_days    = 90,
#'   window_days = 365
#' )
#' }
#'
#' @export
compute_chronic_chemistry <- function(
    df,
    focal_dates,
    tau_days    = 90,
    window_days = 365,
    summary     = c("geom_mean", "arith_mean", "p90"),
    anchor_outside_window = TRUE,
    eps         = 1e-9
) {
  summary <- match.arg(summary)

  # ── Input validation ───────────────────────────────────────────────────────
  checkmate::assert_data_frame(df)
  checkmate::assert_names(
    names(df),
    must.include = c("sample_id", "site_id", "datetime",
                     "analyte", "value", "detected")
  )
  checkmate::assert_number(tau_days, lower = 0.001)
  checkmate::assert_number(window_days, lower = 1)
  checkmate::assert_flag(anchor_outside_window)
  checkmate::assert_number(eps, lower = 0)

  has_imputed <- "imputed" %in% names(df)

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

  # ── Per (focal_date, site_id): compute chronic values ─────────────────────
  results <- purrr::pmap(focal_grid, function(focal_date, site_id) {
    window_start <- focal_date - window_days

    # Samples for this feature
    feat_df <- dplyr::filter(df,
      .data$site_id == .env$site_id
    )

    if (nrow(feat_df) == 0L) {
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
    # anchor is NOT counted in n_samples_in_window
    sample_wt$.n_in_win[sample_wt$sample_id == (if (!is.na(anchor_uid)) anchor_uid else "__none__")] <- 0L

    # Filter chemistry to used samples
    chem_use <- feat_df |>
      dplyr::filter(.data$sample_id %in% use_uids) |>
      dplyr::left_join(sample_wt, by = "sample_id")

    # ── Aggregate per analyte ─────────────────────────────────────────────
    synth_uid <- paste0("chronic_", focal_date, "_", site_id)

    analyte_out <- chem_use |>
      dplyr::group_by(.data$analyte) |>
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
      "No chronic chemistry could be computed.",
      "i" = "Check that {.arg focal_dates} and {.arg df$datetime} overlap \\
             within {.arg window_days} = {window_days} days."
    ))
  }

  dplyr::bind_rows(results) |>
    dplyr::select(
      "focal_date", "site_id", "sample_id",
      "analyte", "value", "detected",
      "n_samples_in_window", "n_imputed_in_window"
    )
}

# ── expand_focal_dates ────────────────────────────────────────────────────────

#' Generate a sequence of focal dates for chronic AmsPAF computation
#'
#' A thin convenience wrapper around [base::seq.Date()] for generating the
#' `focal_dates` vector passed to [compute_chronic_chemistry()]. The most
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
#' \dontrun{
#' # Daily sequence for 2024–2025
#' focal_dates <- expand_focal_dates("2024-01-01", "2025-12-31", by = "day")
#'
#' chr_chem <- compute_chronic_chemistry(imp, focal_dates = focal_dates)
#' }
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
