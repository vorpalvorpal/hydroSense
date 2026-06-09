## ============================================================================
## amspaf_daily -- continuous daily-resolved AmsPAF time series
## ============================================================================
##
## The core idea: grab chemistry is sparse (bi-monthly, weekly, etc.); daily
## AmsPAF requires daily chemistry.  This function interpolates each analyte
## onto a fine date grid, constructs synthetic "one-per-day" samples, and
## runs the existing add_amspaf() engine on them.
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
## add_amspaf(), which would double-integrate the ARA reference.  We want
## .resolve_ref_norm_instant() -- pointwise matching per daily sample.

## ============================================================================
## amspaf_daily
## ============================================================================

#' Continuous daily AmsPAF time series from interpolated grab chemistry
#'
#' Interpolates per-analyte grab chemistry onto a daily date grid and computes
#' AmsPAF for every day within the requested date range.  The result is a tidy
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
#' interpolation purposes, matching the treatment in [add_amspaf()].
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
#'   Chemistry for multiple sites may be stacked; interpolation and AmsPAF are
#'   computed per site.
#' @param temperature Optional daily water temperature data frame for days
#'   without a grab temperature measurement. Required columns: `datetime` (Date
#'   or POSIXct) and `value` (temperature in \eqn{{}^\circ}C). The output of
#'   [estimate_water_temp()] is accepted directly (extra columns are ignored).
#'   When both this argument and grab-sample temperature rows in `df` are
#'   present for the same day, the grab measurement takes priority.
#'   `NULL` (default) means temperature must come from `df` rows alone.
#' @param reference Background reference chemistry for ARA adjustment. Accepts
#'   the same four forms as [add_amspaf()]. Controls **only** whether background
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
#'   the bundled metadata. Passed to [add_amspaf()].
#' @param method SSD method. One of `"multi"` (default, model-averaged) or
#'   `"anzecc"`. Passed to [add_amspaf()].
#' @param guideline_dir Path to the ANZG guideline data folder. Falls back to
#'   `getOption("leachatetools.guideline_dir")`.
#' @param min_analytes Minimum number of SSD-eligible analytes per day for
#'   AmsPAF to be computed. Default `3L`.
#' @param conc_units Character. Concentration units for SSD-eligible rows when
#'   `df` lacks a `units.analyte` column. Passed to [add_amspaf()].
#' @param require_temperature Logical (default `TRUE`). When `TRUE`, any daily
#'   sample with `NH3-N` must also carry a `temperature` value. Passed to
#'   [add_amspaf()]. Set `FALSE` for datasets without ammonia.
#'
#' @return A tibble with one row per (site \eqn{\times} day) for days with
#'   sufficient analyte coverage:
#'   \describe{
#'     \item{`date`}{Date of this daily estimate.}
#'     \item{`site_id`}{Site identifier.}
#'     \item{`amspaf`}{Daily AmsPAF (percentage, 0--100+).}
#'     \item{`n_analytes_used`}{SSD-eligible analytes contributing to AmsPAF.}
#'     \item{`dominant_analyte`}{Analyte with the highest individual PAF.}
#'     \item{`max_paf`}{PAF of the dominant analyte (proportion 0--1).}
#'     \item{`n_measured_analytes`}{SSD-eligible analytes with a direct grab
#'       sample on this day (not interpolated).}
#'     \item{`days_since_last_sample`}{Days since the most recent grab sample
#'       for any SSD-eligible analyte. Helps identify heavily interpolated
#'       regions.}
#'     \item{`analyte_pafs`}{List column of per-analyte diagnostic tibbles
#'       (same structure as from [add_amspaf()]).}
#'   }
#'   An `"ara_summary"` attribute is attached; retrieve it with [ara_summary()].
#'
#' @seealso [add_amspaf()], [time_weighted_aggregate()],
#'   [estimate_water_temp()], [get_silo_air_temp()]
#'
#' @examples
#' \donttest{
#' demo <- leachate_demo()
#' ds  <- subset(demo, site_id == "downstream")
#' out <- amspaf_daily(ds, require_temperature = FALSE)
#' head(out[, c("date", "site_id", "amspaf", "n_measured_analytes",
#'              "days_since_last_sample")])
#' }
#' @export
amspaf_daily <- function(
    df,
    temperature          = NULL,
    reference            = NULL,
    reference_model      = NULL,
    imputation_model     = NULL,
    start                = NULL,
    end                  = NULL,
    by                   = "day",
    interpolation        = c("forward_fill", "linear", "model"),
    leading_edge         = c("drop", "backfill"),
    analyte_metadata     = NULL,
    method               = c("multi", "anzecc"),
    guideline_dir        = getOption("leachatetools.guideline_dir"),
    min_analytes         = 3L,
    conc_units           = NULL,
    require_temperature  = TRUE
) {
  ## --- Validate inputs -------------------------------------------------------
  checkmate::assert_data_frame(df)
  checkmate::assert_names(names(df),
    must.include = c("sample_id", "site_id", "datetime", "analyte", "value"))
  interpolation <- match.arg(interpolation)
  leading_edge  <- match.arg(leading_edge)
  method        <- match.arg(method)
  checkmate::assert_flag(require_temperature)
  checkmate::assert_int(min_analytes, lower = 1L)
  checkmate::assert_string(by, min.chars = 1L)

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

  if (!is.null(temperature)) {
    checkmate::assert_data_frame(temperature)
    checkmate::assert_names(names(temperature),
      must.include = c("datetime", "value"))
  }

  ## Draws not yet supported: interpolating draw-bearing chemistry onto a
  ## daily grid requires the temporal-correlation design from issue #16.
  ## Use add_amspaf() on the draw-carrier frame then time_weighted_aggregate()
  ## to propagate draws through the chronic pipeline instead.
  if ("draw_id" %in% names(df) && !all(is.na(df[["draw_id"]]))) {
    cli::cli_abort(c(
      "{.fn amspaf_daily} does not yet support draws-mode input.",
      "i" = "To propagate imputation uncertainty through the chronic pipeline, \\
             call {.fn add_amspaf} on the draw-carrier frame, then \\
             {.fn time_weighted_aggregate} to aggregate draws over time."
    ))
  }

  ## --- Normalise datetime to Date; ensure detected column -------------------
  df <- dplyr::mutate(df, datetime = as.Date(.data$datetime))
  if (!"detected" %in% names(df)) {
    df <- dplyr::mutate(df, detected = TRUE)
  }

  ## --- Determine date range --------------------------------------------------
  if (is.null(start)) start <- min(df$datetime, na.rm = TRUE)
  if (is.null(end))   end   <- max(df$datetime, na.rm = TRUE)
  start     <- as.Date(start)
  end       <- as.Date(end)
  all_dates <- seq(start, end, by = by)

  ## --- Load analyte metadata once (needed for log-space decisions) -----------
  meta <- .load_analyte_metadata(analyte_metadata)

  tox_analytes <- meta$analyte[!is.na(meta$ssd_available) &
                                meta$ssd_available == TRUE]

  ## --- Per-site processing ---------------------------------------------------
  sites <- unique(df$site_id)

  site_results <- lapply(sites, function(site) {
    site_rows <- dplyr::filter(df, .data$site_id == .env$site)

    ## Step 1: Interpolate each analyte onto the daily grid.  For the "model"
    ## path, co-analytes are forward-filled here and toxicants are overwritten
    ## in Step 1b by the fitted target model.
    base_interp <- if (interpolation == "model") "forward_fill" else interpolation
    daily_long <- .build_daily_chem(
      site_rows     = site_rows,
      dates         = all_dates,
      interpolation = base_interp,
      leading_edge  = leading_edge,
      tox_analytes  = tox_analytes
    )

    if (nrow(daily_long) == 0L) return(NULL)

    ## Step 2: Fill temperature from external series on non-grab days.
    if (!is.null(temperature)) {
      daily_long <- .fill_external_temperature(daily_long, temperature)
    }

    ## Step 1b: model interpolation of toxicants (season-blind impact model).
    impact_tiers <- NULL
    if (interpolation == "model") {
      daily_long   <- .daily_tox_from_model(
        daily_long       = daily_long,
        site_rows        = site_rows,
        reference_model  = reference_model,
        imputation_model = imputation_model,
        conc_units       = conc_units,
        meta             = meta,
        tox_analytes     = tox_analytes
      )
      impact_tiers <- attr(daily_long, "impact_tiers")
    }

    ## Step 3: Compute per-day diagnostics from SSD-eligible rows.
    diag <- .compute_daily_diag(daily_long, tox_analytes, site)

    ## Step 4: Build synthetic long-format samples (no focal_date!).
    synth <- .build_synthetic_samples(daily_long, site)

    list(synth = synth, diag = diag, tiers = impact_tiers, site = site)
  })

  site_results <- Filter(Negate(is.null), site_results)
  if (length(site_results) == 0L) {
    cli::cli_warn("No daily chemistry could be built. Returning empty tibble.")
    return(.empty_daily_result())
  }

  all_synth <- dplyr::bind_rows(lapply(site_results, `[[`, "synth"))
  all_diag  <- dplyr::bind_rows(lapply(site_results, `[[`, "diag"))

  ## Build sample_id -> date lookup before passing to add_amspaf() (add_amspaf
  ## may rearrange rows but sample_id is stable throughout).
  id_date_map <- dplyr::distinct(
    dplyr::select(all_synth, "sample_id", "site_id", ".date")
  )

  ## Remove .date from the df passed to add_amspaf() so it stays unaware of it.
  all_synth_clean <- dplyr::select(all_synth, -".date")

  ## --- Run add_amspaf on the daily synthetic samples -------------------------
  amspaf_out <- add_amspaf(
    df                  = all_synth_clean,
    reference           = reference,
    analyte_metadata    = analyte_metadata,
    method              = method,
    guideline_dir       = guideline_dir,
    min_analytes        = min_analytes,
    conc_units          = conc_units,
    require_temperature = require_temperature
  )

  ara_summ <- attr(amspaf_out, "ara_summary")

  ## Attach the target model's per-analyte impact tier ("model" / "bridge") to
  ## the ARA diagnostics. ara_summary() is keyed by (sample_id, analyte); the
  ## synthetic sample_id maps to (site_id, date) via id_date_map, and the impact
  ## tier is per (site_id, analyte). NA for non-modelled (forward-filled)
  ## analytes; absent entirely for non-"model" interpolation.
  all_tiers <- dplyr::bind_rows(lapply(site_results, function(z) {
    if (is.null(z$tiers) || nrow(z$tiers) == 0L) return(NULL)
    dplyr::mutate(z$tiers, site_id = z$site)
  }))
  if (!is.null(ara_summ) && nrow(all_tiers) > 0L) {
    site_lookup <- dplyr::distinct(id_date_map[, c("sample_id", "site_id")])
    ara_summ <- ara_summ |>
      dplyr::left_join(site_lookup, by = "sample_id") |>
      dplyr::left_join(all_tiers, by = c("site_id", "analyte"))
  }

  ## --- Extract and annotate AmsPAF rows -------------------------------------
  amspaf_rows <- dplyr::filter(amspaf_out, .data$analyte == "AmsPAF")

  if (nrow(amspaf_rows) == 0L) {
    cli::cli_warn(
      "No daily AmsPAF rows produced. \\
       Check {.arg min_analytes} ({min_analytes}) and data coverage."
    )
    result <- .empty_daily_result()
    attr(result, "ara_summary") <- ara_summ
    return(result)
  }

  result <- amspaf_rows |>
    dplyr::left_join(id_date_map, by = c("sample_id", "site_id")) |>
    dplyr::rename(date = ".date", amspaf = "value") |>
    dplyr::left_join(all_diag, by = c("date", "site_id")) |>
    dplyr::select(
      "date", "site_id", "amspaf",
      "n_analytes_used", "dominant_analyte", "max_paf",
      "n_measured_analytes", "days_since_last_sample",
      "analyte_pafs"
    ) |>
    dplyr::arrange(.data$site_id, .data$date)

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
  analytes     <- unique(site_rows$analyte)
  passthru_cols <- intersect(names(site_rows),
                              c("units.analyte", "valence.analyte",
                                "atomic_mass.analyte"))

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

    if (nrow(interp) == 0L) return(NULL)

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
  uniq_mask  <- !duplicated(obs_dates)
  obs_d      <- obs_dates[uniq_mask]
  obs_v      <- obs_values[uniq_mask]
  obs_det    <- obs_detected[uniq_mask]

  ## Sort ascending.
  ord     <- order(obs_d)
  obs_d   <- obs_d[ord]
  obs_v   <- obs_v[ord]
  obs_det <- obs_det[ord]

  ## Work in log-space for toxicants (floor at machine epsilon to avoid log(0)).
  if (log_space) {
    obs_v_work <- log(pmax(obs_v, .Machine$double.eps))
  } else {
    obs_v_work <- obs_v
  }

  n   <- length(target_dates)
  out_val  <- numeric(n)
  out_det  <- logical(n)
  out_meas <- logical(n)
  out_keep <- logical(n)

  ## For each target date find its bracketing observations.
  ## findInterval returns 0 when d < obs_d[1], else the index of the largest
  ## obs_d <= d.
  prev_idx <- findInterval(as.numeric(target_dates), as.numeric(obs_d))

  for (i in seq_len(n)) {
    d   <- target_dates[i]
    pid <- prev_idx[i]

    exact_match <- pid > 0L && obs_d[pid] == d

    if (exact_match) {
      raw        <- obs_v_work[pid]
      out_val[i]  <- if (log_space) exp(raw) else raw
      out_det[i]  <- obs_det[pid]
      out_meas[i] <- TRUE
      out_keep[i] <- TRUE
      next
    }

    if (pid == 0L) {
      ## Before the first observation.
      if (leading_edge == "backfill") {
        raw        <- obs_v_work[1L]
        out_val[i]  <- if (log_space) exp(raw) else raw
        out_det[i]  <- obs_det[1L]
        out_meas[i] <- FALSE
        out_keep[i] <- TRUE
      }
      ## else: drop (out_keep stays FALSE)
      next
    }

    ## We have a previous observation but no exact match.
    nid <- pid + 1L  ## index of the next observation (may be out-of-bounds)

    if (interpolation == "forward_fill" || nid > length(obs_d)) {
      ## Forward-fill, or past the last observation: carry forward.
      raw        <- obs_v_work[pid]
      out_val[i]  <- if (log_space) exp(raw) else raw
      out_det[i]  <- obs_det[pid]
      out_meas[i] <- FALSE
      out_keep[i] <- TRUE
    } else {
      ## Linear / log-linear interpolation between pid and nid.
      frac <- as.numeric(d - obs_d[pid]) /
              as.numeric(obs_d[nid] - obs_d[pid])
      raw  <- obs_v_work[pid] + frac * (obs_v_work[nid] - obs_v_work[pid])
      out_val[i]  <- if (log_space) exp(raw) else raw
      ## A gap between detected/BDL passes detected only if both anchors are.
      out_det[i]  <- obs_det[pid] && obs_det[nid]
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
  new_rows         <- ext[no_existing_grab, , drop = FALSE]

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
#' Operates on SSD-eligible rows only (toxicants drive the AmsPAF; co-analyte
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
      date                  = as.Date(character()),
      site_id               = character(),
      n_measured_analytes   = integer(),
      days_since_last_sample = integer()
    ))
  }

  ## All dates on which any SSD-eligible analyte had a grab measurement.
  grab_dates <- sort(unique(tox_rows$.date[tox_rows$.measured]))

  purrr::map_dfr(daily_dates, function(d) {
    n_meas <- sum(tox_rows$.date == d & tox_rows$.measured, na.rm = TRUE)
    prev   <- grab_dates[grab_dates <= d]
    days_since <- if (length(prev) > 0L) {
      as.integer(d - max(prev))
    } else {
      NA_integer_
    }
    tibble::tibble(
      date                  = d,
      site_id               = site,
      n_measured_analytes   = as.integer(n_meas),
      days_since_last_sample = days_since
    )
  })
}


#' Build synthetic long-format daily samples from interpolated chemistry
#'
#' Assigns `sample_id = "daily_{YYYY-MM-DD}_{site}"` per day.  Keeps `.date`
#' as a column (the caller extracts it before passing to [add_amspaf()]).
#' No `focal_date` column is added -- this is deliberate so [add_amspaf()]
#' uses the instant (pointwise) ARA path, not the chronic integrated path.
#'
#' @param daily_long Output of `.build_daily_chem()` (after temperature fill).
#' @param site Site identifier string.
#' @return Long-format tibble ready for [add_amspaf()] (after removing `.date`).
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
      sample_id = paste0("daily_", format(.data$.date, "%Y-%m-%d"),
                         "_", .env$site),
      site_id   = .env$site,
      datetime  = .data$.date
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
                               conc_units, meta, tox_analytes, daily_long) {
  qdates <- sort(unique(daily_long$.date))

  tm <- tryCatch(
    fit_target_model(
      target           = site_rows,
      reference_model  = reference_model,
      imputation_model = imputation_model,
      conc_units       = conc_units,
      analyte_metadata = meta
    ),
    error = function(e) {
      cli::cli_warn(c(
        "Target model fit failed; falling back to forward-fill for toxicants.",
        "x" = conditionMessage(e)
      ))
      NULL
    }
  )
  if (is.null(tm) || length(tm$models) == 0L) return(NULL)

  modelled <- names(tm$models)

  ## Per-date co-analyte lookup (static across draws; Chunk D may override).
  co <- daily_long[!daily_long$analyte %in% tox_analytes &
                     (is.na(daily_long$detected) | daily_long$detected), , drop = FALSE]
  co_split <- split(
    data.frame(analyte = co$analyte, value = co$value, stringsAsFactors = FALSE),
    as.character(co$.date)
  )

  ## WQ layer data for analytes with a fitted WQ→metal response.
  wq_long <- if (!is.null(tm$pca)) {
    tibble::tibble(sample_id = as.character(co$.date),
                   analyte = co$analyte, value = co$value)
  } else NULL

  ## Normalisation formula lookup (one parsed formula per modelled analyte).
  meta_norm <- meta |>
    dplyr::select("analyte", "normalisation_formula", "coanalytes_required")
  fac_lookup <- stats::setNames(
    lapply(modelled, function(a) {
      row <- meta_norm[meta_norm$analyte == a, , drop = FALSE]
      list(
        parsed = if (nrow(row)) .parse_normalisation_formula(
                   row$normalisation_formula %||% "") else NULL
      )
    }), modelled
  )

  ## Measured (grab) dates — used to set .measured flags on synthetic rows.
  sr_mod <- site_rows[site_rows$analyte %in% modelled &
                        (is.na(site_rows$detected) | site_rows$detected), ,
                      drop = FALSE]
  measured_key <- paste(sr_mod$analyte, as.Date(sr_mod$datetime))

  ## OU bridge: precompute Cholesky factors once per analyte over qdates.
  ## WQ-tier analytes have their own residual d_anchors; others use anchors$S.
  ## Also precompute C_norm_obs at anchor dates for S6 measurement-error scaling.
  ou <- stats::setNames(
    lapply(modelled, function(nm) {
      m      <- tm$models[[nm]]
      has_wq <- !is.null(m$wq_fit) && !is.null(m$d_anchors) &&
                  nrow(m$d_anchors) >= 2L
      anch   <- if (has_wq) m$d_anchors else m$anchors

      if (is.null(anch) || nrow(anch) < 2L) {
        return(list(
          params          = list(theta = 0, sigma2 = 0, gamma = 0, degenerate = TRUE),
          factors         = .ou_bridge_factors(as.Date(character()), qdates, 0, 0, 0),
          c_norm_obs_anch = NULL,
          has_wq          = has_wq
        ))
      }

      params  <- .estimate_ou_params(anch$date, anch$S)
      factors <- .ou_bridge_factors(
        anch$date, qdates, params$theta, params$sigma2, params$gamma
      )

      ## C_norm_obs at anchor dates: needed to scale S6 perturbations correctly.
      ## Model/bridge: C_norm = ref_norm + I.  WQ-tier: C_norm = WQ_pred + d.
      c_norm_obs_anch <- tryCatch({
        if (!has_wq) {
          ref_q   <- .resolve_ref_norm_instant(
            tm$reference_model,
            tibble::tibble(sample_id = as.character(anch$date), datetime = anch$date)
          )
          ref_lkp <- stats::setNames(
            ref_q$ref_norm[ref_q$analyte == nm],
            ref_q$sample_id[ref_q$analyte == nm]
          )
          ref_vec <- as.numeric(ref_lkp[as.character(anch$date)])
          ref_vec[is.na(ref_vec)] <- 0
          pmax(anch$I + ref_vec, 0)
        } else if (!is.null(tm$pca) && !is.null(wq_long)) {
          pc_anch <- .compute_pca_scores(wq_long, tm$pca)
          nd_anch <- dplyr::left_join(
            tibble::tibble(sample_id = as.character(anch$date)),
            pc_anch, by = "sample_id"
          )
          wq_pred <- as.numeric(stats::predict(m$wq_fit, newdata = nd_anch))
          pmax(wq_pred + anch$S, 0)
        } else NULL
      }, error = function(e) NULL)

      list(
        params          = params,
        factors         = factors,
        c_norm_obs_anch = c_norm_obs_anch,
        has_wq          = has_wq
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
        idx <- findInterval(qdates, gd)   # 0 = before first grab
        src <- ifelse(idx == 0L, NA_character_, as.character(gd[pmax(idx, 1L)]))
        stats::setNames(src, as.character(qdates))
      }),
      co_analytes_nm
    )
  } else NULL

  list(
    tm           = tm,
    modelled     = modelled,
    qdates       = qdates,
    co_split     = co_split,
    wq_long      = wq_long,
    fac_lookup   = fac_lookup,
    measured_key = measured_key,
    ou           = ou,
    co_grab_map  = co_grab_map
  )
}


#' Per-draw daily toxicant prediction from a pre-fitted scaffold (issue #16)
#'
#' Predicts normalised and raw concentrations for one draw (or for the
#' deterministic point mode when `eps_paths = NULL`).  Always uses the
#' precomputed static scaffolding from `fdm`; accepts an optionally-perturbed
#' target model (`tm_p`) and OU bridge fluctuations (`eps_paths`) for draw
#' mode.  Pass overrides to `co_split` and `wq_long` (from
#' [.perturb_co_split()]) for S7 co-analyte measurement-error draws.
#'
#' @param fdm Output of [.fit_daily_target()].
#' @param tm_p Target model to predict with (default: `fdm$tm`).  For draw
#'   mode, pass a GAM-perturbed copy from [.perturb_target_model()].
#' @param eps_paths Named list of mean-zero OU bridge fluctuation paths, one
#'   numeric vector of length `length(fdm$qdates)` per analyte.  `NULL` (the
#'   default) means no fluctuation — output is byte-identical to the original
#'   point-mode path.
#' @param co_split Per-date co-analyte lookup (default: `fdm$co_split`).
#'   Chunk D supplies a perturbed version for S7 draws.
#' @param wq_long WQ layer data (default: `fdm$wq_long`).
#' @return Tibble of model rows (`.date`, `value`, `detected`, `.measured`,
#'   `analyte`, `units.analyte`) with `attr("impact_tiers")` attached, or
#'   `NULL` if prediction produced no finite rows.
#' @keywords internal
.predict_daily_tox <- function(fdm,
                                tm_p      = fdm$tm,
                                eps_paths = NULL,
                                co_split  = fdm$co_split,
                                wq_long   = fdm$wq_long) {
  pred <- .resolve_target_impact(tm_p, tibble::tibble(date = fdm$qdates),
                                  fdm$modelled, wq = wq_long)
  if (nrow(pred) == 0L) return(NULL)

  ## Post-hoc OU bridge ε injection (S4 / S5).
  ## Add the mean-zero fluctuation ε_d(t) to C_norm for each analyte, then
  ## re-clamp at 0.  The bridge guarantees ε = 0 at anchor dates, so the
  ## deterministic mean (.interp_residual centre line) is unchanged there.
  if (!is.null(eps_paths)) {
    for (nm in fdm$modelled) {
      eps_nm <- eps_paths[[nm]]
      if (!is.null(eps_nm) && length(eps_nm) > 0L) {
        idx_nm   <- which(pred$analyte == nm)
        date_pos <- match(pred$date[idx_nm], fdm$qdates)
        valid    <- !is.na(date_pos)
        pred$C_norm[idx_nm[valid]] <- pmax(
          pred$C_norm[idx_nm[valid]] + eps_nm[date_pos[valid]],
          0
        )
      }
    }
  }

  ## Reconstruct raw µg/L: C_raw = C_norm / normalisation_factor(co-analytes).
  co_vec_for <- function(d) {
    cd <- co_split[[as.character(d)]]
    if (is.null(cd)) return(numeric(0))
    stats::setNames(cd$value, cd$analyte)
  }
  pred$C_raw <- vapply(seq_len(nrow(pred)), function(i) {
    a      <- pred$analyte[i]; d <- pred$date[i]
    parsed <- fdm$fac_lookup[[a]]$parsed
    if (is.null(parsed)) return(pred$C_norm[i])
    factor <- .apply_normalisation(parsed, 1, co_vec_for(d))
    if (is.na(factor) || factor <= 0) return(NA_real_)
    pred$C_norm[i] / factor
  }, numeric(1L))

  pred_ok <- pred[is.finite(pred$C_raw), , drop = FALSE]
  if (nrow(pred_ok) == 0L) return(NULL)

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


## ── Chunk D: measurement-error perturbation helpers (issue #16, S6 + S7) ─────

#' Perturb anchor residual-state values by grab measurement error (S6)
#'
#' For each modelled analyte, draws one lognormal multiplier per anchor grab
#' (geo-mean = 1, CV = `grab_cv[[nm]]`) and shifts the anchor's S value by
#' `ΔS = C_norm_obs × (mult − 1)`, where `C_norm_obs` was precomputed in
#' [.fit_daily_target()].  WQ-tier analytes' `d_anchors` are treated the same
#' way (their residual `d` responds to C_norm perturbation identically).
#'
#' The updated `tm` copy is then passed to [.predict_daily_tox()] for one
#' draw; the `.interp_residual()` bridge will interpolate between the
#' perturbed anchor values, providing at-anchor spread equal to the grab
#' measurement width and mid-gap spread equal to that plus the OU balloon.
#'
#' @param tm A `target_model` (typically already GAM-perturbed).
#' @param fdm Fitted daily scaffold from [.fit_daily_target()]; provides
#'   `c_norm_obs_anch` and `has_wq` flags.
#' @param grab_cv Named numeric vector of CVs per analyte, or a single scalar
#'   applied to all analytes.  Analytes not present in the vector (and without
#'   a scalar default) are left unperturbed.
#' @return A modified copy of `tm` with perturbed anchor S values.
#' @keywords internal
.perturb_anchors_in_model <- function(tm, fdm, grab_cv) {
  if (is.null(grab_cv) || (length(grab_cv) == 1L && is.na(grab_cv))) return(tm)

  get_cv <- function(nm) {
    if (length(grab_cv) == 1L)           return(as.numeric(grab_cv))
    if (nm %in% names(grab_cv))          return(as.numeric(grab_cv[[nm]]))
    return(NA_real_)
  }

  for (nm in fdm$modelled) {
    cv <- get_cv(nm)
    if (is.na(cv) || cv <= 0) next

    ou_nm   <- fdm$ou[[nm]]
    has_wq  <- isTRUE(ou_nm$has_wq)
    c_norm  <- ou_nm$c_norm_obs_anch
    anch    <- if (has_wq) tm$models[[nm]]$d_anchors else tm$models[[nm]]$anchors

    if (is.null(anch) || is.null(c_norm) || length(c_norm) != nrow(anch)) next

    sigma_ln <- sqrt(log(1 + cv^2))
    mult     <- exp(stats::rnorm(nrow(anch), -sigma_ln^2 / 2, sigma_ln))
    anch_p   <- anch
    anch_p$S <- anch$S + c_norm * (mult - 1)

    if (has_wq) {
      tm$models[[nm]]$d_anchors <- anch_p
    } else {
      tm$models[[nm]]$anchors <- anch_p
    }
  }
  tm
}


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
    if (length(grab_cv_co) == 1L) return(as.numeric(grab_cv_co))
    if (a %in% names(grab_cv_co)) return(as.numeric(grab_cv_co[[a]]))
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
      a    <- co_day$analyte[i]
      cgm  <- fdm$co_grab_map[[a]]
      if (is.null(cgm)) next
      src_d <- cgm[d]
      if (is.na(src_d)) next
      key <- paste(a, src_d, sep = "::")
      m   <- mults[[key]]
      if (!is.null(m) && is.finite(m) && m > 0) co_day$value[i] <- co_day$value[i] * m
    }
    co_day
  })
  names(co_split_d) <- names(fdm$co_split)

  ## Rebuild wq_long (needed for WQ-layer PC-score computation) if present.
  wq_long_d <- if (!is.null(fdm$wq_long)) {
    rows <- lapply(names(co_split_d), function(d) {
      co_d <- co_split_d[[d]]
      if (nrow(co_d) == 0L) return(NULL)
      data.frame(sample_id = d, analyte = co_d$analyte, value = co_d$value,
                 stringsAsFactors = FALSE)
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows) > 0L) tibble::as_tibble(do.call(rbind, rows)) else NULL
  } else NULL

  list(co_split = co_split_d, wq_long = wq_long_d)
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
                                   tox_analytes) {
  fdm <- .fit_daily_target(
    site_rows        = site_rows,
    reference_model  = reference_model,
    imputation_model = imputation_model,
    conc_units       = conc_units,
    meta             = meta,
    tox_analytes     = tox_analytes,
    daily_long       = daily_long
  )
  if (is.null(fdm)) return(daily_long)

  model_rows <- .predict_daily_tox(fdm)
  if (is.null(model_rows)) return(daily_long)

  modelled <- fdm$modelled
  keep <- dplyr::filter(daily_long, !.data$analyte %in% .env$modelled)
  if (!"units.analyte" %in% names(keep)) keep$units.analyte <- NA_character_
  out <- dplyr::bind_rows(keep, model_rows)
  attr(out, "impact_tiers") <- attr(model_rows, "impact_tiers")
  out
}


#' Empty tibble matching the amspaf_daily() return schema
#' @keywords internal
.empty_daily_result <- function() {
  tibble::tibble(
    date                  = as.Date(character()),
    site_id               = character(),
    amspaf                = numeric(),
    n_analytes_used       = integer(),
    dominant_analyte      = character(),
    max_paf               = numeric(),
    n_measured_analytes   = integer(),
    days_since_last_sample = integer(),
    analyte_pafs          = list()
  )
}
