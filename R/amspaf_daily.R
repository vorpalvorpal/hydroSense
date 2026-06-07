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
    if (interpolation == "model") {
      daily_long <- .daily_tox_from_model(
        daily_long       = daily_long,
        site_rows        = site_rows,
        reference_model  = reference_model,
        imputation_model = imputation_model,
        conc_units       = conc_units,
        meta             = meta,
        tox_analytes     = tox_analytes
      )
    }

    ## Step 3: Compute per-day diagnostics from SSD-eligible rows.
    diag <- .compute_daily_diag(daily_long, tox_analytes, site)

    ## Step 4: Build synthetic long-format samples (no focal_date!).
    synth <- .build_synthetic_samples(daily_long, site)

    list(synth = synth, diag = diag)
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


#' Model-based daily toxicant interpolation (season-blind impact model)
#'
#' Fits a [fit_target_model()] on the site's grab chemistry and the supplied
#' `reference_model`, predicts each toxicant's normalised concentration
#' `C_norm = ref_norm + impact` on every daily date, and reconstructs the raw
#' µg/L concentration `C_raw = C_norm / factor` (the normalisation is
#' multiplicative, so the factor is `normalise(1, co-analytes_of_the_day)`).
#' Modelled-toxicant rows in `daily_long` are replaced with these estimates
#' (in µg/L, `units.analyte = "ug/L"`); co-analytes and non-modelled toxicants
#' are left untouched. On any failure the input is returned unchanged.
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
  if (is.null(tm) || length(tm$models) == 0L) return(daily_long)

  modelled <- names(tm$models)
  qdates   <- sort(unique(daily_long$.date))
  pred <- .resolve_target_impact(tm, tibble::tibble(date = qdates), modelled)
  if (nrow(pred) == 0L) return(daily_long)

  ## Per-date co-analyte lookup (named numeric vector keyed by analyte).
  co <- daily_long[!daily_long$analyte %in% tox_analytes &
                     (is.na(daily_long$detected) | daily_long$detected), , drop = FALSE]
  co_split <- split(
    data.frame(analyte = co$analyte, value = co$value, stringsAsFactors = FALSE),
    as.character(co$.date)
  )
  co_vec_for <- function(d) {
    cd <- co_split[[as.character(d)]]
    if (is.null(cd)) return(numeric(0))
    stats::setNames(cd$value, cd$analyte)
  }

  ## Normalisation factor = normalise(1, co-analytes) per analyte (multiplicative).
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

  ## Measured (grab) dates per modelled analyte, for the .measured flag.
  sr_mod <- site_rows[site_rows$analyte %in% modelled &
                        (is.na(site_rows$detected) | site_rows$detected), ,
                      drop = FALSE]
  measured_key <- paste(sr_mod$analyte, as.Date(sr_mod$datetime))

  ## Reconstruct raw µg/L per (date, analyte).
  pred$C_raw <- vapply(seq_len(nrow(pred)), function(i) {
    a  <- pred$analyte[i]; d <- pred$date[i]
    parsed <- fac_lookup[[a]]$parsed
    if (is.null(parsed)) return(pred$C_norm[i])
    factor <- .apply_normalisation(parsed, 1, co_vec_for(d))
    if (is.na(factor) || factor <= 0) return(NA_real_)
    pred$C_norm[i] / factor
  }, numeric(1L))

  pred_ok <- pred[is.finite(pred$C_raw), , drop = FALSE]
  model_rows <- tibble::tibble(
    .date         = pred_ok$date,
    value         = pred_ok$C_raw,
    detected      = TRUE,
    .measured     = paste(pred_ok$analyte, pred_ok$date) %in% measured_key,
    analyte       = pred_ok$analyte,
    units.analyte = "ug/L"
  )

  ## Replace modelled-toxicant rows; keep co-analytes + non-modelled toxicants.
  keep <- dplyr::filter(daily_long, !.data$analyte %in% .env$modelled)
  if (!"units.analyte" %in% names(keep)) keep$units.analyte <- NA_character_
  dplyr::bind_rows(keep, model_rows)
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
