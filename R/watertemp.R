#' Estimate water temperature from air temperature measurements
#'
#' Fits a simple linear regression between observed air temperature and
#' observed water temperature, then predicts water temperature for a set of
#' target dates. The result can be added to a chemistry data frame as
#' `analyte = "temperature"` rows, enabling the NH₃-N pH/temperature
#' normalisation in [add_amspaf()].
#'
#' @section Air temperature variable:
#' This function requires **mean daily air temperature** (°C). Mean daily air
#' temperature is the standard predictor for daily mean water temperature in
#' the literature and has the best empirical performance for streams and ponds.
#' If your weather station or BOM data only provides daily maximum and minimum,
#' compute the mean as `(Tmax + Tmin) / 2`.
#'
#' @section Why water temperature matters:
#' The ANZG NH₃-N guideline values are derived at a specific reference
#' condition (pH 7.0, 20 °C total ammonia-N). The normalisation formula
#' adjusts measured NH₃-N to these conditions before SSD lookup — and since
#' the un-ionised fraction changes steeply with temperature and pH (roughly
#' doubling for each 5 °C increase near 20 °C), using the wrong temperature
#' has a large effect on the normalised concentration. A default value is not
#' provided: you must supply site-appropriate temperature estimates.
#'
#' @param air_temp_df Data frame of air temperature observations. Required
#'   columns:
#'   - `datetime` (Date or POSIXct) — measurement date
#'   - `air_temp_mean_C` (numeric) — mean daily air temperature in Celsius
#' @param water_temp_obs Data frame of observed water temperatures, used to
#'   calibrate the air-water regression. Required columns:
#'   - `datetime` (Date or POSIXct)
#'   - `water_temp_C` (numeric) — observed water temperature in Celsius
#'   At least 10 paired observations are recommended. Cannot be `NULL` —
#'   a regression requires training data. If you have no water temperature
#'   observations at all, collect a season of spot measurements alongside
#'   air temperatures before proceeding.
#' @param target_dates Date vector of dates for which to predict water
#'   temperature. Defaults to all dates present in `air_temp_df`. Dates
#'   with no air temperature observation are excluded with a warning.
#' @param lag_days Integer. Apply a lag to air temperature before regressing:
#'   the water temperature on day `t` is predicted from the air temperature on
#'   day `t - lag_days`. A lag of 1–3 days is appropriate for streams; 0
#'   (default) is appropriate for ponds or when using aggregated data.
#' @param site_id Character. Value to use for the `site_id` column in the
#'   returned chemistry rows. Default `NA_character_`.
#'
#' @return A tibble suitable for binding onto a chemistry data frame, with
#'   columns:
#'   - `datetime` (Date)
#'   - `analyte` — `"temperature"`
#'   - `value` — predicted water temperature (°C)
#'   - `detected` — `TRUE`
#'   - `site_id` — from `site_id` argument
#'   - `sample_id` — `NA_character_` (set by caller if needed)
#'   Attach `attr(result, "lm_fit")` — the fitted `lm` object for inspection.
#'
#' @examples
#' \dontrun{
#' air  <- tibble::tibble(
#'   datetime       = seq(as.Date("2020-01-01"), as.Date("2023-12-31"), by = "day"),
#'   air_temp_mean_C = rnorm(length(datetime), mean = 15, sd = 8)
#' )
#' wt_obs <- tibble::tibble(
#'   datetime      = sample(air$datetime, 80),
#'   water_temp_C  = air$air_temp_mean_C[match(datetime, air$datetime)] * 0.85 + 2 + rnorm(80, 0, 1)
#' )
#' wt <- estimate_water_temp(air, wt_obs)
#' # Add to chemistry df:
#' # chemistry <- dplyr::bind_rows(chemistry, dplyr::mutate(wt, sample_id = ...))
#' }
#'
#' @export
estimate_water_temp <- function(
    air_temp_df,
    water_temp_obs,
    target_dates  = NULL,
    lag_days      = 0L,
    site_id       = NA_character_
) {
  checkmate::assert_data_frame(air_temp_df)
  checkmate::assert_names(names(air_temp_df),
    must.include = c("datetime", "air_temp_mean_C"))
  checkmate::assert_data_frame(water_temp_obs)
  checkmate::assert_names(names(water_temp_obs),
    must.include = c("datetime", "water_temp_C"))
  checkmate::assert_int(lag_days, lower = 0L)
  checkmate::assert_character(site_id, len = 1L)

  # Normalise dates
  air_df <- dplyr::mutate(air_temp_df,
    .date = as.Date(.data$datetime)
  ) |>
    dplyr::select(".date", "air_temp_mean_C") |>
    dplyr::distinct(.data$.date, .keep_all = TRUE) |>
    dplyr::arrange(.data$.date)

  wt_df <- dplyr::mutate(water_temp_obs,
    .date = as.Date(.data$datetime)
  ) |>
    dplyr::select(".date", "water_temp_C")

  n_wt <- nrow(wt_df)
  if (n_wt < 5L) {
    cli::cli_abort(c(
      "{.arg water_temp_obs} has only {n_wt} row{?s}.",
      "i" = "At least 5 paired observations are required to fit a regression. \\
             Recommend >= 10 covering different seasons."
    ))
  }
  if (n_wt < 10L) {
    cli::cli_warn(c(
      "!" = "{.arg water_temp_obs} has only {n_wt} observation{?s}.",
      "i" = "A regression calibrated on fewer than 10 paired observations \\
             may have poor predictive accuracy."
    ))
  }

  # Apply lag: shift air temp data by lag_days relative to water temp
  if (lag_days > 0L) {
    air_lagged <- dplyr::mutate(air_df, .date = .data$.date + lag_days)
  } else {
    air_lagged <- air_df
  }

  # Merge training data on matching dates
  train_df <- dplyr::inner_join(wt_df, air_lagged, by = ".date")

  if (nrow(train_df) < 5L) {
    cli::cli_abort(c(
      "After date matching{if (lag_days > 0L) paste0(' (lag = ', lag_days, ' d)') else ''}, \\
       only {nrow(train_df)} paired observations are available.",
      "i" = "Ensure {.arg air_temp_df} and {.arg water_temp_obs} cover overlapping \\
             date ranges."
    ))
  }

  # Fit linear regression: water_temp_C ~ air_temp_mean_C
  fit <- stats::lm(water_temp_C ~ air_temp_mean_C, data = train_df)

  r2   <- summary(fit)$r.squared
  rmse <- sqrt(mean(stats::residuals(fit)^2))
  cli::cli_inform(c(
    "i" = "Air-water temperature regression: R² = {round(r2, 3)}, \\
           RMSE = {round(rmse, 2)} °C \\
           (n = {nrow(train_df)} paired observations{if (lag_days > 0L) paste0(', lag = ', lag_days, ' d') else ''})."
  ))

  if (r2 < 0.70) {
    cli::cli_warn(c(
      "!" = "Air-water temperature regression R² = {round(r2, 3)} (below 0.70).",
      "i" = "Consider collecting more paired observations or checking for \\
             unusual thermal conditions at your site."
    ))
  }

  # Build prediction data frame
  if (is.null(target_dates)) {
    pred_dates <- air_df$.date
  } else {
    pred_dates <- as.Date(target_dates)
  }

  # Apply lag to prediction dates
  pred_air_dates <- if (lag_days > 0L) {
    pred_dates - lag_days
  } else {
    pred_dates
  }

  # Locate air temp for each prediction date
  pred_df <- tibble::tibble(
    datetime         = pred_dates,
    .air_date        = pred_air_dates
  ) |>
    dplyr::left_join(
      dplyr::rename(air_df, .air_date = ".date"),
      by = ".air_date"
    )

  n_missing_air <- sum(is.na(pred_df$air_temp_mean_C))
  if (n_missing_air > 0L) {
    cli::cli_warn(c(
      "!" = "{n_missing_air} target date{?s} have no matching air temperature \\
             observation and will be dropped from the output."
    ))
    pred_df <- dplyr::filter(pred_df, !is.na(.data$air_temp_mean_C))
  }

  # Predict
  predicted_wt <- stats::predict(fit, newdata = pred_df)

  result <- tibble::tibble(
    datetime  = pred_df$datetime,
    analyte   = "temperature",
    value     = as.numeric(predicted_wt),
    detected  = TRUE,
    site_id   = site_id,
    sample_id = NA_character_
  )

  attr(result, "lm_fit") <- fit
  result
}
