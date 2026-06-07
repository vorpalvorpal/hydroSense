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
#'   At least 10 paired observations are recommended; fewer than 5 triggers an
#'   error unless `fallback` is also supplied. `NULL` is accepted when `fallback`
#'   is provided.
#' @param target_dates Date vector of dates for which to predict water
#'   temperature. Defaults to all dates present in `air_temp_df`. Dates
#'   with no air temperature observation are excluded with a warning.
#' @param lag_days Integer. Apply a lag to air temperature before regressing:
#'   the water temperature on day `t` is predicted from the air temperature on
#'   day `t - lag_days`. A lag of 1–3 days is appropriate for streams; 0
#'   (default) is appropriate for ponds or when using aggregated data.
#' @param site_id Character. Value to use for the `site_id` column in the
#'   returned chemistry rows. Default `NA_character_`.
#' @param seasonal One of `"auto"` (default), `"off"`, or `"on"`. Controls
#'   whether a day-of-year seasonal term is added to the air-temperature
#'   regression (see *Model selection*). `"auto"` fits both the air-only and
#'   air + season models and keeps whichever has the lower AICc; `"off"`
#'   forces the air-only model (the legacy behaviour); `"on"` forces the
#'   seasonal model whenever it is eligible.
#' @param seasonal_min_n Integer. Minimum number of paired observations before
#'   the seasonal model is even considered (default 8). The seasonal model adds
#'   two parameters, so a small buffer above that is needed for AICc to behave.
#' @param seasonal_min_quarters Integer 1-4. Minimum number of distinct
#'   calendar quarters the training observations must span before the seasonal
#'   model is considered (default 3). Calibration data confined to one or two
#'   seasons cannot anchor an annual cycle and would extrapolate badly.
#' @param fallback Air-to-water relationship used when there are too few paired
#'   observations to fit a regression (fewer than 5 rows in `water_temp_obs`,
#'   or `water_temp_obs = NULL`). Accepted forms:
#'   - `"identity"` — `water = air` (simplest default for small streams).
#'   - `c(intercept, slope)` — `water = intercept + slope * air`.
#'   - A function `function(air) ...` — arbitrary transformation.
#'   `NULL` (default) preserves the existing behaviour: fewer than 5 pairs is
#'   a hard error.
#' @param select Antecedence window for the air-temperature regressor, in days.
#'   The regressor is the trailing `select`-day mean of air temperature.
#'   - `"auto"` (default) — select the window by leave-one-out cross-validation
#'     (PRESS statistic, exact for OLS) with the **1-SE rule**: the shortest
#'     window whose CV RMSE is within one SE of the best is chosen, to avoid
#'     over-fitting to noise. Requires at least `auto_min_n` pairs; falls back
#'     to 1 (or to the `water_body_type` default) when pairs are insufficient.
#'   - A positive integer — use that fixed window (1 = daily air temp,
#'     equivalent to the original behaviour).
#' @param auto_min_n Integer. Minimum number of paired observations required
#'   before the `select = "auto"` window search is attempted (default 20).
#'   Below this threshold the window falls back to 1 (or `water_body_type`).
#' @param water_body_type Character. Physical type of the water body, used as
#'   the fallback window when `select = "auto"` but `n_pairs < auto_min_n`.
#'   One of `"stream_small"` (1 d), `"stream_large"` (5 d), `"estuary"` (7 d),
#'   `"pond"` (10 d), `"lake"` (30 d), `"reservoir"` (60 d), `"groundwater"`
#'   (300 d). `NULL` (default) falls back to window = 1. Ignored when the CV
#'   search runs.
#'
#' @section Model selection:
#' **Window selection** (antecedence timescale): water bodies with more thermal
#' mass integrate air temperature over longer periods. The function searches
#' log-spaced trailing windows from 1 to 300 days, fitting a simple OLS for
#' each, and selects by leave-one-out CV RMSE (PRESS statistic — closed form,
#' no refitting loop). The **1-SE rule** then chooses the shortest window within
#' one standard error of the best, avoiding spurious preference for a long window
#' that barely beats shorter alternatives. All window candidates have the same
#' number of parameters, so AICc cannot distinguish them; CV is required.
#'
#' **Seasonal term** (thermal hysteresis): water temperature lags air
#' temperature seasonally — at the same air temperature, water tends to be
#' cooler in spring and warmer in autumn. Adding a **first-harmonic day-of-year
#' term** — `sin(2*pi*doy/365.25)` and `cos(2*pi*doy/365.25)` — lets the
#' regression represent it. Under `seasonal = "auto"` both the window-regressor-
#' only and window-regressor + harmonics models are compared by **AICc**; the
#' lower-AICc model is used. AICc is appropriate here because the two candidates
#' differ in complexity (2 vs 4 parameters). The seasonal model is only eligible
#' when there are at least `seasonal_min_n` observations spanning at least
#' `seasonal_min_quarters` quarters.
#'
#' @return A tibble suitable for binding onto a chemistry data frame, with
#'   columns:
#'   - `datetime` (Date)
#'   - `analyte` — `"temperature"`
#'   - `value` — predicted water temperature (°C)
#'   - `detected` — `TRUE`
#'   - `site_id` — from `site_id` argument
#'   - `sample_id` — `NA_character_` (set by caller if needed)
#'   Attributes attached for inspection:
#'   - `attr(result, "lm_fit")` — the selected fitted `lm` object (`NULL` for fallback).
#'   - `attr(result, "model")` — label of the selected model (`"fallback"` when
#'     no pairs were available).
#'   - `attr(result, "seasonal_used")` — `TRUE` if the seasonal model won.
#'   - `attr(result, "model_comparison")` — data frame of AICc / R² per
#'     candidate model and which was selected.
#'   - `attr(result, "window_selected")` — chosen antecedence window in days
#'     (`NA` for fallback path).
#'   - `attr(result, "window_comparison")` — data frame of LOO-CV RMSE per
#'     candidate window (`NULL` when window was fixed or search was skipped).
#'
#' @examples
#' set.seed(1)
#' air <- tibble::tibble(
#'   datetime        = seq(as.Date("2020-01-01"), as.Date("2022-12-31"), by = "day")
#' )
#' air$air_temp_mean_C <- rnorm(nrow(air), mean = 15, sd = 8)
#' wt_obs <- tibble::tibble(
#'   datetime     = sample(air$datetime, 80),
#'   water_temp_C = NA_real_
#' )
#' wt_obs$water_temp_C <-
#'   air$air_temp_mean_C[match(wt_obs$datetime, air$datetime)] * 0.85 + 2 +
#'   rnorm(80, 0, 1)
#' wt <- estimate_water_temp(air, wt_obs)
#' attr(wt, "model")  # which model was selected (air-only vs air + season)
#'
#' @seealso [get_silo_air_temp()] to source `air_temp_df` from SILO for an
#'   Australian location; [add_amspaf()], which requires the resulting water
#'   `temperature` rows for ammonia.
#' @export
estimate_water_temp <- function(
    air_temp_df,
    water_temp_obs        = NULL,
    target_dates          = NULL,
    lag_days              = 0L,
    site_id               = NA_character_,
    seasonal              = c("auto", "off", "on"),
    seasonal_min_n        = 8L,
    seasonal_min_quarters = 3L,
    fallback              = NULL,
    select                = "auto",
    auto_min_n            = 20L,
    water_body_type       = NULL
) {
  seasonal <- match.arg(seasonal)
  checkmate::assert_data_frame(air_temp_df)
  checkmate::assert_names(names(air_temp_df),
    must.include = c("datetime", "air_temp_mean_C"))
  if (!is.null(water_temp_obs)) {
    checkmate::assert_data_frame(water_temp_obs)
    checkmate::assert_names(names(water_temp_obs),
      must.include = c("datetime", "water_temp_C"))
  }
  checkmate::assert_int(lag_days, lower = 0L)
  checkmate::assert_character(site_id, len = 1L)
  checkmate::assert_int(seasonal_min_n, lower = 7L)
  checkmate::assert_int(seasonal_min_quarters, lower = 1L, upper = 4L)
  checkmate::assert_int(auto_min_n, lower = 1L)
  if (!is.null(fallback)) {
    checkmate::assert(
      checkmate::check_string(fallback, pattern = "^identity$"),
      checkmate::check_numeric(fallback, len = 2, any.missing = FALSE),
      checkmate::check_function(fallback),
      .var.name = "fallback"
    )
  }
  if (is.character(select)) {
    select <- match.arg(select, "auto")
  } else {
    checkmate::assert_int(as.integer(select), lower = 1L)
    select <- as.integer(select)
  }
  if (!is.null(water_body_type)) {
    checkmate::assert_choice(water_body_type, names(.WATER_BODY_WINDOWS))
  }

  # Normalise dates
  air_df <- dplyr::mutate(air_temp_df,
    .date = as.Date(.data$datetime)
  ) |>
    dplyr::select(".date", "air_temp_mean_C") |>
    dplyr::distinct(.data$.date, .keep_all = TRUE) |>
    dplyr::arrange(.data$.date)

  if (!is.null(water_temp_obs)) {
    wt_df <- dplyr::mutate(water_temp_obs,
      .date = as.Date(.data$datetime)
    ) |>
      dplyr::select(".date", "water_temp_C")
    n_wt <- nrow(wt_df)
  } else {
    wt_df <- NULL
    n_wt  <- 0L
  }

  use_fallback <- n_wt < 5L
  if (use_fallback && is.null(fallback)) {
    cli::cli_abort(c(
      "{.arg water_temp_obs} has only {n_wt} row{?s}.",
      "i" = "At least 5 paired observations are required to fit a regression. \\
             Recommend >= 10 covering different seasons.",
      "i" = "Alternatively, supply a {.arg fallback} relationship \\
             ({.val \"identity\"}, {.code c(intercept, slope)}, or a function)."
    ))
  }
  if (!use_fallback && n_wt < 10L) {
    cli::cli_warn(c(
      "!" = "{.arg water_temp_obs} has only {n_wt} observation{?s}.",
      "i" = "A regression calibrated on fewer than 10 paired observations \\
             may have poor predictive accuracy."
    ))
  }

  # \u2500\u2500 FALLBACK PATH \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
  if (use_fallback) {
    pred_fn <- .build_fallback_fn(fallback)

    if (is.null(target_dates)) {
      pred_dates <- air_df$.date
    } else {
      pred_dates <- as.Date(target_dates)
    }
    pred_air_dates <- if (lag_days > 0L) pred_dates - lag_days else pred_dates

    pred_df <- tibble::tibble(
      datetime  = pred_dates,
      .air_date = pred_air_dates
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

    cli::cli_inform(c("i" = "Using {.arg fallback} air-water relationship (no paired observations)."))

    result <- tibble::tibble(
      datetime  = pred_df$datetime,
      analyte   = "temperature",
      value     = as.numeric(pred_fn(pred_df$air_temp_mean_C)),
      detected  = TRUE,
      site_id   = site_id,
      sample_id = NA_character_
    )
    attr(result, "lm_fit")        <- NULL
    attr(result, "model")         <- "fallback"
    attr(result, "seasonal_used") <- FALSE
    attr(result, "model_comparison") <- data.frame(
      model    = c("air_only", "air_plus_season"),
      aicc     = c(NA_real_, NA_real_),
      r2       = c(NA_real_, NA_real_),
      selected = c(FALSE, FALSE)
    )
    attr(result, "window_selected")  <- NA_integer_
    attr(result, "window_comparison") <- NULL
    return(result)
  }

  # \u2500\u2500 REGRESSION PATH \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500

  # Apply lag: shift air temp data by lag_days relative to water temp. Rolling
  # mean is computed on the original air series before the lag shift, so the
  # w-day mean ending on the air date is then aligned to the target water date.
  air_for_roll <- air_df   # rolling mean always computed on actual air dates
  air_lagged <- if (lag_days > 0L) {
    dplyr::mutate(air_df, .date = .data$.date + lag_days)
  } else {
    air_df
  }

  # \u2500\u2500 Window selection \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
  window_comp_df <- NULL
  if (is.character(select) && select == "auto") {
    if (n_wt >= auto_min_n) {
      cv_rows <- lapply(.WINDOW_CANDIDATES, function(w) {
        roll_vals <- .trailing_rollmean(air_for_roll$air_temp_mean_C, w)
        air_roll  <- dplyr::mutate(air_lagged, .roll = roll_vals)
        tr <- dplyr::inner_join(wt_df, air_roll, by = ".date") |>
              dplyr::filter(!is.na(.data$.roll))
        if (nrow(tr) < 5L) return(NULL)
        fit_w <- stats::lm(water_temp_C ~ .roll, data = tr)
        cv    <- .loo_rmse_and_se(fit_w)
        list(window = w, loo_rmse = cv$rmse, loo_rmse_se = cv$se, n_pairs = nrow(tr))
      })
      cv_rows <- Filter(Negate(is.null), cv_rows)

      if (length(cv_rows) == 0L) {
        selected_window <- 1L
      } else {
        cv_df     <- dplyr::bind_rows(cv_rows)
        best_idx  <- which.min(cv_df$loo_rmse)
        # 1-SE rule: prefer the shortest window whose CV RMSE is within one SE
        # of the best (SE in RMSE units via delta method). This avoids chasing
        # marginal gains into physically implausible long windows.
        threshold  <- cv_df$loo_rmse[best_idx] + cv_df$loo_rmse_se[best_idx]
        eligible   <- cv_df$window[cv_df$loo_rmse <= threshold]
        selected_window <- min(eligible)
        moved_off_min   <- selected_window != cv_df$window[best_idx]
        n_sel <- cv_df$n_pairs[cv_df$window == selected_window]
        r_sel <- round(cv_df$loo_rmse[cv_df$window == selected_window], 2)
        cli::cli_inform(c(
          "i" = "Selected antecedence window: {selected_window} d \\
                 (LOO-CV RMSE = {r_sel} \u00b0C, n = {n_sel}\\
                 {if (moved_off_min) '; 1-SE rule moved choice off minimum' else ''})."
        ))
        window_comp_df <- data.frame(
          window      = cv_df$window,
          loo_rmse    = cv_df$loo_rmse,
          loo_rmse_se = cv_df$loo_rmse_se,
          n_pairs     = cv_df$n_pairs,
          selected    = cv_df$window == selected_window
        )
      }
    } else {
      selected_window <- if (!is.null(water_body_type)) {
        .WATER_BODY_WINDOWS[[water_body_type]]
      } else {
        1L
      }
      type_note <- if (!is.null(water_body_type)) {
        paste0(" (", water_body_type, " default)")
      } else {
        ""
      }
      cli::cli_inform(c(
        "i" = "Window search skipped: {n_wt} pair{?s} < auto_min_n = {auto_min_n}; \\
               using window = {selected_window} d{type_note}."
      ))
    }
  } else {
    selected_window <- select   # fixed integer window
  }

  # \u2500\u2500 Build training data with the selected rolling-mean regressor \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
  roll_vals_train <- .trailing_rollmean(air_for_roll$air_temp_mean_C, selected_window)
  air_lagged_roll <- dplyr::mutate(air_lagged, .roll = roll_vals_train)

  train_df <- dplyr::inner_join(wt_df, air_lagged_roll, by = ".date") |>
              dplyr::filter(!is.na(.data$.roll))

  if (nrow(train_df) < 5L) {
    extra <- if (selected_window > 1L) {
      paste0(" and ", selected_window, "-day rolling mean")
    } else {
      ""
    }
    cli::cli_abort(c(
      "After date matching{if (lag_days > 0L) paste0(' (lag = ', lag_days, ' d)') else ''} \\
       {extra}, only {nrow(train_df)} paired observations are available.",
      "i" = "Ensure {.arg air_temp_df} and {.arg water_temp_obs} cover overlapping \\
             date ranges."
    ))
  }

  # Day-of-year first-harmonic terms (cyclic seasonal signal).
  train_df <- .add_doy_harmonics(train_df, ".date")

  # Candidate 1: rolling-mean air temp only (the baseline; always fitted).
  fit_air <- stats::lm(water_temp_C ~ .roll, data = train_df)

  # Candidate 2: rolling-mean air + first-harmonic day-of-year.
  n_train  <- nrow(train_df)
  quarters <- length(unique((as.integer(format(train_df$.date, "%m")) - 1L) %/% 3L))
  seasonal_eligible <-
    n_train >= seasonal_min_n && quarters >= seasonal_min_quarters

  fit_seas <- NULL
  if (seasonal != "off" && seasonal_eligible) {
    fit_seas <- stats::lm(
      water_temp_C ~ .roll + sin_doy + cos_doy, data = train_df)
  }

  # Select by AICc (small-sample-corrected). Prefer the seasonal model only if
  # it genuinely improves AICc (or seasonal = "on" forces it when eligible).
  aicc_air  <- .aicc(fit_air)
  aicc_seas <- if (!is.null(fit_seas)) .aicc(fit_seas) else NA_real_
  use_seasonal <- !is.null(fit_seas) && is.finite(aicc_seas) &&
    (seasonal == "on" || aicc_seas < aicc_air)

  fit       <- if (use_seasonal) fit_seas else fit_air
  model_lbl <- if (use_seasonal) "air + season (day-of-year harmonic)" else "air only"

  r2   <- summary(fit)$r.squared
  rmse <- sqrt(mean(stats::residuals(fit)^2))
  cli::cli_inform(c(
    "i" = "Selected air-water model: {model_lbl}. R\u00b2 = {round(r2, 3)}, \\
           RMSE = {round(rmse, 2)} \u00b0C, AICc = {round(if (use_seasonal) aicc_seas else aicc_air, 1)} \\
           (n = {n_train}{if (lag_days > 0L) paste0(', lag = ', lag_days, ' d') else ''})."
  ))
  if (!is.null(fit_seas)) {
    cli::cli_inform(c(
      " " = "AICc comparison (lower preferred): air only = {round(aicc_air, 1)}, \\
             air + season = {round(aicc_seas, 1)}."
    ))
  } else if (seasonal != "off" && !seasonal_eligible) {
    cli::cli_inform(c(
      "i" = "Seasonal model not considered: needs >= {seasonal_min_n} observations \\
             across >= {seasonal_min_quarters} quarters (have {n_train} across {quarters})."
    ))
  }

  if (r2 < 0.70) {
    cli::cli_warn(c(
      "!" = "Air-water temperature model R\u00b2 = {round(r2, 3)} (below 0.70).",
      "i" = "Consider collecting more paired observations or checking for \\
             unusual thermal conditions at your site."
    ))
  }

  # \u2500\u2500 Prediction \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
  if (is.null(target_dates)) {
    pred_dates <- air_df$.date
  } else {
    pred_dates <- as.Date(target_dates)
  }

  # Rolling mean for predictions is computed on the FULL air series so target
  # dates beyond the training range have full antecedent history.
  pred_air_dates <- if (lag_days > 0L) pred_dates - lag_days else pred_dates
  roll_vals_pred  <- .trailing_rollmean(air_for_roll$air_temp_mean_C, selected_window)
  air_for_roll_df <- dplyr::mutate(air_for_roll, .roll = roll_vals_pred)

  pred_df <- tibble::tibble(
    datetime  = pred_dates,
    .air_date = pred_air_dates
  ) |>
    dplyr::left_join(
      dplyr::rename(air_for_roll_df, .air_date = ".date"),
      by = ".air_date"
    )

  # A target date can be missing either because air data is absent or because
  # the rolling window cannot yet be computed (first w-1 days of air series).
  n_missing_air <- sum(is.na(pred_df$.roll))
  if (n_missing_air > 0L) {
    cli::cli_warn(c(
      "!" = "{n_missing_air} target date{?s} have no matching air temperature \\
             observation and will be dropped from the output."
    ))
    pred_df <- dplyr::filter(pred_df, !is.na(.data$.roll))
  }

  # Predict (seasonal terms indexed to the target/water date)
  pred_df      <- .add_doy_harmonics(pred_df, "datetime")
  predicted_wt <- stats::predict(fit, newdata = pred_df)

  result <- tibble::tibble(
    datetime  = pred_df$datetime,
    analyte   = "temperature",
    value     = as.numeric(predicted_wt),
    detected  = TRUE,
    site_id   = site_id,
    sample_id = NA_character_
  )

  attr(result, "lm_fit")        <- fit            # chosen model (back-compat)
  attr(result, "model")         <- model_lbl
  attr(result, "seasonal_used") <- use_seasonal
  attr(result, "model_comparison") <- data.frame(
    model    = c("air_only", "air_plus_season"),
    aicc     = c(aicc_air, aicc_seas),
    r2       = c(summary(fit_air)$r.squared,
                 if (!is.null(fit_seas)) summary(fit_seas)$r.squared else NA_real_),
    selected = c(!use_seasonal, use_seasonal)
  )
  attr(result, "window_selected")  <- selected_window
  attr(result, "window_comparison") <- window_comp_df
  result
}

# Append first-harmonic day-of-year terms for a Date column. Day-of-year is
# cyclic (Dec 31 ~ Jan 1), so it enters a linear model as sin/cos of the annual
# angle rather than as a raw linear term.
.add_doy_harmonics <- function(df, date_col) {
  doy <- as.integer(format(as.Date(df[[date_col]]), "%j"))
  ang <- 2 * pi * doy / 365.25
  df$sin_doy <- sin(ang)
  df$cos_doy <- cos(ang)
  df
}

# Small-sample-corrected Akaike Information Criterion. k counts estimated
# parameters including the residual variance. Returns Inf when n is too small
# for the correction to be defined (n <= k + 1), so an over-parameterised model
# is never selected.
.aicc <- function(fit) {
  n <- length(stats::residuals(fit))
  k <- length(stats::coef(fit)) + 1L
  if (n - k - 1L <= 0L) return(Inf)
  stats::AIC(fit) + (2 * k * (k + 1)) / (n - k - 1L)
}

# Build a scalar prediction function from the three fallback forms.
.build_fallback_fn <- function(fallback) {
  if (is.function(fallback)) return(fallback)
  if (is.numeric(fallback)) {
    a <- fallback[1]; b <- fallback[2]
    return(function(air) a + b * air)
  }
  function(air) air   # "identity"
}

# Trailing w-day rolling mean using base stats::filter (no extra dependency).
# Returns NA for the first w-1 positions where history is insufficient.
.trailing_rollmean <- function(x, w) {
  if (w == 1L) return(x)
  as.numeric(stats::filter(x, rep(1 / w, w), sides = 1))
}

# LOO-CV RMSE and its SE via the PRESS statistic (closed form for OLS — no
# refitting loop). e_i = r_i / (1 - h_i) is the exact leave-one-out residual.
# SE of RMSE is approximated via the delta method: SE(sqrt(MSE)) ≈ SE(MSE)/(2*RMSE).
.loo_rmse_and_se <- function(fit) {
  r  <- stats::residuals(fit)
  h  <- stats::hatvalues(fit)
  e2 <- (r / (1 - h))^2
  mse    <- mean(e2)
  se_mse <- sqrt(stats::var(e2) / length(e2))
  rmse   <- sqrt(mse)
  list(rmse = rmse, se = se_mse / (2 * rmse))
}

# Log-spaced candidate antecedence windows (days) for multi-timescale selection.
.WINDOW_CANDIDATES <- c(1, 2, 3, 4, 5, 7, 10, 15, 20, 25, 30, 40, 50, 60,
                         80, 90, 100, 120, 140, 160, 180, 220, 260, 300)

# Per-water-body-type default antecedence windows, used as fallback when
# select = "auto" but n_pairs < auto_min_n. Reflects known thermal-mass scales.
.WATER_BODY_WINDOWS <- c(
  stream_small =   1L,   # 1st/2nd-order stream; tracks daily air
  stream_large =   5L,   # larger river; modest thermal buffering
  estuary      =   7L,   # tidally mixed; roughly weekly scale
  pond         =  10L,   # shallow still water
  lake         =  30L,   # deep lake; month-scale thermal lag
  reservoir    =  60L,   # large dam / storage; 1–2 month memory
  groundwater  = 300L    # effectively annual-mean air temperature
)

# ── SILO air-temperature lookup ───────────────────────────────────────────────

#' Fetch daily mean air temperature from SILO for an Australian location
#'
#' Retrieves daily air temperature from the SILO Data Drill (a ~5 km gridded,
#' spatially interpolated climate surface covering Australia, 1889–present) for
#' a single latitude/longitude, and returns **daily mean air temperature**
#' (`(Tmax + Tmin) / 2`, °C) in exactly the shape [estimate_water_temp()]
#' expects as its `air_temp_df`. The typical workflow is:
#'
#' ```
#' air   <- get_silo_air_temp(lat, lon, start, end)   # SILO mean air temp
#' wt    <- estimate_water_temp(air, water_temp_obs)  # calibrate water = f(air)
#' chem  <- dplyr::bind_rows(chem, dplyr::mutate(wt, sample_id = ...))
#' ```
#'
#' This wraps [weatherOz::get_data_drill()]; the **weatherOz** package must be
#' installed (it is listed under `Suggests`). Results are cached on disk so
#' repeat calls for the same grid cell and date range do not re-hit the API.
#'
#' @section Attribution:
#' SILO data are © State of Queensland (Department of Environment, Science and
#' Innovation) and released under CC-BY 4.0. Cite SILO when you publish results
#' derived from this function.
#'
#' @section API key:
#' SILO requires an API key, which is simply your email address. By default it
#' is auto-detected via [weatherOz::get_key()] (from `.Renviron`/`.Rprofile`);
#' see the weatherOz documentation for one-time setup, or pass `api_key`
#' directly.
#'
#' @param latitude,longitude Numeric decimal-degree coordinates of the point of
#'   interest. Must fall within the SILO grid (approximately latitude -44 to
#'   -10, longitude 112 to 154). Snapped to the 0.05° grid by SILO.
#' @param start_date,end_date Start and end of the (inclusive) date range, as
#'   `Date` objects or `"YYYY-MM-DD"` strings.
#' @param api_key Character SILO API key (your email address). Default `NULL`
#'   defers to [weatherOz::get_key()].
#' @param cache Logical; cache the result on disk under
#'   `tools::R_user_dir("leachatetools", "cache")/silo`. Default `TRUE`.
#' @param refresh Logical; if `TRUE`, ignore and overwrite any cached result.
#'   Default `FALSE`.
#'
#' @return A tibble with one row per day:
#'   - `datetime` (`Date`)
#'   - `air_temp_mean_C` (numeric, °C)
#'   Ready to pass as `air_temp_df` to [estimate_water_temp()].
#'
#' @seealso [estimate_water_temp()]
#'
#' @examples
#' \dontrun{
#' air <- get_silo_air_temp(
#'   latitude = -33.87, longitude = 151.21,
#'   start_date = "2020-01-01", end_date = "2023-12-31"
#' )
#' }
#' @export
get_silo_air_temp <- function(
    latitude,
    longitude,
    start_date,
    end_date,
    api_key = NULL,
    cache   = TRUE,
    refresh = FALSE
) {
  if (!requireNamespace("weatherOz", quietly = TRUE)) {
    cli::cli_abort(c(
      "The {.pkg weatherOz} package is required to fetch SILO climate data.",
      "i" = "Install it with {.run install.packages(\"weatherOz\")}."
    ))
  }
  checkmate::assert_number(latitude,  lower = -44, upper = -10)
  checkmate::assert_number(longitude, lower = 112, upper = 154)
  checkmate::assert_flag(cache)
  checkmate::assert_flag(refresh)
  start_date <- as.Date(start_date)
  end_date   <- as.Date(end_date)
  if (is.na(start_date) || is.na(end_date)) {
    cli::cli_abort("{.arg start_date} and {.arg end_date} must be valid dates.")
  }
  if (end_date < start_date) {
    cli::cli_abort("{.arg end_date} ({end_date}) is before {.arg start_date} ({start_date}).")
  }

  # ── Disk cache keyed by grid cell (0.05°) + date range ──────────────────────
  cache_path <- NULL
  if (cache) {
    cache_dir <- file.path(tools::R_user_dir("leachatetools", "cache"), "silo")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    key <- sprintf(
      "silo_%+07.2f_%+07.2f_%s_%s.qs",
      latitude, longitude, format(start_date, "%Y%m%d"), format(end_date, "%Y%m%d")
    )
    cache_path <- file.path(cache_dir, key)
    if (!refresh && file.exists(cache_path)) {
      return(qs2::qs_read(cache_path))
    }
  }

  if (is.null(api_key)) api_key <- weatherOz::get_key(service = "SILO")

  dt <- weatherOz::get_data_drill(
    longitude  = longitude,
    latitude   = latitude,
    start_date = start_date,
    end_date   = end_date,
    values     = c("max_temp", "min_temp"),
    api_key    = api_key
  )

  # weatherOz returns the daily temperature columns as `air_tmax`/`air_tmin`
  # (older releases used `max_temp`/`min_temp`); accept either.
  tmax_col <- intersect(c("air_tmax", "max_temp"), names(dt))[1]
  tmin_col <- intersect(c("air_tmin", "min_temp"), names(dt))[1]
  if (is.na(tmax_col) || is.na(tmin_col) || !"date" %in% names(dt)) {
    cli::cli_abort(
      "SILO response is missing expected columns ({.field date}, \\
       {.field air_tmax}/{.field max_temp}, {.field air_tmin}/{.field min_temp})."
    )
  }

  # weatherOz returns `date` as Date (or YYYYMMDD); parse defensively.
  raw_date <- dt[["date"]]
  dts <- suppressWarnings(as.Date(raw_date))
  if (anyNA(dts)) {
    dts <- as.Date(as.character(raw_date), format = "%Y%m%d")
  }

  out <- tibble::tibble(
    datetime        = dts,
    air_temp_mean_C = (as.numeric(dt[[tmax_col]]) + as.numeric(dt[[tmin_col]])) / 2
  )

  if (cache && !is.null(cache_path)) qs2::qs_save(out, cache_path)
  out
}


# ── SILO rainfall lookup ──────────────────────────────────────────────────────

#' Fetch daily rainfall from SILO for an Australian location
#'
#' Retrieves daily rainfall from the SILO Data Drill (a ~5 km gridded,
#' spatially interpolated climate surface covering Australia, 1889–present) for
#' a single latitude/longitude.  The result is the input hydrology series
#' accepted by [fit_reference_model()] when no gauge record is available.
#'
#' This function is the rainfall sibling of [get_silo_air_temp()] and uses the
#' same API, cache, and key mechanism.
#'
#' @section Attribution:
#' SILO data are © State of Queensland (Department of Environment, Science and
#' Innovation) and released under CC-BY 4.0. Cite SILO when you publish results
#' derived from this function.
#'
#' @section API key:
#' SILO requires an API key, which is simply your email address. By default it
#' is auto-detected via [weatherOz::get_key()] (from `.Renviron`/`.Rprofile`).
#'
#' @param latitude,longitude Numeric decimal-degree coordinates of the point of
#'   interest. Must fall within the SILO grid (approximately latitude -44 to
#'   -10, longitude 112 to 154). Snapped to the 0.05° grid by SILO.
#' @param start_date,end_date Start and end of the (inclusive) date range, as
#'   `Date` objects or `"YYYY-MM-DD"` strings.
#' @param api_key Character SILO API key (your email address). Default `NULL`
#'   defers to [weatherOz::get_key()].
#' @param cache Logical; cache the result on disk under
#'   `tools::R_user_dir("leachatetools", "cache")/silo`. Default `TRUE`.
#' @param refresh Logical; if `TRUE`, ignore and overwrite any cached result.
#'   Default `FALSE`.
#'
#' @return A tibble with one row per day:
#'   - `date` (`Date`)
#'   - `rainfall_mm` (numeric, mm/day)
#'   Ready to pass as `hydro` to [fit_reference_model()] with
#'   `hydro_type = "rainfall"`.
#'
#' @seealso [fit_reference_model()], [get_silo_air_temp()]
#'
#' @examples
#' \dontrun{
#' rain <- get_silo_rainfall(
#'   latitude   = -33.87, longitude  = 151.21,
#'   start_date = "2020-01-01", end_date = "2023-12-31"
#' )
#' }
#' @export
get_silo_rainfall <- function(
    latitude,
    longitude,
    start_date,
    end_date,
    api_key = NULL,
    cache   = TRUE,
    refresh = FALSE
) {
  if (!requireNamespace("weatherOz", quietly = TRUE)) {
    cli::cli_abort(c(
      "The {.pkg weatherOz} package is required to fetch SILO climate data.",
      "i" = "Install it with {.run install.packages(\"weatherOz\")}."
    ))
  }
  checkmate::assert_number(latitude,  lower = -44, upper = -10)
  checkmate::assert_number(longitude, lower = 112, upper = 154)
  checkmate::assert_flag(cache)
  checkmate::assert_flag(refresh)
  start_date <- as.Date(start_date)
  end_date   <- as.Date(end_date)
  if (is.na(start_date) || is.na(end_date)) {
    cli::cli_abort("{.arg start_date} and {.arg end_date} must be valid dates.")
  }
  if (end_date < start_date) {
    cli::cli_abort("{.arg end_date} ({end_date}) is before {.arg start_date} ({start_date}).")
  }

  cache_path <- NULL
  if (cache) {
    cache_dir <- file.path(tools::R_user_dir("leachatetools", "cache"), "silo")
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    key <- sprintf(
      "silo_rain_%+07.2f_%+07.2f_%s_%s.qs",
      latitude, longitude, format(start_date, "%Y%m%d"), format(end_date, "%Y%m%d")
    )
    cache_path <- file.path(cache_dir, key)
    if (!refresh && file.exists(cache_path)) {
      return(qs2::qs_read(cache_path))
    }
  }

  if (is.null(api_key)) api_key <- weatherOz::get_key(service = "SILO")

  dt <- weatherOz::get_data_drill(
    longitude  = longitude,
    latitude   = latitude,
    start_date = start_date,
    end_date   = end_date,
    values     = "rain",
    api_key    = api_key
  )

  if (!"rainfall" %in% names(dt) || !"date" %in% names(dt)) {
    cli::cli_abort(
      "SILO response is missing expected columns ({.field date}, {.field rainfall})."
    )
  }

  raw_date <- dt[["date"]]
  dts <- suppressWarnings(as.Date(raw_date))
  if (anyNA(dts)) dts <- as.Date(as.character(raw_date), format = "%Y%m%d")

  out <- tibble::tibble(
    date         = dts,
    rainfall_mm  = as.numeric(dt[["rainfall"]])
  )

  if (cache && !is.null(cache_path)) qs2::qs_save(out, cache_path)
  out
}
