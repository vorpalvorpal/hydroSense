# Estimate water temperature from air temperature measurements

Fits a simple linear regression between observed air temperature and
observed water temperature, then predicts water temperature for a set of
target dates. The result can be added to a chemistry data frame as
`analyte = "temperature"` rows, enabling the NH₃-N pH/temperature
normalisation in
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md).

## Usage

``` r
estimate_water_temp(
  air_temp_df,
  water_temp_obs = NULL,
  target_dates = NULL,
  lag_days = 0L,
  site_id = NA_character_,
  seasonal = c("auto", "off", "on"),
  seasonal_min_n = 8L,
  seasonal_min_quarters = 3L,
  fallback = NULL,
  select = "auto",
  auto_min_n = 20L,
  water_body_type = NULL
)
```

## Arguments

- air_temp_df:

  Data frame of air temperature observations. Required columns:

  - `datetime` (Date or POSIXct) — measurement date

  - `air_temp_mean_C` (numeric) — mean daily air temperature in Celsius

- water_temp_obs:

  Data frame of observed water temperatures, used to calibrate the
  air-water regression. Required columns:

  - `datetime` (Date or POSIXct)

  - `water_temp_C` (numeric) — observed water temperature in Celsius At
    least 10 paired observations are recommended; fewer than 5 triggers
    an error unless `fallback` is also supplied. `NULL` is accepted when
    `fallback` is provided.

- target_dates:

  Date vector of dates for which to predict water temperature. Defaults
  to all dates present in `air_temp_df`. Dates with no air temperature
  observation are excluded with a warning.

- lag_days:

  Integer. Apply a lag to air temperature before regressing: the water
  temperature on day `t` is predicted from the air temperature on day
  `t - lag_days`. A lag of 1–3 days is appropriate for streams; 0
  (default) is appropriate for ponds or when using aggregated data.

- site_id:

  Character. Value to use for the `site_id` column in the returned
  chemistry rows. Default `NA_character_`.

- seasonal:

  One of `"auto"` (default), `"off"`, or `"on"`. Controls whether a
  day-of-year seasonal term is added to the air-temperature regression
  (see *Model selection*). `"auto"` fits both the air-only and air +
  season models and keeps whichever has the lower AICc; `"off"` forces
  the air-only model (the legacy behaviour); `"on"` forces the seasonal
  model whenever it is eligible.

- seasonal_min_n:

  Integer. Minimum number of paired observations before the seasonal
  model is even considered (default 8). The seasonal model adds two
  parameters, so a small buffer above that is needed for AICc to behave.

- seasonal_min_quarters:

  Integer 1-4. Minimum number of distinct calendar quarters the training
  observations must span before the seasonal model is considered
  (default 3). Calibration data confined to one or two seasons cannot
  anchor an annual cycle and would extrapolate badly.

- fallback:

  Air-to-water relationship used when there are too few paired
  observations to fit a regression (fewer than 5 rows in
  `water_temp_obs`, or `water_temp_obs = NULL`). Accepted forms:

  - `"identity"` — `water = air` (simplest default for small streams).

  - `c(intercept, slope)` — `water = intercept + slope * air`.

  - A function `function(air) ...` — arbitrary transformation. `NULL`
    (default) preserves the existing behaviour: fewer than 5 pairs is a
    hard error.

- select:

  Antecedence window for the air-temperature regressor, in days. The
  regressor is the trailing `select`-day mean of air temperature.

  - `"auto"` (default) — select the window by leave-one-out
    cross-validation (PRESS statistic, exact for OLS) with the **1-SE
    rule**: the shortest window whose CV RMSE is within one SE of the
    best is chosen, to avoid over-fitting to noise. Requires at least
    `auto_min_n` pairs; falls back to 1 (or to the `water_body_type`
    default) when pairs are insufficient.

  - A positive integer — use that fixed window (1 = daily air temp,
    equivalent to the original behaviour).

- auto_min_n:

  Integer. Minimum number of paired observations required before the
  `select = "auto"` window search is attempted (default 20). Below this
  threshold the window falls back to 1 (or `water_body_type`).

- water_body_type:

  Character. Physical type of the water body, used as the fallback
  window when `select = "auto"` but `n_pairs < auto_min_n`. One of
  `"stream_small"` (1 d), `"stream_large"` (5 d), `"estuary"` (7 d),
  `"pond"` (10 d), `"lake"` (30 d), `"reservoir"` (60 d),
  `"groundwater"` (300 d). `NULL` (default) falls back to window = 1.
  Ignored when the CV search runs.

## Value

A tibble suitable for binding onto a chemistry data frame, with columns:

- `datetime` (Date)

- `analyte` — `"temperature"`

- `value` — predicted water temperature (°C)

- `detected` — `TRUE`

- `site_id` — from `site_id` argument

- `sample_id` — `NA_character_` (set by caller if needed) Attributes
  attached for inspection:

- `attr(result, "lm_fit")` — the selected fitted `lm` object (`NULL` for
  fallback).

- `attr(result, "model")` — label of the selected model (`"fallback"`
  when no pairs were available).

- `attr(result, "seasonal_used")` — `TRUE` if the seasonal model won.

- `attr(result, "model_comparison")` — data frame of AICc / R² per
  candidate model and which was selected.

- `attr(result, "window_selected")` — chosen antecedence window in days
  (`NA` for fallback path).

- `attr(result, "window_comparison")` — data frame of LOO-CV RMSE per
  candidate window (`NULL` when window was fixed or search was skipped).

## Air temperature variable

This function requires **mean daily air temperature** (°C). Mean daily
air temperature is the standard predictor for daily mean water
temperature in the literature and has the best empirical performance for
streams and ponds. If your weather station or BOM data only provides
daily maximum and minimum, compute the mean as `(Tmax + Tmin) / 2`.

## Why water temperature matters

The ANZG NH₃-N guideline values are derived at a specific reference
condition (pH 7.0, 20 °C total ammonia-N). The normalisation formula
adjusts measured NH₃-N to these conditions before SSD lookup — and since
the un-ionised fraction changes steeply with temperature and pH (roughly
doubling for each 5 °C increase near 20 °C), using the wrong temperature
has a large effect on the normalised concentration. A default value is
not provided: you must supply site-appropriate temperature estimates.

## Model selection

**Window selection** (antecedence timescale): water bodies with more
thermal mass integrate air temperature over longer periods. The function
searches log-spaced trailing windows from 1 to 300 days, fitting a
simple OLS for each, and selects by leave-one-out CV RMSE (PRESS
statistic — closed form, no refitting loop). The **1-SE rule** then
chooses the shortest window within one standard error of the best,
avoiding spurious preference for a long window that barely beats shorter
alternatives. All window candidates have the same number of parameters,
so AICc cannot distinguish them; CV is required.

**Seasonal term** (thermal hysteresis): water temperature lags air
temperature seasonally — at the same air temperature, water tends to be
cooler in spring and warmer in autumn. Adding a **first-harmonic
day-of-year term** — `sin(2*pi*doy/365.25)` and `cos(2*pi*doy/365.25)` —
lets the regression represent it. Under `seasonal = "auto"` both the
window-regressor- only and window-regressor + harmonics models are
compared by **AICc**; the lower-AICc model is used. AICc is appropriate
here because the two candidates differ in complexity (2 vs 4
parameters). The seasonal model is only eligible when there are at least
`seasonal_min_n` observations spanning at least `seasonal_min_quarters`
quarters.

## See also

[`get_silo_air_temp()`](https://vorpalvorpal.github.io/leachatetools/reference/get_silo_air_temp.md)
to source `air_temp_df` from SILO for an Australian location;
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md),
which requires the resulting water `temperature` rows for ammonia.

## Examples

``` r
set.seed(1)
air <- tibble::tibble(
  datetime        = seq(as.Date("2020-01-01"), as.Date("2022-12-31"), by = "day")
)
air$air_temp_mean_C <- rnorm(nrow(air), mean = 15, sd = 8)
wt_obs <- tibble::tibble(
  datetime     = sample(air$datetime, 80),
  water_temp_C = NA_real_
)
wt_obs$water_temp_C <-
  air$air_temp_mean_C[match(wt_obs$datetime, air$datetime)] * 0.85 + 2 +
  rnorm(80, 0, 1)
wt <- estimate_water_temp(air, wt_obs)
#> ℹ Selected antecedence window: 1 d (LOO-CV RMSE = 1.09 °C, n = 80).
#> ℹ Selected air-water model: air + season (day-of-year harmonic). R² = 0.976,
#>   RMSE = 1.01 °C, AICc = 239.4 (n = 80).
#>   AICc comparison (lower preferred): air only = 243.7, air + season = 239.4.
attr(wt, "model")  # which model was selected (air-only vs air + season)
#> [1] "air + season (day-of-year harmonic)"
```
