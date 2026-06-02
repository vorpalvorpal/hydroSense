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
  water_temp_obs,
  target_dates = NULL,
  lag_days = 0L,
  site_id = NA_character_
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
    least 10 paired observations are recommended. Cannot be `NULL` — a
    regression requires training data. If you have no water temperature
    observations at all, collect a season of spot measurements alongside
    air temperatures before proceeding.

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

## Value

A tibble suitable for binding onto a chemistry data frame, with columns:

- `datetime` (Date)

- `analyte` — `"temperature"`

- `value` — predicted water temperature (°C)

- `detected` — `TRUE`

- `site_id` — from `site_id` argument

- `sample_id` — `NA_character_` (set by caller if needed) Attach
  `attr(result, "lm_fit")` — the fitted `lm` object for inspection.

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

## See also

[`get_silo_air_temp()`](https://vorpalvorpal.github.io/leachatetools/reference/get_silo_air_temp.md)
to source `air_temp_df` from SILO for an Australian location;
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md),
which requires the resulting water `temperature` rows for ammonia.

## Examples

``` r
if (FALSE) { # \dontrun{
air  <- tibble::tibble(
  datetime       = seq(as.Date("2020-01-01"), as.Date("2023-12-31"), by = "day"),
  air_temp_mean_C = rnorm(length(datetime), mean = 15, sd = 8)
)
wt_obs <- tibble::tibble(
  datetime      = sample(air$datetime, 80),
  water_temp_C  = air$air_temp_mean_C[match(datetime, air$datetime)] * 0.85 + 2 + rnorm(80, 0, 1)
)
wt <- estimate_water_temp(air, wt_obs)
# Add to chemistry df:
# chemistry <- dplyr::bind_rows(chemistry, dplyr::mutate(wt, sample_id = ...))
} # }
```
