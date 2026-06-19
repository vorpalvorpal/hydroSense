# Fetch daily mean air temperature from SILO for an Australian location

Retrieves daily air temperature from the SILO Data Drill (a ~5 km
gridded, spatially interpolated climate surface covering Australia,
1889–present) for a single latitude/longitude, and returns **daily mean
air temperature** (`(Tmax + Tmin) / 2`, °C) in exactly the shape
[`estimate_water_temp()`](https://vorpalvorpal.github.io/hydroSense/reference/estimate_water_temp.md)
expects as its `air_temp_df`. The typical workflow is:

## Usage

``` r
get_silo_air_temp(
  latitude,
  longitude,
  start_date,
  end_date,
  api_key = NULL,
  cache = TRUE,
  refresh = FALSE
)
```

## Arguments

- latitude, longitude:

  Numeric decimal-degree coordinates of the point of interest. Must fall
  within the SILO grid (approximately latitude -44 to -10, longitude 112
  to 154). Snapped to the 0.05° grid by SILO.

- start_date, end_date:

  Start and end of the (inclusive) date range, as `Date` objects or
  `"YYYY-MM-DD"` strings.

- api_key:

  Character SILO API key (your email address). Default `NULL` defers to
  [`weatherOz::get_key()`](https://docs.ropensci.org/weatherOz/reference/get_key.html).

- cache:

  Logical; cache the result on disk under
  `tools::R_user_dir("hydroSense", "cache")/silo`. Default `TRUE`.

- refresh:

  Logical; if `TRUE`, ignore and overwrite any cached result. Default
  `FALSE`.

## Value

A tibble with one row per day:

- `datetime` (`Date`)

- `air_temp_mean_C` (numeric, °C) Ready to pass as `air_temp_df` to
  [`estimate_water_temp()`](https://vorpalvorpal.github.io/hydroSense/reference/estimate_water_temp.md).

## Details

    air   <- get_silo_air_temp(lat, lon, start, end)   # SILO mean air temp
    wt    <- estimate_water_temp(air, water_temp_obs)  # calibrate water = f(air)
    chem  <- dplyr::bind_rows(chem, dplyr::mutate(wt, sample_id = ...))

This wraps
[`weatherOz::get_data_drill()`](https://docs.ropensci.org/weatherOz/reference/get_data_drill.html);
the **weatherOz** package must be installed (it is listed under
`Suggests`). Results are cached on disk so repeat calls for the same
grid cell and date range do not re-hit the API.

## Attribution

SILO data are © State of Queensland (Department of Environment, Science
and Innovation) and released under CC-BY 4.0. Cite SILO when you publish
results derived from this function.

## API key

SILO requires an API key, which is simply your email address. By default
it is auto-detected via
[`weatherOz::get_key()`](https://docs.ropensci.org/weatherOz/reference/get_key.html)
(from `.Renviron`/`.Rprofile`); see the weatherOz documentation for
one-time setup, or pass `api_key` directly.

## See also

[`estimate_water_temp()`](https://vorpalvorpal.github.io/hydroSense/reference/estimate_water_temp.md)

## Examples

``` r
if (FALSE) { # \dontrun{
air <- get_silo_air_temp(
  latitude = -33.87, longitude = 151.21,
  start_date = "2020-01-01", end_date = "2023-12-31"
)
} # }
```
