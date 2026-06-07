# Fetch daily rainfall from SILO for an Australian location

Retrieves daily rainfall from the SILO Data Drill (a ~5 km gridded,
spatially interpolated climate surface covering Australia, 1889–present)
for a single latitude/longitude. The result is the input hydrology
series accepted by
[`fit_reference_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_reference_model.md)
when no gauge record is available.

## Usage

``` r
get_silo_rainfall(
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
  `tools::R_user_dir("leachatetools", "cache")/silo`. Default `TRUE`.

- refresh:

  Logical; if `TRUE`, ignore and overwrite any cached result. Default
  `FALSE`.

## Value

A tibble with one row per day:

- `date` (`Date`)

- `rainfall_mm` (numeric, mm/day) Ready to pass as `hydro` to
  [`fit_reference_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_reference_model.md)
  with `hydro_type = "rainfall"`.

## Details

This function is the rainfall sibling of
[`get_silo_air_temp()`](https://vorpalvorpal.github.io/leachatetools/reference/get_silo_air_temp.md)
and uses the same API, cache, and key mechanism.

## Attribution

SILO data are © State of Queensland (Department of Environment, Science
and Innovation) and released under CC-BY 4.0. Cite SILO when you publish
results derived from this function.

## API key

SILO requires an API key, which is simply your email address. By default
it is auto-detected via
[`weatherOz::get_key()`](https://docs.ropensci.org/weatherOz/reference/get_key.html)
(from `.Renviron`/`.Rprofile`).

## See also

[`fit_reference_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_reference_model.md),
[`get_silo_air_temp()`](https://vorpalvorpal.github.io/leachatetools/reference/get_silo_air_temp.md)

## Examples

``` r
if (FALSE) { # \dontrun{
rain <- get_silo_rainfall(
  latitude   = -33.87, longitude  = 151.21,
  start_date = "2020-01-01", end_date = "2023-12-31"
)
} # }
```
