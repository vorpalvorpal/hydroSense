# Build the model data frame for one analyte

Computes hydro features at the observation dates and assembles
`(y, doy, hydro_short, hydro_long)` where `y = log(value_norm + eps)`.

## Usage

``` r
.build_analyte_model_df(
  obs,
  hydro,
  tau_short,
  tau_long,
  hydro_type,
  eps = 1e-09
)
```

## Arguments

- obs:

  Tibble `(date, value_norm)` for one analyte.

- hydro:

  Daily hydro series.

- tau_short, tau_long:

  Reservoir memory constants (days).

- hydro_type:

  Hydro type string.

- eps:

  Guard for log transform.

## Value

Tibble `(y, doy, hydro_short, hydro_long)` sorted by date.
