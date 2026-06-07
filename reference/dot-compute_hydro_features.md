# Compute hydrology features (hydro_short, hydro_long) at target dates

Dispatches to
[`.compute_api()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-compute_api.md)
for rainfall or
[`.compute_antecedent_mean()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-compute_antecedent_mean.md)
for stage/discharge, at both the short and long window lengths.

## Usage

``` r
.compute_hydro_features(
  hydro,
  target_dates,
  window_short,
  window_long,
  hydro_type = "rainfall"
)
```

## Arguments

- hydro:

  Data frame with columns `date` (Date) and `value` (numeric).

- target_dates:

  Date vector.

- window_short, window_long:

  Memory lengths (days).

- hydro_type:

  `"rainfall"`, `"stage"`, or `"discharge"`.

## Value

Tibble `(date, hydro_short, hydro_long)`.
