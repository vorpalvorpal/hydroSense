# Compute hydrology features (hydro_short, hydro_long) at target dates

Dispatches to
[`.compute_api()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-compute_api.md)
for rainfall or
[`.compute_antecedent_mean()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-compute_antecedent_mean.md)
for stage/discharge, at both the short and long tau values.

## Usage

``` r
.compute_hydro_features(
  hydro,
  target_dates,
  tau_short,
  tau_long,
  hydro_type = "rainfall"
)
```

## Arguments

- hydro:

  Data frame with columns `date` (Date) and `value` (numeric).

- target_dates:

  Date vector.

- tau_short, tau_long:

  Memory constants (days).

- hydro_type:

  `"rainfall"`, `"stage"`, or `"discharge"`.

## Value

Tibble `(date, hydro_short, hydro_long)`.
