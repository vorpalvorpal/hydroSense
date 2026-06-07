# Antecedent mean for stage or discharge

Rolling mean over `[t - window_days, t]`. Returns `NA` for target dates
with no hydro record in the window.

## Usage

``` r
.compute_antecedent_mean(hydro_values, hydro_dates, target_dates, window_days)
```

## Arguments

- hydro_values:

  Numeric vector of daily stage/discharge values.

- hydro_dates:

  Date vector matching `hydro_values`.

- target_dates:

  Date vector; mean is evaluated at each.

- window_days:

  Integer memory length (days).

## Value

Numeric vector, same length as `target_dates`.
