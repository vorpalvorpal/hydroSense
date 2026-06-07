# Antecedent Precipitation Index (API) for rainfall

Computes `API(t) = sum_{i=0}^{window} P(t-i) * k^i` where the decay
constant `k = exp(-1 / window_days)` — meaning rainfall `window_days`
days ago carries weight `1/e ≈ 0.37`. Returns 0 for target dates with no
hydro record in the preceding window.

## Usage

``` r
.compute_api(hydro_values, hydro_dates, target_dates, window_days)
```

## Arguments

- hydro_values:

  Numeric vector of daily rainfall values.

- hydro_dates:

  Date vector matching `hydro_values`.

- target_dates:

  Date vector; API is evaluated at each.

- window_days:

  Integer memory length (days).

## Value

Numeric vector, same length as `target_dates`.
