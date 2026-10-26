# Chronic aggregation for one point over all of its focal dates

Chronic aggregation for one point over all of its focal dates

## Usage

``` r
.twa_site(
  x,
  f,
  r_idx,
  tau_days,
  window_days,
  summary,
  anchor_outside_window,
  eps,
  draws
)
```

## Arguments

- x:

  List of the point's row vectors: `sample_id`, `date` (numeric days),
  `analyte`, `value`, `imputed` (logical), `draw_id` (or `NULL`).

- f:

  Numeric focal dates (days); `r_idx` their focal-grid rows.

- tau_days, window_days, summary, anchor_outside_window, eps:

  As in
  [`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md).

- draws:

  Active draw domain
  ([`.draw_domain()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-draw_domain.md));
  `integer(0)` for point input.

## Value

Data frame `(r, analyte, draw, value, n_in, n_imp)` or `NULL`.
