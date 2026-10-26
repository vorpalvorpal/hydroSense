# Vectorised chronic aggregation over a focal grid

Vectorised chronic aggregation over a focal grid

## Usage

``` r
.twa_compute(
  df,
  focal_grid,
  tau_days,
  window_days,
  summary,
  anchor_outside_window,
  eps,
  has_imputed,
  draws
)
```

## Arguments

- df:

  Long chemistry with a `.date` column (see
  [`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md)).

- focal_grid:

  Data frame `(focal_date, site_id)`; one output block per row, in row
  order.

- tau_days, window_days, summary, anchor_outside_window, eps:

  As in
  [`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md).

- has_imputed:

  Logical; `df` carries an `imputed` column.

- draws:

  Active draw domain
  ([`.draw_domain()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-draw_domain.md));
  `integer(0)` for point input.

## Value

`NULL` when nothing could be aggregated, else a data frame
`(r, analyte, draw, value, n_in, n_imp)` ordered by focal-grid row `r`,
then analyte (C locale), then draw.
