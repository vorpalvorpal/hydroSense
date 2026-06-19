# Summarise a draw-carrier frame to posterior median and credible interval

Collapses per-cell posterior draws to a point estimate and credible
interval. Returns `df` unchanged (identity) when the input carries no
draws (point frame), preserving the degradation guarantee — no CI
columns appear for callers that never produce draws.

## Usage

``` r
summarise_draws(df, interval = 0.9, central = c("median", "mean"))
```

## Arguments

- df:

  Long-format data frame following the draw-carrier contract.

- interval:

  Width of the credible interval. Default `0.90` yields 5th / 95th
  percentile bounds.

- central:

  Central-tendency statistic: `"median"` (default) or `"mean"`.

## Value

`df` collapsed to one row per cell, with columns `value` (central
estimate), `value_lower`, `value_upper` (interval bounds), and
`n_draws`. Returns `df` unchanged when the input carries no draws.

## Details

**Per-draw diagnostic columns** (`draw_id`, `dominant_analyte`,
`max_paf`, `analyte_pafs`) are dropped in draws mode; they are available
by passing `return = "draws"` to
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)
or
[`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md).

Exact cells (`draw_id = NA`) each collapse to a degenerate interval:
`value_lower = value_upper = value`, `n_draws = 1`.

## See also

[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md),
[`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Collapse msPAF draws to median + 90 % CI
mspaf_draws |>
  add_mspaf(return = "draws") |>
  summarise_draws(interval = 0.90)
} # }
```
