# Per-analyte and pooled LOO coverage over a fitted target model

Per-analyte and pooled LOO coverage over a fitted target model

## Usage

``` r
.loo_anchor_coverage(
  target_model,
  interval = 0.9,
  block = 1L,
  n_fit_min = 8L,
  scale = 1
)
```

## Arguments

- target_model:

  A `target_model` (uses each model's `$anchors$S` / `$d_anchors$S` and
  `$tier`).

- interval, block, n_fit_min, scale:

  Passed to
  [`.loo_coverage_series()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-loo_coverage_series.md).

## Value

A tibble: one row per analyte (`analyte, tier, n, coverage, mean_width`)
plus a final `(pooled)` row (n-weighted coverage, mean width).
