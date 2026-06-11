# Align a residual path (values on the smoother's clipped grid) to query dates

Returns `NA` for query dates outside the analyte's clipped grab span –
those rows are dropped by
[`.resolve_target_impact()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-resolve_target_impact.md)
(per-analyte clipping).

## Usage

``` r
.residual_on_qdates(grid_dates, path, qdates)
```
