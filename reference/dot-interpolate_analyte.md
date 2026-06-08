# Interpolate one analyte onto a target date vector

Core per-analyte interpolation. Returns a tibble with `.date`, `value`,
`detected`, `.measured`.

## Usage

``` r
.interpolate_analyte(
  obs_dates,
  obs_values,
  obs_detected,
  target_dates,
  interpolation,
  leading_edge,
  log_space = FALSE
)
```
