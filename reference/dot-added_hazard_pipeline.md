# Run the chronic added-hazard msPAF

Engine of `mspaf_pipeline(ara_form = "added_hazard")`; arguments are
validated by
[`.ah_validate()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-ah_validate.md)
first. Per target and day it forms the added hazard
`D = max((msT - msR) / (1 - msR), 0)` and the total `msT0` (the target
alone, below-detection analyte-days at concentration 0, on the same
analytes as `D`); both are time-weighted over the days where `D` exists.

## Usage

``` r
.added_hazard_pipeline(
  target,
  reference,
  daily_args,
  chronic,
  focal_dates,
  focal_by,
  tau,
  window,
  also_total,
  contributions,
  min_coverage
)
```

## Value

See
[`mspaf_pipeline()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_pipeline.md).
