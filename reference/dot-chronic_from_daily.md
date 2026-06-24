# Aggregate a daily msPAF frame into chronic time-weighted msPAF

Reshapes the
[`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)
output into the long `analyte = "msPAF"` schema
[`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md)
requires, then aggregates. In point mode the daily `mspaf` line is
smeared into a single chronic value per (focal date × site); in draws
mode each envelope column (`mspaf_ignorable` / `mspaf_informative`) is
propagated per draw and collapsed with
[`.summarise_bracket()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-summarise_bracket.md).

## Usage

``` r
.chronic_from_daily(
  daily,
  draws_mode,
  gap_uncertainty,
  focal,
  tau,
  window,
  chronic_summary,
  interval,
  central
)
```
