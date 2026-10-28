# Aggregate daily msPAF contributions into chronic contributions

Each analyte's daily contribution series (zero on days it did not enter
that day's msPAF) is passed through
[`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md)
with the same arguments, focal grid and daily schedule as the chronic
msPAF itself, so the chronic contributions sum to the chronic value.

## Usage

``` r
.chronic_contributions(
  daily,
  focal,
  tau,
  window,
  chronic_summary,
  also_total = FALSE
)
```

## Arguments

- daily:

  [`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)
  result carrying a `"contributions"` attribute.

## Value

Tibble `(site_id, focal_date, analyte, group, contribution)`, plus a
leading `variant` column when `also_total`.
