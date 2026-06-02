# Generate a sequence of focal dates for chronic AmsPAF computation

A thin convenience wrapper around
[`base::seq.Date()`](https://rdrr.io/r/base/seq.Date.html) for
generating the `focal_dates` vector passed to
[`time_weighted_aggregate()`](https://vorpalvorpal.github.io/leachatetools/reference/time_weighted_aggregate.md).
The most common use is a daily sequence for time-series analysis.

## Usage

``` r
expand_focal_dates(start, end, by = "day")
```

## Arguments

- start:

  Start date (character `"YYYY-MM-DD"` or `Date`).

- end:

  End date (character `"YYYY-MM-DD"` or `Date`).

- by:

  Increment. Passed to
  [`base::seq.Date()`](https://rdrr.io/r/base/seq.Date.html). Common
  values: `"day"`, `"week"`, `"month"`. Default `"day"`.

## Value

A `Date` vector from `start` to `end` at the specified increment.

## Examples

``` r
if (FALSE) { # \dontrun{
# Daily sequence for 2024–2025
focal_dates <- expand_focal_dates("2024-01-01", "2025-12-31", by = "day")

chr_chem <- time_weighted_aggregate(imp, focal_dates = focal_dates)
} # }
```
