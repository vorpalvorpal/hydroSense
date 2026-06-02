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
# Monthly focal dates across 2024, e.g. to feed time_weighted_aggregate().
expand_focal_dates("2024-01-01", "2024-12-31", by = "month")
#>  [1] "2024-01-01" "2024-02-01" "2024-03-01" "2024-04-01" "2024-05-01"
#>  [6] "2024-06-01" "2024-07-01" "2024-08-01" "2024-09-01" "2024-10-01"
#> [11] "2024-11-01" "2024-12-01"
```
