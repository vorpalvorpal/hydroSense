# Build synthetic long-format daily samples from interpolated chemistry

Assigns `sample_id = "daily_{YYYY-MM-DD}_{site}"` per day. Keeps `.date`
as a column (the caller extracts it before passing to
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md)).
No `focal_date` column is added – this is deliberate so
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md)
uses the instant (pointwise) ARA path, not the chronic integrated path.

## Usage

``` r
.build_synthetic_samples(daily_long, site)
```

## Arguments

- daily_long:

  Output of
  [`.build_daily_chem()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-build_daily_chem.md)
  (after temperature fill).

- site:

  Site identifier string.

## Value

Long-format tibble ready for
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md)
(after removing `.date`).
