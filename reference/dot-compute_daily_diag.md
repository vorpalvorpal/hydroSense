# Compute per-day diagnostics: n_measured_analytes, days_since_last_sample

Operates on SSD-eligible rows only (toxicants drive the AmsPAF;
co-analyte sampling frequency is generally higher and not the
bottleneck).

## Usage

``` r
.compute_daily_diag(daily_long, tox_analytes, site)
```

## Arguments

- daily_long:

  Output of
  [`.build_daily_chem()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-build_daily_chem.md)
  (possibly augmented).

- tox_analytes:

  SSD-eligible analyte names.

- site:

  Site identifier for the output `site_id` column.

## Value

Tibble `(date, site_id, n_measured_analytes, days_since_last_sample)`.
