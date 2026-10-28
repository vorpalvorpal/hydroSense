# Re-key an add_mspaf() contributions table from synthetic sample to date

Re-key an add_mspaf() contributions table from synthetic sample to date

## Usage

``` r
.daily_contributions(ctb, id_date_map)
```

## Arguments

- ctb:

  The `"contributions"` attribute of an
  [`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)
  result.

- id_date_map:

  `(sample_id, site_id, .date)` lookup.

## Value

Tibble keyed by `site_id`, `date`, ordered by site, date, analyte.
