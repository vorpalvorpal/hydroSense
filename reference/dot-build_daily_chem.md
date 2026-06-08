# Interpolate per-analyte grab chemistry onto a daily date grid

For each analyte present in `site_rows`, produces one row per date in
`dates` using either forward-fill or linear interpolation. Adds `.date`
and `.measured` columns for downstream diagnostics.

## Usage

``` r
.build_daily_chem(site_rows, dates, interpolation, leading_edge, tox_analytes)
```

## Arguments

- site_rows:

  Long-format chemistry for one site.

- dates:

  Date vector giving the target daily grid.

- interpolation:

  `"forward_fill"` or `"linear"`.

- leading_edge:

  `"drop"` or `"backfill"`.

- tox_analytes:

  SSD-eligible analyte names (for log-space decision).

## Value

Long-format tibble with `.date`, `analyte`, `value`, `detected`,
`.measured`, and any pass-through columns (`units.analyte`, etc.).
