# Pivot and collapse per-ion censored (non-detect) flags to match [`collapse_species`](https://vorpalvorpal.github.io/hydroSense/reference/collapse_species.md)'s ion panel

Internal helper used by
[`add_lmf`](https://vorpalvorpal.github.io/hydroSense/reference/add_lmf.md)
to build the censored-flag side table consumed by
[`compute_lmf_for_sample`](https://vorpalvorpal.github.io/hydroSense/reference/compute_lmf_for_sample.md)
(Issue 5: give non-detects detection-limit-scale uncertainty instead of
a fixed RSD). A collapsed total (`total_N_`, `total_alk_`) is censored
only if every present constituent species was a non-detect.

## Usage

``` r
collapse_censored(df_meq, id_cols)
```

## Arguments

- df_meq:

  Long-format dataframe with columns `analyte`, `detected`, plus any
  columns named in `id_cols`.

- id_cols:

  Character vector of column names to preserve as row identifiers during
  the wide pivot.

## Value

Wide-format tibble of logicals, one row per unique combination of
`id_cols`, one column per collapsed panel ion.
