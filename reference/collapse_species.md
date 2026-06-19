# Pivot and collapse transforming ion species to conserved totals

Internal helper used by
[`add_lmf`](https://vorpalvorpal.github.io/hydroSense/reference/add_lmf.md),
[`build_reference_endmember`](https://vorpalvorpal.github.io/hydroSense/reference/build_reference_endmember.md),
and
[`build_leachate_endmember`](https://vorpalvorpal.github.io/hydroSense/reference/build_leachate_endmember.md)
to avoid code duplication. Pivots long-format meq data to wide format
and replaces the individual N species (NH4-N, NO3-N, NO2-N) with
`total_N_` and the carbonate species (CO3, HCO3) with `total_alk_`.

## Usage

``` r
collapse_species(df_meq, id_cols)
```

## Arguments

- df_meq:

  Long-format dataframe containing at minimum columns `name.analyte` and
  `value`, plus any columns named in `id_cols`.

- id_cols:

  Character vector of column names to preserve as row identifiers during
  the wide pivot.

## Value

Wide-format tibble with one row per unique combination of `id_cols`
values. Individual N and carbonate species columns are replaced by
`total_N_` and `total_alk_` respectively.
