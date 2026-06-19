# Per-analyte PAF breakdown from add_mspaf()

After calling
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md),
this accessor returns the per-analyte contribution breakdown behind each
msPAF value: the ARA-adjusted concentration, single-substance PAF, MOA
group and reference source for every assessed (sample × draw × analyte).
It replaces the former `analyte_pafs` list-column with a flat, tidy
frame (filter/join directly).

## Usage

``` r
analyte_pafs(x)
```

## Arguments

- x:

  A data frame returned by
  [`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md).

## Value

A tibble with columns `site_id`, `sample_id`, `draw_id` (draws mode
only), `analyte`, `C_adj`, `PAF`, `moa_group`, `ref_source`. Returns
`NULL` (with a message) if the attribute is absent.

## Details

The attribute is stored by
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)
and is dropped by most dplyr verbs, so read it before further wrangling.

## See also

[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md),
[`ara_summary()`](https://vorpalvorpal.github.io/hydroSense/reference/ara_summary.md)
