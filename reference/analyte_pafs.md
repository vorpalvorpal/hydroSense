# Per-analyte PAF breakdown from add_amspaf()

After calling
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md),
this accessor returns the per-analyte contribution breakdown behind each
AmsPAF value: the ARA-adjusted concentration, single-substance PAF, MOA
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
  [`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md).

## Value

A tibble with columns `site_id`, `sample_id`, `draw_id` (draws mode
only), `analyte`, `C_adj`, `PAF`, `moa_group`, `ref_source`. Returns
`NULL` (with a message) if the attribute is absent.

## Details

The attribute is stored by
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md)
and is dropped by most dplyr verbs, so read it before further wrangling.

## See also

[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md),
[`ara_summary()`](https://vorpalvorpal.github.io/leachatetools/reference/ara_summary.md)
