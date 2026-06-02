# Combine adjusted+PAF'd tox rows for one sample into a single AmsPAF row

Performs Concentration Addition within each mode-of-action group and
Independent Action across groups, then assembles the one-row diagnostic
tibble. Operates on a *single* sample's rows.

## Usage

``` r
.amspaf_combine(tox, has_imputed = FALSE)
```

## Arguments

- tox:

  Rows for one sample from
  [`.amspaf_add_paf()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-amspaf_add_paf.md)
  (needs `C_adj`, `hc50`, `sigma`, `moa_group`, `PAF`; `analyte`,
  `ref_source` carried into the `analyte_pafs` diagnostic if present).

- has_imputed:

  Logical; whether the input carried an `imputed` column.

## Value

A one-row tibble: `value`, `n_analytes_used`, `n_analytes_imputed`,
`dominant_analyte`, `max_paf`, `analyte_pafs` (does NOT add `sample_id`
or `dropped_analytes` — the caller attaches those).
