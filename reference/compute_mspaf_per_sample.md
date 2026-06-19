# Compute msPAF for each sample in a per-feature data block

Internal workhorse called by
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)
for each `site_id` group. Runs in three phases so the (relatively
expensive) SSD evaluation is **batched across samples** rather than
called once per (sample × analyte):

1.  chemistry normalisation + ARA shift, vectorised across all rows
    ([`.mspaf_adjust()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-mspaf_adjust.md));

2.  one vectorised
    [`ssdtools::ssd_hp()`](https://bcgov.github.io/ssdtools/reference/ssd_hp.html)
    call per analyte across every sample
    ([`.mspaf_add_paf()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-mspaf_add_paf.md));

3.  CA/IA mixture combination via grouped reductions over (sample, draw,
    MOA group).

## Usage

``` r
compute_mspaf_per_sample(
  sample_data,
  ref_table,
  ssd_params,
  min_analytes,
  method,
  guideline_dir,
  ara_enabled = TRUE
)
```

## Arguments

- sample_data:

  Per-feature long-format df (may include co-analyte rows such as pH,
  DOC alongside toxicant rows).

- ref_table:

  Tibble `(analyte, ref_norm)` from `prep_ref$ref_table`.

- ssd_params:

  Tibble from
  [`derive_ssd_params()`](https://vorpalvorpal.github.io/hydroSense/reference/derive_ssd_params.md)
  (incl. the `fit` column).

- min_analytes:

  Minimum analytes required.

- method:

  SSD method (used only for the rare NULL-fit fallback).

- guideline_dir:

  Path to ANZG XLSX folder (NULL-fit fallback only).

- ara_enabled:

  Logical; whether the caller supplied a reference. Controls the
  `ref_source` diagnostic (see
  [`.mspaf_adjust()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-mspaf_adjust.md)).

## Value

Tibble with one row per sample that passes `min_analytes`, columns:
`sample_id`, `value`, `n_analytes_used`, `n_analytes_imputed`,
`dominant_analyte`, `max_paf`, `analyte_pafs`, `dropped_analytes`.
