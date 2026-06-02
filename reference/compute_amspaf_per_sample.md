# Compute AmsPAF for each sample in a per-feature data block

Internal workhorse called by
[`add_amspaf()`](https://www.kedumba.com.au/leachatetools/reference/add_amspaf.md)
for each `site_id` group. Runs in three phases so the (relatively
expensive) SSD evaluation is **batched across samples** rather than
called once per (sample × analyte):

1.  per-sample chemistry normalisation + ARA shift
    ([`.amspaf_adjust()`](https://www.kedumba.com.au/leachatetools/reference/dot-amspaf_adjust.md));

2.  one vectorised
    [`ssdtools::ssd_hp()`](https://bcgov.github.io/ssdtools/reference/ssd_hp.html)
    call per analyte across every sample
    ([`.amspaf_add_paf()`](https://www.kedumba.com.au/leachatetools/reference/dot-amspaf_add_paf.md));

3.  per-sample CA/IA mixture combination
    ([`.amspaf_combine()`](https://www.kedumba.com.au/leachatetools/reference/dot-amspaf_combine.md)).

## Usage

``` r
compute_amspaf_per_sample(
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
  [`derive_ssd_params()`](https://www.kedumba.com.au/leachatetools/reference/derive_ssd_params.md)
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
  [`.amspaf_adjust()`](https://www.kedumba.com.au/leachatetools/reference/dot-amspaf_adjust.md)).

## Value

Tibble with one row per sample that passes `min_analytes`, columns:
`sample_id`, `value`, `n_analytes_used`, `n_analytes_imputed`,
`dominant_analyte`, `max_paf`, `analyte_pafs`, `dropped_analytes`.
