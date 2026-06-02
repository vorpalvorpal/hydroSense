# Compute AmsPAF for a single sample

Thin wrapper over the shared AmsPAF helpers
([`.amspaf_adjust()`](https://www.kedumba.com.au/leachatetools/reference/dot-amspaf_adjust.md),
[`.amspaf_add_paf()`](https://www.kedumba.com.au/leachatetools/reference/dot-amspaf_add_paf.md),
[`.amspaf_combine()`](https://www.kedumba.com.au/leachatetools/reference/dot-amspaf_combine.md))
that processes one sample in isolation. Retained so the normalisation /
ARA / CA / IA steps can be driven (and unit-tested) for a single sample;
the batched
[`compute_amspaf_per_sample()`](https://www.kedumba.com.au/leachatetools/reference/compute_amspaf_per_sample.md)
path uses the same helpers so behaviour is guaranteed identical.

## Usage

``` r
compute_amspaf_one_sample(
  sample_rows,
  ref_table,
  ssd_params,
  min_analytes,
  method,
  guideline_dir,
  has_imputed = FALSE,
  ara_enabled = TRUE
)
```

## Arguments

- sample_rows:

  Long-format chemistry rows for one sample (one row per analyte; may
  include co-analyte rows used for normalisation). Must carry a
  `detected` column.

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

- has_imputed:

  Logical; whether the input carried an `imputed` column (controls
  `n_analytes_imputed` accounting).

- ara_enabled:

  Logical; whether the caller supplied a reference. Controls the
  `ref_source` diagnostic (see
  [`.amspaf_adjust()`](https://www.kedumba.com.au/leachatetools/reference/dot-amspaf_adjust.md)).

## Value

A one-row tibble (or zero-row tibble if the sample fails `min_analytes`)
with columns `value`, `n_analytes_used`, `n_analytes_imputed`,
`dominant_analyte`, `max_paf`, `analyte_pafs`, `dropped_analytes`.
