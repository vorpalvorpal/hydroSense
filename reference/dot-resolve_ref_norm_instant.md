# Resolve ref_norm for instantaneous target samples (one date per sample)

Resolve ref_norm for instantaneous target samples (one date per sample)

## Usage

``` r
.resolve_ref_norm_instant(ref_model, target_df)
```

## Arguments

- ref_model:

  A `reference_model` object.

- target_df:

  Tibble with columns `sample_id`, `datetime` (Date).

## Value

Tibble `(sample_id, analyte, ref_norm, ref_tier)`.
