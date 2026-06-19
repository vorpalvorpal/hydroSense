# Resolve reference norms for a block of target samples

Dispatches to
[`.resolve_ref_norm_instant()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-resolve_ref_norm_instant.md)
or
[`.resolve_ref_norm_chronic()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-resolve_ref_norm_chronic.md)
depending on whether `df` contains a `focal_date` column. Returns a
tibble `(sample_id, analyte, ref_norm, ref_tier)` that
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)
uses instead of the static `(analyte, ref_norm)` table produced by
[`prepare_reference()`](https://vorpalvorpal.github.io/hydroSense/reference/prepare_reference.md).

## Usage

``` r
.resolve_ref_norm(ref_model, df, tau_days = 90, window_days = 365)
```

## Arguments

- ref_model:

  A `reference_model` object.

- df:

  Target chemistry data frame (as passed to
  [`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)).

- tau_days:

  Exponential-decay parameter for chronic integration.

- window_days:

  Window length for chronic integration.

## Value

Tibble `(sample_id, analyte, ref_norm, ref_tier)`.
