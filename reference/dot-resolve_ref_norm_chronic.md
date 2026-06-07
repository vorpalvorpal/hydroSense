# Resolve ref_norm for chronic target samples (window-integrated)

For each chronic sample (`focal_date` column), integrates the predicted
reference at daily resolution over
`[focal_date - window_days, focal_date]` using the same
exponential-decay kernel as
[`time_weighted_aggregate()`](https://vorpalvorpal.github.io/leachatetools/reference/time_weighted_aggregate.md).
Tier 1 does not apply to chronic targets (no single reference grab is a
proxy for an entire integration window).

## Usage

``` r
.resolve_ref_norm_chronic(ref_model, target_df, tau_days, window_days)
```

## Arguments

- ref_model:

  A `reference_model` object.

- target_df:

  Tibble with columns `sample_id`, `focal_date`.

- tau_days:

  Exponential-decay parameter (days).

- window_days:

  Look-back window (days).

## Value

Tibble `(sample_id, analyte, ref_norm, ref_tier)`.
