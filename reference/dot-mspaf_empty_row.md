# Construct an empty (zero-row) msPAF result tibble

Construct an empty (zero-row) msPAF result tibble

## Usage

``` r
.mspaf_empty_row(with_sample_id = FALSE)
```

## Arguments

- with_sample_id:

  Logical; if `TRUE`, prepend a `sample_id` column (used by the batched
  [`compute_mspaf_per_sample()`](https://vorpalvorpal.github.io/hydroSense/reference/compute_mspaf_per_sample.md)
  path).

## Value

A zero-row tibble with the msPAF output schema.
