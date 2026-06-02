# Construct an empty (zero-row) AmsPAF result tibble

Construct an empty (zero-row) AmsPAF result tibble

## Usage

``` r
.amspaf_empty_row(with_sample_id = FALSE)
```

## Arguments

- with_sample_id:

  Logical; if `TRUE`, prepend a `sample_id` column (used by the batched
  [`compute_amspaf_per_sample()`](https://www.kedumba.com.au/leachatetools/reference/compute_amspaf_per_sample.md)
  path).

## Value

A zero-row tibble with the AmsPAF output schema.
