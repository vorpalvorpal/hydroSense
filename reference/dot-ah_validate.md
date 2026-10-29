# Validate the added-hazard arguments of mspaf_pipeline()

Raises every argument error of the `ara_form = "added_hazard"` path
before anything is computed, and resolves `ara` and `impute`.

## Usage

``` r
.ah_validate(
  reference,
  impute,
  impute_missing,
  daily_args,
  chronic_summary,
  ara,
  ara_missing,
  min_coverage,
  tau,
  window
)
```

## Value

A list `(ara, reference)`; `reference` is `NULL` when `ara = "none"`.
