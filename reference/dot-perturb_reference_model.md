# Perturb all GAM fits inside a reference_model

Draws independent posterior-coefficient samples for each analyte's
`gamm_fit`. Used when `mspaf_daily(reference = NULL)` (total-
concentration mode), where `ref_norm` enters the result directly and its
uncertainty should be reflected. When background is subtracted via the
same reference model, the reference-GAM contribution cancels and this
function should NOT be called.

## Usage

``` r
.perturb_reference_model(rm)
```

## Arguments

- rm:

  A `reference_model` from
  [`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md).

## Value

A copy of `rm` with each `gamm_fit` replaced by a posterior draw.
