# Build a per-feature reference end-member for LMF calibration

Computes mean (`R_i`) and standard deviation (`sigma_R_i`) for each ion
in the collapsed LMF panel from reference samples within the calibration
window. Called once per feature by
[`add_lmf`](https://vorpalvorpal.github.io/hydroSense/reference/add_lmf.md).

## Usage

``` r
build_reference_endmember(reference_feature_uuid, cal_start, cal_end)
```

## Arguments

- reference_feature_uuid:

  UUID of the matched reference feature.

- cal_start:

  Ideal calibration window start date.

- cal_end:

  Ideal calibration window end date. Shifted backwards if reference data
  do not extend this far.

## Value

`NULL` if no reference data are found. Otherwise a named list with
elements:

- stats:

  Tibble with columns `ion`, `n_ref`, `R`, `sigma_R`. Only ions with \>=
  3 observations are included.

- window_start:

  Actual calibration window start used.

- window_end:

  Actual calibration window end used.

- n_samples:

  Number of sample events contributing.
