# Build a system-wide pooled reference end-member

Pools data from all reference features to compute per-ion mean and
standard deviation. Used exclusively for computing ion informativeness
scores in
[`add_lmf`](https://vorpalvorpal.github.io/hydroSense/reference/add_lmf.md);
per-feature reference end-members are built by
[`build_reference_endmember`](https://vorpalvorpal.github.io/hydroSense/reference/build_reference_endmember.md).

## Usage

``` r
build_pooled_reference_endmember(calibration_start, calibration_end)
```

## Arguments

- calibration_start:

  Calibration window start date.

- calibration_end:

  Calibration window end date. Shifted backwards if reference data do
  not extend this far.

## Value

Tibble with columns `ion`, `n_ref`, `R`, `sigma_R`. Only ions with \>= 3
observations across all reference features are included.
