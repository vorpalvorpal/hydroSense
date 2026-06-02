# Build the Cl-anchored leachate end-member for LMF calibration

Constructs the leachate end-member L by computing the mean ratio of each
ion to Cl across valid leachate samples, then anchoring to the median Cl
concentration. This makes the end-member robust to EC variability
between leachate samples while preserving the stable compositional
fingerprint.

## Usage

``` r
build_leachate_endmember(
  calibration_start,
  calibration_end,
  min_leachate_total_n_mgl,
  min_leachate_samples
)
```

## Arguments

- calibration_start:

  Calibration window start date.

- calibration_end:

  Calibration window end date. Shifted backwards if leachate data do not
  extend this far.

- min_leachate_total_n_mgl:

  Minimum total N (mg/L) for a leachate sample to be considered valid.
  Applied in original mg/L units before meq conversion. Filters out
  anomalous samples not representing genuine leachate.

- min_leachate_samples:

  Minimum number of valid leachate samples required. Stops with an error
  if not met.

## Value

A named list with elements:

- L_values:

  Tibble with columns `ion`, `mean_ratio`, `L` (the end-member
  concentration in meq/L for each ion).

- cl_anchor:

  Median Cl concentration (meq/L) used to anchor the end-member.

- n_samples:

  Number of valid leachate samples used.

- f_included:

  Logical; `TRUE` if F was included in the end-member.

- window_start:

  Actual calibration window start used.

- window_end:

  Actual calibration window end used.
