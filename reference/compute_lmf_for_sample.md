# Compute LMF for a single sample

Internal function called once per sample by
[`add_lmf`](https://www.kedumba.com.au/leachatetools/reference/add_lmf.md).
Computes the inverse-variance weighted mixing fraction and associated
quality metrics, applying admission and quality gates before returning.

## Usage

``` r
compute_lmf_for_sample(
  sample_wide,
  endmembers,
  high_info_ions,
  min_high_info_ions,
  rsd_default,
  max_sigma,
  max_chi2_df,
  uuid_lmf,
  ref_window_start,
  ref_window_end,
  leach_window_start,
  leach_window_end,
  cl_anchor,
  robust_iterations = 3L,
  robust_threshold_k = 1.5,
  verbose = FALSE,
  datetime_sample = NA_character_
)
```

## Arguments

- sample_wide:

  One-row wide-format tibble with one column per measured (collapsed)
  panel ion.

- endmembers:

  Tibble with columns `ion`, `R`, `sigma_R`, `n_ref`, `L`. One row per
  ion with usable end-member data for this feature.

- high_info_ions:

  Character vector of ion names classified as high-information at
  calibration time.

- min_high_info_ions:

  Minimum count of `high_info_ions` that must be present and non-NA in
  this sample.

- rsd_default:

  Default relative analytical uncertainty (fraction).

- max_sigma:

  Maximum permitted `sigma_lmf`.

- max_chi2_df:

  Maximum permitted chi-squared per degree of freedom.

- uuid_lmf:

  UUID of the LMF analyte entry in `analyteDF`.

- ref_window_start, ref_window_end:

  Calibration window dates for the reference end-member; stored in
  guideline comment strings.

- leach_window_start, leach_window_end:

  Calibration window dates for the leachate end-member; stored in
  guideline comment strings.

- cl_anchor:

  Cl concentration (meq/L) used to anchor the leachate end-member;
  stored in guideline comment strings.

- robust_iterations:

  Number of Huber reweighting passes. Inherited from
  [`add_lmf`](https://www.kedumba.com.au/leachatetools/reference/add_lmf.md).

- robust_threshold_k:

  MAD multiplier for Huber downweighting. Inherited from
  [`add_lmf`](https://www.kedumba.com.au/leachatetools/reference/add_lmf.md).

- verbose:

  If `TRUE`, emits a per-ion diagnostic table via
  [`cli::cli_inform()`](https://cli.r-lib.org/reference/cli_abort.html)
  showing observed concentration, end-member values, per-ion mixing
  fraction, original and robust weights. Inherited from
  [`add_lmf`](https://www.kedumba.com.au/leachatetools/reference/add_lmf.md).
  Default `FALSE`.

- datetime_sample:

  The sample datetime as a character string, used to label the per-ion
  diagnostic table when `verbose = TRUE`. Default `NA_character_`.

## Value

A one-row tibble. On success: `value` is the robust LMF estimate,
`quantified = TRUE`. On failure: `value = NA`, `quantified = FALSE`,
`lmf_reason` carries a descriptive reason code. All paths return
identical column structure.
