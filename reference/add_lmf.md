# Compute the Leachate Mixing Fraction (LMF) for water quality samples

Appends LMF rows to a long-format water quality dataframe. LMF is a
source-detection index estimating what fraction of a sample's chemistry
can be attributed to the leachate end-member under a two-component
mixing model. See the file-level header for full methodological detail.

## Usage

``` r
add_lmf(
  df,
  leachate_data = NULL,
  reference_data = NULL,
  lmf_analyte_uuid = NA_character_,
  calibration_window_years = 5,
  min_leachate_total_n = NULL,
  min_leachate_total_n_units = NULL,
  rsd_default = 0.05,
  min_ref_samples = 10,
  min_leachate_samples = 10,
  max_sigma_lsi = 10,
  max_chi2_per_df = Inf,
  informativeness_threshold = 0.2,
  min_high_info_ions = 3L,
  robust_iterations = 3L,
  robust_threshold_k = 1.5,
  verbose = FALSE
)
```

## Arguments

- df:

  Long-format dataframe. Required columns: `sample_id`, `site_id`,
  `analyte`, `value`, `detected`, `datetime`.

- leachate_data:

  Optional long-format dataframe in the same structure as `df` providing
  leachate end-member chemistry. Required columns: `sample_id`,
  `analyte`, `value`, `detected`. When supplied, *all* samples are used
  to build the leachate end-member regardless of date and regardless of
  total-N content (the total-N quality filter is bypassed — supply
  curated data). Samples missing the required LMF panel analytes are
  silently dropped. `NULL` (default) uses the standard leachate feature
  detection logic (requires dashboard infrastructure: `feature_df()`,
  `data_df()`).

- reference_data:

  Optional long-format dataframe in the same structure as `df` providing
  reference site chemistry for end-member calibration. Required columns:
  `sample_id`, `analyte`, `value`, `detected`. When supplied, *all*
  samples are used to build the reference end-member regardless of date;
  the calibration window is ignored. The same end-member is applied
  globally to all features. `NULL` (default) uses the standard
  per-feature matched reference site logic (requires dashboard
  infrastructure: `feature_sfc()`, `data_df()`).

- lmf_analyte_uuid:

  UUID assigned to LMF rows in the `uuid.analyte` column of the output.
  Set to a fixed UUID if integrating with a data system that tracks
  analytes by UUID. Default `NA_character_`.

- calibration_window_years:

  Number of years of calibration data to use. Window is centred on the
  input data's date range and shifted backwards if future data are
  unavailable. Default `5`.

- min_leachate_total_n:

  Minimum total N (sum of NH4-N + NO3-N + NO2-N) for a leachate sample
  to qualify for end-member calibration. Excludes non-representative
  samples (dilution events, treatment upsets). Total N is used rather
  than NH4-N alone because aerated leachate features may have nitrified
  NH4 to NO3. Numeric or `units` object; bare numeric requires
  `min_leachate_total_n_units`. Default `NULL` → 20 mg/L.

- min_leachate_total_n_units:

  Character unit string for `min_leachate_total_n` when it is bare
  numeric, e.g. `"mg/L"`. Ignored when `min_leachate_total_n` is a
  `units` object or `NULL`.

- rsd_default:

  Default relative analytical uncertainty, applied as
  `sigma_meas = rsd_default * |x|` with a floor at `rsd_default * |R|`.
  Default `0.05` (5% RSD).

- min_ref_samples:

  Minimum reference samples required in the calibration window. Features
  matched to a reference site with fewer samples are skipped silently.
  Default `10`.

- min_leachate_samples:

  Minimum valid leachate samples required for end-member construction.
  The function stops with an informative error if not met. Default `10`.

- max_sigma_lsi:

  Maximum permitted `sigma_lmf` in percentage points (uncertainty on the
  mixing fraction estimate). Samples above this threshold return `NA`
  with reason code `"insufficient_precision"`. Default `10` (i.e., ±10
  percentage points).

- max_chi2_per_df:

  Maximum permitted chi-squared per degree of freedom (computed on the
  *original* inverse-variance weights, before robust reweighting).
  Retained as a diagnostic output in all cases. Default `Inf` (no hard
  gate — the robust reweighting handles outlier ions directly, making a
  hard chi2 gate redundant for most samples).

- informativeness_threshold:

  Threshold on `sigma_R / |L - R|` for classifying ions as
  high-information. Computed once from calibration data, not per-sample.
  Lower values are more stringent. Default `0.20`.

- min_high_info_ions:

  Minimum number of high-information ions that must be measured in a
  sample for LMF to be computed. Samples below this threshold return
  `NA` with reason code `"insufficient_high_info_ions"`. Default `3L`.

- robust_iterations:

  Number of Huber M-estimator reweighting passes. In each pass, ions
  whose residuals exceed `robust_threshold_k * MAD(residuals)` are
  downweighted proportionally. Three iterations is sufficient for
  convergence in practice. Set to `0L` to disable robust reweighting
  entirely. Default `3L`.

- robust_threshold_k:

  Threshold multiplier applied to the unweighted median absolute
  deviation (MAD) of per-ion residuals. Ions with `|r_i| > k * MAD` have
  their estimation weight reduced by a factor of `k * MAD / |r_i|`.
  Smaller values downweight more aggressively; `k = 1.5` is fairly
  tight, `k = 2.5` is close to the conventional Huber constant. Default
  `1.5`.

- verbose:

  If `TRUE`, prints an ion informativeness table via
  [`cli::cli_inform()`](https://cli.r-lib.org/reference/cli_abort.html)
  showing each ion's score and high/low classification. Useful when
  tuning `informativeness_threshold`. Default `FALSE`.

## Value

The input `df` with LMF rows appended. Each LMF row carries `value` (the
robust LMF estimate as a percentage, 0 = pure reference, 100 = pure
leachate), `name.analyte = "LMF"`, `uuid.feature` (from the input),
`lmf_naive` (the non-robust estimate for comparison), `lmf_reason` (`NA`
on success, reason code on failure), `n_ions_used`,
`n_ions_downweighted` (count of ions whose weight was reduced by robust
reweighting), `sigma_lmf`, and `chi2_per_df` (diagnostic; computed on
original weights). Columns present in `df` but not produced by the LMF
computation are `NA` in the appended rows.

## References

Christophersen N, Hooper RP (1992) Water Resources Research
28(1):99-107. Christensen JB et al. (2001) Applied Geochemistry
16(7-8):659-718. Kjeldsen P et al. (2002) Crit Rev Environ Sci Technol
32(4).

## See also

[`build_leachate_endmember`](https://vorpalvorpal.github.io/hydroSense/reference/build_leachate_endmember.md),
[`build_reference_endmember`](https://vorpalvorpal.github.io/hydroSense/reference/build_reference_endmember.md),
[`compute_lmf_for_sample`](https://vorpalvorpal.github.io/hydroSense/reference/compute_lmf_for_sample.md)

## Examples

``` r
# \donttest{
# Quantify each downstream sample's leachate-mixing fraction from its
# major-ion signature, calibrated against the leachate and reference
# end-members in the bundled demo data.
demo <- leachate_demo()
out <- add_lmf(
  df             = subset(demo, site_id == "downstream"),
  leachate_data  = subset(demo, site_id == "leachate"),
  reference_data = subset(demo, site_id == "reference")
)
subset(out, analyte == "LMF", c("sample_id", "value"))
#> # A tibble: 6 × 2
#>   sample_id value
#>   <chr>     <dbl>
#> 1 DS-01      15.1
#> 2 DS-02      14.9
#> 3 DS-03      14.2
#> 4 DS-04      14.9
#> 5 DS-05      14.7
#> 6 DS-06      15.0
# }
```
