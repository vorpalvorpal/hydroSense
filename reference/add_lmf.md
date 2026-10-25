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
  calibration = c("window", "trailing"),
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
  robust_threshold_k = 2.5,
  verbose = FALSE
)
```

## Arguments

- df:

  Long-format dataframe. Required columns: `sample_id`, `site_id`,
  `analyte`, `value`, `detected`, `datetime`, and a unit column
  (`units.analyte` or sampleTidy's `units`). Analyte names may be
  hydroSense's or sampleTidy's (see
  [`analyte_aliases()`](https://vorpalvorpal.github.io/hydroSense/reference/analyte_aliases.md));
  `valence.analyte` and `atomic_mass.analyte` are optional (see
  [`to_meq()`](https://vorpalvorpal.github.io/hydroSense/reference/to_meq.md)).
  The LMF rows are appended to `df` as supplied.

- leachate_data:

  Required. Long-format dataframe in the same structure as `df`
  providing leachate end-member chemistry. Required columns:
  `sample_id`, `analyte`, `value`, `detected` (and `datetime` for
  `calibration = "trailing"`). The total-N quality filter
  (`min_leachate_total_n`) applies in both calibration modes: leachate
  samples below it are ignored. Samples missing the required LMF panel
  analytes are silently dropped.

- reference_data:

  Required. Long-format dataframe in the same structure as `df`
  providing reference site chemistry for end-member calibration.
  Required columns: `sample_id`, `analyte`, `value`, `detected` (and
  `datetime` for `calibration = "trailing"`). The same end-member is
  applied to every site in `df`.

- lmf_analyte_uuid:

  UUID assigned to LMF rows in the `uuid.analyte` column of the output.
  Set to a fixed UUID if integrating with a data system that tracks
  analytes by UUID. Default `NA_character_`.

- calibration:

  How the end-members are calibrated. `"window"` (default) builds one
  reference and one leachate end-member from *all* of `reference_data` /
  `leachate_data`, regardless of date, and applies them to every sample.
  `"trailing"` calibrates each sample separately, from only the rows of
  `leachate_data` and `reference_data` whose `datetime` lies in
  `[t - calibration_window_years, t]` (inclusive at both ends, where `t`
  is the sample's own `datetime`). Bounds are compared as whole calendar
  days, each `datetime`'s day taken in its own time zone (UTC when it
  has none), so clock time never matters. A sample's LMF therefore
  depends on data up to its own date only, and appending later samples
  or later end-member data never changes it. A sample whose window has
  no usable reference or leachate data (none, too few reference samples
  to estimate spread, too few samples for `min_ref_samples` or
  `min_leachate_samples`, or too few high-information ions) gets an `NA`
  LMF with `lmf_reason` starting `"insufficient_calibration_data"`.
  `"trailing"` requires a `datetime` column in both end-member data
  frames.

- calibration_window_years:

  Length in years of the trailing calibration window. Used only when
  `calibration = "trailing"`. Default `5`.

- min_leachate_total_n:

  Minimum total N (sum of NH3-N + NO3-N + NO2-N, each converted to mg/L
  by its row's unit) for a leachate sample to qualify for end-member
  calibration. Excludes non-representative samples (dilution events,
  treatment upsets); samples below it do not contribute to the leachate
  end-member at all. Total N is used rather than NH3-N alone because
  aerated leachate may have nitrified NH3 to NO3. Numeric or `units`
  object; bare numeric requires `min_leachate_total_n_units`. Default
  `NULL` (20 mg/L).

- min_leachate_total_n_units:

  Character unit string for `min_leachate_total_n` when it is bare
  numeric, e.g. `"mg/L"`. Ignored when `min_leachate_total_n` is a
  `units` object or `NULL`.

- rsd_default:

  Default relative analytical uncertainty, applied as
  `sigma_meas = rsd_default * |x|` with a floor at `rsd_default * |R|`.
  Default `0.05` (5% RSD).

- min_ref_samples:

  Minimum number of reference samples (distinct `sample_id`s carrying at
  least one LMF panel ion) in the calibration data – in the sample's
  trailing window when `calibration = "trailing"`. Below it the LMF is
  `NA` with `lmf_reason` `"insufficient_calibration_data"`. Default
  `10`.

- min_leachate_samples:

  Minimum number of valid leachate samples (total N at least
  `min_leachate_total_n`) in the calibration data (in the trailing
  window for `calibration = "trailing"`); below it the LMF is `NA` with
  `lmf_reason` `"insufficient_calibration_data"`. Also the number of
  valid samples with fluoride measured needed for F to join the panel.
  Default `10`.

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
  whose studentized residuals exceed `robust_threshold_k` are
  downweighted proportionally. Three iterations is sufficient for
  convergence in practice. Set to `0L` to disable robust reweighting
  entirely. Default `3L`.

- robust_threshold_k:

  Huber threshold in standard-deviation units on studentized residuals
  (`r_i / sqrt(var_f_i)`). Ions with `|sr_i| > k` have their estimation
  weight reduced by a factor of `k / |sr_i|`. Smaller values downweight
  more aggressively; `k = 2.5` is close to the conventional Huber
  constant. Default `2.5`.

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
reweighting), `n_high_info_downweighted` (count of downweighted ions
that were high-information tracers; a robust estimate carried by the
low-information majority against sensitive tracers signals a
mis-attributed source, not a transformation artifact), `sigma_lmf`, and
`chi2_per_df` (diagnostic; computed on original weights). Columns
present in `df` but not produced by the LMF computation are `NA` in the
appended rows.

## References

Christophersen N, Hooper RP (1992) Water Resources Research
28(1):99-107. Christensen JB et al. (2001) Applied Geochemistry
16(7-8):659-718. Kjeldsen P et al. (2002) Crit Rev Environ Sci Technol
32(4).

## See also

[`compute_lmf_for_sample`](https://vorpalvorpal.github.io/hydroSense/reference/compute_lmf_for_sample.md),
[`to_meq`](https://vorpalvorpal.github.io/hydroSense/reference/to_meq.md),
[`analyte_aliases`](https://vorpalvorpal.github.io/hydroSense/reference/analyte_aliases.md)

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
  reference_data = subset(demo, site_id == "reference"),
  # the demo has 6 samples per end-member; the defaults ask for 10
  min_ref_samples      = 5,
  min_leachate_samples = 5
)
subset(out, analyte == "LMF", c("sample_id", "value"))
#> # A tibble: 6 × 2
#>   sample_id value
#>   <chr>     <dbl>
#> 1 DS-01      15.0
#> 2 DS-02      14.8
#> 3 DS-03      14.7
#> 4 DS-04      14.9
#> 5 DS-05      14.7
#> 6 DS-06      14.9
# }
```
