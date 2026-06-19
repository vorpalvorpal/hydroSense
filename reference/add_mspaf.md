# Compute the Adjusted multi-substance PAF (msPAF) for water quality samples

Appends msPAF rows to a long-format water quality dataframe. msPAF
estimates the fraction of aquatic species potentially affected by the
combined toxicant mixture, adjusted for local geogenic background via
the Added Risk Approach. See the file-level header for full
methodological detail.

## Usage

``` r
add_mspaf(
  df,
  reference = NULL,
  analyte_metadata = NULL,
  method = c("multi", "anzecc"),
  guideline_dir = getOption("hydroSense.guideline_dir"),
  min_analytes = 3,
  ref_summary = c("geom_mean", "median", "arith_mean", "p80", "p90", "p95"),
  conc_units = NULL,
  require_temperature = TRUE,
  tau = 90,
  tau_units = "day",
  window = 365,
  window_units = "day",
  return = c("summary", "draws"),
  interval = 0.9,
  central = c("median", "mean")
)
```

## Arguments

- df:

  Long-format monitoring dataframe. Required columns: `sample_id`,
  `site_id`, `analyte`, `value`. Toxicant concentrations must ultimately
  be in µg/L for SSD lookup; supply them either via a `units.analyte`
  column (one unit string per row, e.g. `"mg/L"`) or via the
  `conc_units` argument (applied uniformly to all SSD-eligible rows).
  Optional but recommended: `datetime` (propagated to msPAF rows if
  present), `detected` (assumed `TRUE` if absent), `imputed` (logical;
  if present, `n_analytes_imputed` is populated in output). Driver
  analytes needed for chemistry normalisation (e.g. pH, DOC) should be
  present as rows in `df`.

- reference:

  Background reference chemistry for the ARA adjustment. Accepts four
  forms:

  - A `reference_model` from
    [`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md)
    — contemporaneous (temporal) ARA; the model predicts what the
    reference site would show at each target sample's exact moment using
    hydrology and seasonality.

  - A `prepared_reference` object from
    [`prepare_reference()`](https://vorpalvorpal.github.io/hydroSense/reference/prepare_reference.md)
    — normalisation has already been applied; used directly (static
    ARA).

  - A raw long-format data frame (same schema as `df`) — will be passed
    to
    [`prepare_reference()`](https://vorpalvorpal.github.io/hydroSense/reference/prepare_reference.md)
    internally (static ARA).

  - `NULL` (default) — no ARA adjustment; raw concentrations assessed
    directly against SSDs.

- analyte_metadata:

  Data frame of analyte metadata, or `NULL` to load the bundled
  `inst/extdata/anzecc_analyte_metadata.csv`. Passed to
  [`prepare_reference()`](https://vorpalvorpal.github.io/hydroSense/reference/prepare_reference.md)
  and
  [`derive_ssd_params()`](https://vorpalvorpal.github.io/hydroSense/reference/derive_ssd_params.md).

- method:

  SSD method passed to
  [`ssd_hc50()`](https://vorpalvorpal.github.io/hydroSense/reference/ssd_hc50.md)
  and
  [`ssd_paf()`](https://vorpalvorpal.github.io/hydroSense/reference/ssd_paf.md).
  `"multi"` (default) fits all 6 BCANZ distributions and model-averages;
  `"anzecc"` uses the per-analyte distribution matching the original
  ANZG derivation.

- guideline_dir:

  Path to the "guideline data" folder containing ANZG XLSX files. Falls
  back to `getOption("hydroSense.guideline_dir")`.

- min_analytes:

  Minimum number of analytes with fitted SSDs required to compute msPAF
  for a sample. Default `3`.

- ref_summary:

  Summary statistic for the reference distribution when `reference` is a
  raw data frame. Passed through to
  [`prepare_reference()`](https://vorpalvorpal.github.io/hydroSense/reference/prepare_reference.md).
  Default `"geom_mean"` — the maximum-likelihood central tendency for
  log-normal concentrations, and a PICT-consistent estimate of the
  "typical" exposure the resident community has adapted to. Other
  options: `"median"`, `"arith_mean"`, `"p80"`, `"p90"`, `"p95"`.

- conc_units:

  Character. Unit string (e.g. `"mg/L"`, `"ug/L"`) applied uniformly to
  all SSD-eligible rows in `df` when `df` has no `units.analyte` column.
  Ignored when `df` already carries `units.analyte`. Required when `df`
  lacks `units.analyte` and toxicant concentrations are not already in
  µg/L.

- require_temperature:

  Logical (default `TRUE`). When `TRUE`, any sample that reports an
  `NH3-N` measurement **must** also carry a water `temperature` row (the
  ammonia un-ionised-fraction normalisation is undefined without it); a
  missing temperature is a hard error rather than a silent drop of
  ammonia. Supply temperature via direct measurement, or derive it with
  [`estimate_water_temp()`](https://vorpalvorpal.github.io/hydroSense/reference/estimate_water_temp.md)
  (optionally fed by
  [`get_silo_air_temp()`](https://vorpalvorpal.github.io/hydroSense/reference/get_silo_air_temp.md)).
  Set `FALSE` only for datasets that do not assess ammonia.

- tau, tau_units:

  Exponential-decay half-life for chronic window integration when
  `reference` is a `reference_model`. Default `90` days.

- window, window_units:

  Look-back window length for chronic integration. Default `365` days.

- return:

  Output mode for draw-carrier input (see
  [`summarise_draws()`](https://vorpalvorpal.github.io/hydroSense/reference/summarise_draws.md)).
  `"summary"` (default) collapses posterior draws to a central estimate
  plus a credible interval (`value`, `value_lower`, `value_upper`,
  `n_draws`); `"draws"` returns the raw per-draw msPAF rows
  (`draw_id 1..N`) for external risk models or further composition (e.g.
  into
  [`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md)).
  For point (non-draw) input both modes return byte-identical output
  with no interval columns.

- interval:

  Width of the credible interval when `return = "summary"`. Default
  `0.90` (5th/95th percentile bounds).

- central:

  Central-tendency statistic when `return = "summary"`: `"median"`
  (default) or `"mean"`.

## Value

The input `df` with msPAF rows appended. Each msPAF row carries: `value`
(msPAF as a percentage, 0–100+), `detected = TRUE`, `analyte = "msPAF"`,
`n_analytes_used` (integer), `n_analytes_imputed` (integer, 0 if
`imputed` column absent), `dominant_analyte` (character), and `max_paf`
(numeric).

The result carries two attributes (read them before further dplyr
wrangling, which drops attributes): `"analyte_pafs"` — the per-analyte
PAF breakdown as a flat tibble (`site_id`, `sample_id`, `draw_id` in
draws mode, `analyte`, `C_adj`, `PAF`, `moa_group`, `ref_source`),
retrieved with
[`analyte_pafs()`](https://vorpalvorpal.github.io/hydroSense/reference/analyte_pafs.md);
and `"ara_summary"` — per-(sample × analyte) ARA diagnostics, retrieved
with
[`ara_summary()`](https://vorpalvorpal.github.io/hydroSense/reference/ara_summary.md).
(`analyte_pafs` was formerly a per-row list-column; it is now a flat
attribute — issue \#30.)

Tier breaks are not provided by this package — msPAF is a continuous
risk metric and the threshold at which a community is "impacted" depends
on the assessment context. See
[`vignette("chronic-mspaf-interpretation")`](https://vorpalvorpal.github.io/hydroSense/articles/chronic-mspaf-interpretation.md)
for guidance.

## Details

The function accepts either per-sample or chronic-integrated chemistry
(from
[`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md)).
It does not need to know which — the distinction is entirely in the
input data. Similarly, `reference` may be a raw long-format chemistry
data frame or a pre-built
[`prepare_reference()`](https://vorpalvorpal.github.io/hydroSense/reference/prepare_reference.md)
object.

## References

De Zwart D, Posthuma L (2005) Environmental Toxicology and Chemistry
24(10):2665-2676.

## See also

[`ssd_paf()`](https://vorpalvorpal.github.io/hydroSense/reference/ssd_paf.md),
[`ssd_hc50()`](https://vorpalvorpal.github.io/hydroSense/reference/ssd_hc50.md),
[`prepare_reference()`](https://vorpalvorpal.github.io/hydroSense/reference/prepare_reference.md),
[`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md),
[`ara_summary()`](https://vorpalvorpal.github.io/hydroSense/reference/ara_summary.md),
[`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md),
[`prescreen_analytes()`](https://vorpalvorpal.github.io/hydroSense/reference/prescreen_analytes.md),
[`impute_chemistry()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_chemistry.md)

## Examples

``` r
# \donttest{
# Per-sample multi-substance PAF for the impacted site, with local
# background subtracted via the reference site. Uses the bundled SSD data.
demo <- leachate_demo()
ds  <- subset(demo, site_id == "downstream")
ref <- subset(demo, site_id == "reference")
out <- add_mspaf(ds, reference = ref)
subset(out, analyte == "msPAF", c("sample_id", "value"))
#> # A tibble: 6 × 2
#>   sample_id value
#>   <chr>     <dbl>
#> 1 DS-01      92.7
#> 2 DS-02      95.1
#> 3 DS-03      96.0
#> 4 DS-04      94.2
#> 5 DS-05      95.8
#> 6 DS-06      95.9
# }
```
