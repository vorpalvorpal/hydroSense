# Compute the Adjusted multi-substance PAF (AmsPAF) for water quality samples

Appends AmsPAF rows to a long-format water quality dataframe. AmsPAF
estimates the fraction of aquatic species potentially affected by the
combined toxicant mixture, adjusted for local geogenic background via
the Added Risk Approach. See the file-level header for full
methodological detail.

## Usage

``` r
add_amspaf(
  df,
  reference = NULL,
  analyte_metadata = NULL,
  method = c("multi", "anzecc"),
  guideline_dir = getOption("leachatetools.guideline_dir"),
  min_analytes = 3,
  ref_summary = c("geom_mean", "median", "arith_mean", "p80", "p90", "p95"),
  require_temperature = TRUE
)
```

## Arguments

- df:

  Long-format monitoring dataframe. Required columns: `sample_id`,
  `site_id`, `analyte`, `value` (concentrations in µg/L). Optional but
  recommended: `datetime` (propagated to AmsPAF rows if present),
  `detected` (assumed `TRUE` if absent), `imputed` (logical; if present,
  `n_analytes_imputed` is populated in output). Driver analytes needed
  for chemistry normalisation (e.g. pH, DOC) should be present as rows
  in `df`.

- reference:

  Background reference chemistry for the ARA adjustment. Accepts three
  forms:

  - A `prepared_reference` object from
    [`prepare_reference()`](https://www.kedumba.com.au/leachatetools/reference/prepare_reference.md)
    — normalisation has already been applied; used directly.

  - A raw long-format data frame (same schema as `df`) — will be passed
    to
    [`prepare_reference()`](https://www.kedumba.com.au/leachatetools/reference/prepare_reference.md)
    internally.

  - `NULL` (default) — no ARA adjustment; raw concentrations assessed
    directly against SSDs.

- analyte_metadata:

  Data frame of analyte metadata, or `NULL` to load the bundled
  `inst/extdata/anzecc_analyte_metadata.csv`. Passed to
  [`prepare_reference()`](https://www.kedumba.com.au/leachatetools/reference/prepare_reference.md)
  and
  [`derive_ssd_params()`](https://www.kedumba.com.au/leachatetools/reference/derive_ssd_params.md).

- method:

  SSD method passed to
  [`ssd_hc50()`](https://www.kedumba.com.au/leachatetools/reference/ssd_hc50.md)
  and
  [`ssd_paf()`](https://www.kedumba.com.au/leachatetools/reference/ssd_paf.md).
  `"multi"` (default) fits all 6 BCANZ distributions and model-averages;
  `"anzecc"` uses the per-analyte distribution matching the original
  ANZG derivation.

- guideline_dir:

  Path to the "guideline data" folder containing ANZG XLSX files. Falls
  back to `getOption("leachatetools.guideline_dir")`.

- min_analytes:

  Minimum number of analytes with fitted SSDs required to compute AmsPAF
  for a sample. Default `3`.

- ref_summary:

  Summary statistic for the reference distribution when `reference` is a
  raw data frame. Passed through to
  [`prepare_reference()`](https://www.kedumba.com.au/leachatetools/reference/prepare_reference.md).
  Default `"geom_mean"` — the maximum-likelihood central tendency for
  log-normal concentrations, and a PICT-consistent estimate of the
  "typical" exposure the resident community has adapted to. Other
  options: `"median"`, `"arith_mean"`, `"p80"`, `"p90"`, `"p95"`.

- require_temperature:

  Logical (default `TRUE`). When `TRUE`, any sample that reports an
  `NH3-N` measurement **must** also carry a water `temperature` row (the
  ammonia un-ionised-fraction normalisation is undefined without it); a
  missing temperature is a hard error rather than a silent drop of
  ammonia. Supply temperature via direct measurement, or derive it with
  [`estimate_water_temp()`](https://www.kedumba.com.au/leachatetools/reference/estimate_water_temp.md)
  (optionally fed by
  [`get_silo_air_temp()`](https://www.kedumba.com.au/leachatetools/reference/get_silo_air_temp.md)).
  Set `FALSE` only for datasets that do not assess ammonia.

## Value

The input `df` with AmsPAF rows appended. Each AmsPAF row carries:
`value` (AmsPAF as a percentage, 0–100+), `detected = TRUE`,
`analyte = "AmsPAF"`, `n_analytes_used` (integer), `n_analytes_imputed`
(integer, 0 if `imputed` column absent), `dominant_analyte` (character),
`max_paf` (numeric), `analyte_pafs` (list column of per-analyte
diagnostic tibbles, each with `analyte`, `C_adj`, `PAF`, `moa_group`,
and `ref_source` — one of `"disabled"`, `"matched"`, `"unmatched"`
recording how the ARA reference was resolved for that analyte).

Tier breaks are not provided by this package — AmsPAF is a continuous
risk metric and the threshold at which a community is "impacted" depends
on the assessment context. See
[`vignette("chronic-amspaf-interpretation")`](https://www.kedumba.com.au/leachatetools/articles/chronic-amspaf-interpretation.md)
for guidance.

## Details

The function accepts either per-sample or chronic-integrated chemistry
(from
[`time_weighted_aggregate()`](https://www.kedumba.com.au/leachatetools/reference/time_weighted_aggregate.md)).
It does not need to know which — the distinction is entirely in the
input data. Similarly, `reference` may be a raw long-format chemistry
data frame or a pre-built
[`prepare_reference()`](https://www.kedumba.com.au/leachatetools/reference/prepare_reference.md)
object.

## References

De Zwart D, Posthuma L (2005) Environmental Toxicology and Chemistry
24(10):2665-2676.

## See also

[`ssd_paf()`](https://www.kedumba.com.au/leachatetools/reference/ssd_paf.md),
[`ssd_hc50()`](https://www.kedumba.com.au/leachatetools/reference/ssd_hc50.md),
[`prepare_reference()`](https://www.kedumba.com.au/leachatetools/reference/prepare_reference.md),
[`time_weighted_aggregate()`](https://www.kedumba.com.au/leachatetools/reference/time_weighted_aggregate.md),
[`prescreen_analytes()`](https://www.kedumba.com.au/leachatetools/reference/prescreen_analytes.md),
[`impute_chemistry()`](https://www.kedumba.com.au/leachatetools/reference/impute_chemistry.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Long-format monitoring data: one row per sample x analyte.
obs <- tibble::tibble(
  sample_id = c("S1", "S1", "S1"),
  site_id   = "downstream",
  analyte   = c("Cu", "Zn", "temperature"),
  value     = c(3.2, 18, 19)
)
ref <- tibble::tibble(
  sample_id = "R1", site_id = "upstream",
  analyte   = c("Cu", "Zn"), value = c(1.1, 6)
)
options(leachatetools.guideline_dir = "path/to/guideline data")
add_amspaf(obs, reference = ref)
} # }
```
