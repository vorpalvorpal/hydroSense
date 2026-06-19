# Pre-screen analytes by detection frequency

Computes detection frequency (proportion of samples with
`detected == TRUE`) per analyte and returns the names that meet a
minimum threshold. Use this before
[`impute_chemistry()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_chemistry.md)
to drop analytes that were almost never detected — imputing such
analytes from near-zero priors adds noise without ecological signal.

## Usage

``` r
prescreen_analytes(
  df,
  k = 0.05,
  protect = NULL,
  potency_keep = TRUE,
  potency_frac = 1,
  group_by_feature = FALSE,
  conc_units = NULL,
  analyte_metadata = NULL,
  return = c("vector", "table")
)
```

## Arguments

- df:

  Long-format chemistry data frame. Required columns: `sample_id`
  (character), `analyte` (character), `detected` (logical). When
  `group_by_feature = TRUE`, also requires `site_id` (character).

- k:

  Minimum detection frequency (proportion, 0–1). Analytes with
  `n_detected / n_samples < k` are excluded unless protected. Default
  `0.05` (5 %).

- protect:

  Optional character vector of additional analyte names to protect from
  prescreen exclusion, on top of the metadata-derived co-analytes.
  Default `NULL`. Pass the `required_vars` from your downstream
  imputation step here.

- potency_keep:

  Logical. Enable the potency escape hatch (keep a frequency-failing
  analyte whose concentration reaches its guideline value). Default
  `TRUE`. Set `FALSE` for a frequency-only prescreen.

- potency_frac:

  Numeric (\>= 0). Fraction of the 95 % guideline value
  (`dgv_95pct_ug_L`) a detected concentration must reach to rescue an
  analyte. Default `1` (must reach the guideline). Lower it (e.g. `0.5`)
  to be more precautionary.

- group_by_feature:

  Logical. If `TRUE`, detection frequency is computed per `site_id` and
  an analyte is included only if it passes in *all* features (worst-case
  feature, pooled counts). If `FALSE` (default), frequency is pooled
  across all samples.

- conc_units:

  Character. Unit string (e.g. `"mg/L"`, `"ug/L"`) applied to toxicant
  concentrations in `df` when it has no `units.analyte` column. Only
  matters when `potency_keep = TRUE`, because the potency threshold
  comparison needs values in µg/L (matching `dgv_95pct_ug_L`). Ignored
  when `potency_keep = FALSE` or `df` already carries a `units.analyte`
  column.

- analyte_metadata:

  Data frame of analyte metadata, or `NULL` to load the bundled
  `inst/extdata/anzecc_analyte_metadata.csv`. Used to identify
  co-analytes that must be protected from exclusion.

- return:

  Either `"vector"` (default) to return a character vector of included
  analyte names, or `"table"` to return a tibble with one row per
  analyte showing detection statistics and inclusion flag.

## Value

When `return = "vector"`: a character vector of passing analyte names.
The vector carries an attribute `"excluded"` listing analytes that did
not pass the threshold (non-protected only), so callers can record what
was dropped.

When `return = "table"`: a tibble with columns `analyte`, `n_samples`,
`n_detected`, `detect_freq`, `limiting_site`, `protected` (logical),
`potency_kept` (logical; rescued by the potency escape hatch),
`included` (logical). When `group_by_feature = TRUE`, `n_samples` /
`n_detected` / `detect_freq` describe the *limiting* site (the
worst-case feature that drove the inclusion decision) and
`limiting_site` names it; when `group_by_feature = FALSE`, counts are
pooled across all samples and `limiting_site` is `NA`.

## Details

Analytes listed in `coanalytes_required` in the bundled metadata (e.g.
pH, DOC, Ca, Mg, hardness, temperature) are **automatically protected**
from exclusion regardless of detection frequency, because they are
needed for chemistry normalisation in
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md).
Additional analytes can be protected via the `protect` argument (typical
use: pass the `required_vars` you intend to use in
[`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md)
here so those vars survive prescreen). Protected analytes that fall
below the threshold are reported separately so the caller is aware.

**Potency escape hatch.** Frequency alone can screen out a
rare-but-potent toxicant (e.g. a pesticide detected in 2 % of samples
but at ecotoxicologically significant concentrations). With
`potency_keep = TRUE` (default), an analyte that fails the frequency
threshold is still kept if any detected concentration reaches
`potency_frac` times its 95 %-species- protection guideline value
(`dgv_95pct_ug_L` in the metadata). This needs a numeric `value` column
in the same units as the guideline (µg/L); only analytes that carry a
guideline value (toxicants) can be rescued this way, so major ions and
analytes with no guideline are unaffected.

## Examples

``` r
included <- prescreen_analytes(leachate_demo(), k = 0.05)
included                    # analytes retained
#>  [1] "CO3-CaCO3"   "Ca"          "Cl"          "Cu"          "DOC"        
#>  [6] "F"           "HCO3-CaCO3"  "K"           "Mg"          "NH3-N"      
#> [11] "NO2-N"       "NO3-N"       "Na"          "SO4²⁻"       "Zn"         
#> [16] "hardness"    "pH"          "temperature"
#> attr(,"excluded")
#> character(0)
attr(included, "excluded")  # see what was dropped (and why)
#> character(0)
```
