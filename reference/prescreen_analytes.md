# Pre-screen analytes by detection frequency

Computes detection frequency (proportion of samples with
`detected == TRUE`) per analyte and returns the names that meet a
minimum threshold. Use this before
[`impute_chemistry()`](https://www.kedumba.com.au/leachatetools/reference/impute_chemistry.md)
to drop analytes that were almost never detected — imputing such
analytes from near-zero priors adds noise without ecological signal.

## Usage

``` r
prescreen_analytes(
  df,
  k = 0.05,
  protect = NULL,
  group_by_feature = FALSE,
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

- group_by_feature:

  Logical. If `TRUE`, detection frequency is computed per `site_id` and
  an analyte is included only if it passes in *all* features (worst-case
  feature, pooled counts). If `FALSE` (default), frequency is pooled
  across all samples.

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
[`add_amspaf()`](https://www.kedumba.com.au/leachatetools/reference/add_amspaf.md).
Additional analytes can be protected via the `protect` argument (typical
use: pass the `required_vars` you intend to use in
[`fit_imputation_model()`](https://www.kedumba.com.au/leachatetools/reference/fit_imputation_model.md)
here so those vars survive prescreen). Protected analytes that fall
below the threshold are reported separately so the caller is aware.

## Examples

``` r
if (FALSE) { # \dontrun{
included <- prescreen_analytes(chemistry, k = 0.05)
attr(included, "excluded")  # see what was dropped
chem_f <- dplyr::filter(chemistry, analyte %in% included)
} # }
```
