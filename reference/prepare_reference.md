# Prepare reference chemistry for AmsPAF background subtraction

Applies chemistry normalisation (if formulas are populated in the
metadata) and computes a per-analyte central-tendency summary from
reference-site data, optionally with a bootstrap confidence interval.
The resulting object is passed as the `reference` argument to
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md).

## Usage

``` r
prepare_reference(
  reference_data,
  analyte_metadata = NULL,
  summary = c("geom_mean", "median", "arith_mean", "p80", "p90", "p95"),
  bootstrap_ci = FALSE,
  n_boot = 1000L,
  eps = 1e-09
)
```

## Arguments

- reference_data:

  Long-format chemistry data frame for the reference (background)
  site(s). Same schema as the input to
  [`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md):
  `analyte`, `value`, `detected`. BDL (`detected == FALSE`) observations
  contribute `0` to the summary statistic.

- analyte_metadata:

  Data frame of analyte metadata, or `NULL` to load the bundled
  `inst/extdata/anzecc_analyte_metadata.csv`. Accepts either a data
  frame or a file path string. Must contain columns `analyte`,
  `coanalytes_required`, and `normalisation_formula`.

- summary:

  Summary statistic for the reference distribution. One of `"geom_mean"`
  (default), `"median"`, `"arith_mean"`, `"p80"`, `"p90"`, `"p95"`.

- bootstrap_ci:

  Logical. If `TRUE`, compute a 95% bootstrap CI on the reference
  summary statistic for each analyte (1,000 replicates by default). Adds
  `ref_lower`, `ref_upper`, and `n_boot_valid` columns to `$ref_table`.
  Default `FALSE`.

- n_boot:

  Number of bootstrap replicates if `bootstrap_ci = TRUE`. Default
  `1000L`.

- eps:

  Small positive guard added inside the log for geometric-mean
  computation, to handle BDL contributions of `0`. Default `1e-9`.

## Value

A list of class `"prepared_reference"` with elements:

- `$ref_table`: tibble with columns `analyte`, `ref_norm` (normalised
  summary concentration), `n_obs` (count of observations contributing).
  If `bootstrap_ci = TRUE`, also `ref_lower` / `ref_upper` (95% CI) and
  `n_boot_valid` (number of bootstrap draws that yielded a finite
  summary; CIs built from far fewer than `n_boot` draws are flagged with
  a warning).

- `$dropped`: character vector of analytes excluded due to no reference
  observations.

- `$summary`: the summary statistic used.

## Details

This is a pure function — it has no side effects and no internal cache.
In the chronic pipeline, call it once after computing chronic chemistry
for the reference feature(s), then pass the object into
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md)
for every focal date. In the per-sample pipeline, call it once on the
raw reference chemistry.

**Summary statistic.** The default is `"geom_mean"` — the geometric mean
of all detected observations. This is preferred over a fixed quantile
because:

- it is the maximum-likelihood central tendency for log-normal
  concentrations (which is how aquatic concentration data typically
  distribute);

- it uses all observations rather than a single ranked point, so it is
  more robust to small reference datasets;

- it is PICT-consistent: the resident community has adapted to the
  integrated typical exposure over time, not to a particular upper
  quantile.

BDL observations contribute `0` to the geometric mean via an
\\\epsilon\\-shifted log: `exp(mean(log(value + eps)))`. Other summaries
available: `"median"`, `"arith_mean"`, `"p80"`, `"p90"`, `"p95"`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Default: geometric mean (recommended)
prep_ref <- prepare_reference(ref_df)

# 80th percentile
prep_ref <- prepare_reference(ref_df, summary = "p80")

# With bootstrap CI
prep_ref <- prepare_reference(ref_df, bootstrap_ci = TRUE)
prep_ref$ref_table
} # }
```
