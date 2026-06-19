# Prepare reference chemistry for msPAF background subtraction

Applies chemistry normalisation (if formulas are populated in the
metadata) and computes a per-analyte central-tendency summary from
reference-site data, optionally with a bootstrap confidence interval.
The resulting object is passed as the `reference` argument to
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md).

## Usage

``` r
prepare_reference(
  reference_data,
  analyte_metadata = NULL,
  summary = c("geom_mean", "median", "arith_mean", "p80", "p90", "p95"),
  bootstrap_ci = FALSE,
  n_boot = 1000L,
  conc_units = NULL,
  eps = 1e-09
)
```

## Arguments

- reference_data:

  Long-format chemistry data frame for the reference (background)
  site(s). Same schema as the input to
  [`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md):
  `analyte`, `value`, `detected`. Toxicant concentrations must be in
  µg/L before normalisation; supply them either via a `units.analyte`
  column or via the `conc_units` argument. BDL (`detected == FALSE`)
  observations contribute `0` to the summary statistic.

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

- conc_units:

  Character. Unit string (e.g. `"mg/L"`, `"ug/L"`) applied uniformly to
  all SSD-eligible rows in `reference_data` when it has no
  `units.analyte` column. Ignored when `reference_data` carries
  `units.analyte`. Required when the data lacks `units.analyte` and
  toxicant concentrations are not already in µg/L.

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
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)
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
ref <- subset(leachate_demo(), site_id == "reference")

# Default: geometric mean (recommended)
prep_ref <- prepare_reference(ref)
prep_ref$ref_table
#> # A tibble: 18 × 3
#>    analyte     ref_norm n_obs
#>    <chr>          <dbl> <int>
#>  1 CO3-CaCO3     4.11       6
#>  2 Ca           21.8        6
#>  3 Cl           18.0        6
#>  4 Cu            0.257      6
#>  5 DOC           2.02       6
#>  6 F             0.201      6
#>  7 HCO3-CaCO3   54.3        6
#>  8 K             3.01       6
#>  9 Mg            6.08       6
#> 10 NH3-N       150.         6
#> 11 NO2-N         0.0302     6
#> 12 NO3-N         0.790      6
#> 13 Na           12.1        6
#> 14 SO4²⁻        12.0        6
#> 15 Zn            1.75       6
#> 16 hardness     80.3        6
#> 17 pH            7.62       6
#> 18 temperature  15.6        6

# Or a higher percentile of the local background distribution:
prepare_reference(ref, summary = "p80")$ref_table
#> # A tibble: 18 × 3
#>    analyte     ref_norm n_obs
#>    <chr>          <dbl> <int>
#>  1 CO3-CaCO3     4.13       6
#>  2 Ca           22.1        6
#>  3 Cl           18.5        6
#>  4 Cu            0.262      6
#>  5 DOC           2.07       6
#>  6 F             0.206      6
#>  7 HCO3-CaCO3   54.4        6
#>  8 K             3.04       6
#>  9 Mg            6.28       6
#> 10 NH3-N       159.         6
#> 11 NO2-N         0.0308     6
#> 12 NO3-N         0.804      6
#> 13 Na           12.3        6
#> 14 SO4²⁻        12.3        6
#> 15 Zn            1.80       6
#> 16 hardness     81.3        6
#> 17 pH            7.64       6
#> 18 temperature  16.3        6
```
