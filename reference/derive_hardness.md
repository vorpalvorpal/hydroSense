# Three-way reconciliation between Ca, Mg, and hardness

Water hardness (total, expressed as CaCO3-equivalent mg/L) is
stoichiometric in Ca and Mg:

## Usage

``` r
derive_hardness(df, tolerance = 0.05, verbose = TRUE)
```

## Arguments

- df:

  Long-format chemistry data frame with columns `sample_id`, `analyte`,
  `value`, `detected`. Required columns `site_id` and `datetime` are
  propagated to new rows from the existing per-sample metadata.

- tolerance:

  Relative tolerance (proportion of measured hardness) for consistency
  check when all three are present. Default `0.05` (5%).

- verbose:

  Logical. If `TRUE`, prints a summary of how many rows were derived and
  how many inconsistencies were detected. Default `TRUE`.

## Value

The input `df` with derived rows appended for `Ca`, `Mg`, or `hardness`
wherever exactly two of the three were available. Derived rows are
tagged with `imputed = TRUE` and `imputed_kind = "derived"` if those
columns exist on the input (otherwise they are added).

## Details

\$\$hardness = 2.497 \cdot Ca + 4.118 \cdot Mg\$\$

where Ca and Mg are in mg/L. This helper fills in any missing member of
{Ca, Mg, hardness} when the other two are available, and warns if all
three are present but inconsistent.

**Per-sample logic:**

- All three present:

  Check consistency. If \\\|\hat{hardness} - hardness\| / hardness \>
  tolerance\\, emit a per-sample warning and keep the user-supplied
  values.

- Exactly two present:

  Compute the third exactly from stoichiometry and append as a new row
  with `detected = TRUE`.

- One or none present:

  Leave alone — fill via imputation if needed.

**Recommended pipeline use:**

Call twice — once **before** imputation to fill samples where the third
member can be derived from raw measurements, and again **after**
[`impute_coanalytes()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_coanalytes.md)
to fill hardness for samples whose Ca and Mg were just imputed.
Idempotent if all three are already consistent.

## Examples

``` r
# Derive hardness from Ca + Mg where it is not measured directly.
chem <- subset(leachate_demo(),
               site_id == "downstream" & analyte %in% c("Ca", "Mg"))
out <- derive_hardness(chem)
#> ℹ derive_hardness: 6 rows derived (hardness=6); 0 samples flagged inconsistent.
subset(out, analyte == "hardness")
#> # A tibble: 6 × 11
#>   sample_id site_id    datetime   analyte  value detected units.analyte
#>   <chr>     <chr>      <date>     <chr>    <dbl> <lgl>    <chr>        
#> 1 DS-01     downstream 2024-01-15 hardness  295. TRUE     NA           
#> 2 DS-02     downstream 2024-03-15 hardness  302. TRUE     NA           
#> 3 DS-03     downstream 2024-05-15 hardness  306. TRUE     NA           
#> 4 DS-04     downstream 2024-07-15 hardness  300. TRUE     NA           
#> 5 DS-05     downstream 2024-09-15 hardness  296. TRUE     NA           
#> 6 DS-06     downstream 2024-11-15 hardness  303. TRUE     NA           
#> # ℹ 4 more variables: valence.analyte <dbl>, atomic_mass.analyte <dbl>,
#> #   imputed <lgl>, imputed_kind <chr>
```
