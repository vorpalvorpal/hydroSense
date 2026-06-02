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
[`impute_coanalytes()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_coanalytes.md)
to fill hardness for samples whose Ca and Mg were just imputed.
Idempotent if all three are already consistent.

## Examples

``` r
if (FALSE) { # \dontrun{
# Fill missing hardness where Ca and Mg are both measured
chem2 <- derive_hardness(chem)

# Call twice in a typical pipeline: pre- and post-imputation
chem      <- derive_hardness(chem)
chem_imp  <- impute_chemistry(chem, model)
chem_imp2 <- impute_coanalytes(chem_imp, model)
chem_imp3 <- derive_hardness(chem_imp2)
} # }
```
