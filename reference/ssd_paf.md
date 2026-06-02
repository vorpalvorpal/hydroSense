# Estimate the fraction of species potentially affected at a concentration.

Estimate the fraction of species potentially affected at a
concentration.

## Usage

``` r
ssd_paf(
  analyte,
  conc_ug_L,
  method = c("multi", "anzecc"),
  hardness_mg_L = NULL,
  hardness_cv = 0.05,
  guideline_dir = getOption("leachatetools.guideline_dir"),
  nboot = 0L,
  level = 0.95
)
```

## Arguments

- analyte:

  Character. Analyte name (key in .SSD_NAME_MAP). Supply "NO3-N"
  together with `hardness_mg_L` for automatic class selection, or supply
  the explicit class name ("NO3-N_soft" etc.).

- conc_ug_L:

  Numeric. Concentration in µg/L (after any external physicochemical
  corrections).

- method:

  Character. "multi" (default) fits all 6 BCANZ distributions and
  model-averages; "anzecc" uses the per-analyte distribution that best
  matches the original ANZG derivation.

- hardness_mg_L:

  Numeric or NULL. Required for NO3-N analyte. Hardness in mg/L CaCO3 at
  the time of measurement.

- hardness_cv:

  Numeric. CV of hardness measurement for probabilistic class weighting.
  Currently reserved — hard cutoffs used. Default 0.05.

- guideline_dir:

  Character. Path to the "guideline data" folder containing ANZG XLSX
  files. Falls back to getOption("leachatetools.guideline_dir").

- nboot:

  Integer. Bootstrap replicates for CI. 0 = no CI.

- level:

  Numeric. Confidence level (default 0.95).

## Value

Named list: \$analyte character \$conc_ug_L numeric \$method character
\$pct numeric — % species affected, NA if no model \$lower numeric —
lower CI % \$upper numeric — upper CI % \$note character vector —
caveats

## Examples

``` r
if (FALSE) { # \dontrun{
# Fraction of species affected by 9.3 mg/L total ammonia-N at the
# pH 7.0 / 20 degC index condition (requires the ANZG guideline data):
options(leachatetools.guideline_dir = "path/to/guideline data")
ssd_paf("NH3-N", conc_ug_L = 9321)

# NO3-N needs hardness for automatic soft/moderate/hard SSD selection:
ssd_paf("NO3-N", conc_ug_L = 50000, hardness_mg_L = 90)
} # }
```
