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

  Numeric. CV of the hardness measurement, used to weight the three
  NO3-N hardness-class SSDs probabilistically (a log-normal blend that
  smooths the soft/moderate/hard boundaries). Default 0.05 (5%). Set to
  0 to recover hard class cutoffs.

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
# Fraction of species affected by 9.3 mg/L total ammonia-N at the pH 7.0 /
# 20 degC index condition. Uses the package's bundled SSD data; set
# options(leachatetools.guideline_dir=) to fit from the ANZG XLSX files.
ssd_paf("NH3-N", conc_ug_L = 9321)
#> $analyte
#> [1] "NH3-N"
#> 
#> $conc_ug_L
#> [1] 9321
#> 
#> $method
#> [1] "multi"
#> 
#> $pct
#> [1] 49.99988
#> 
#> $lower
#> [1] NA
#> 
#> $upper
#> [1] NA
#> 
#> $note
#> character(0)
#> 

# NO3-N needs hardness for automatic soft/moderate/hard SSD selection:
ssd_paf("NO3-N", conc_ug_L = 50000, hardness_mg_L = 90)
#> $analyte
#> [1] "NO3-N"
#> 
#> $conc_ug_L
#> [1] 50000
#> 
#> $method
#> [1] "multi"
#> 
#> $pct
#> [1] 74.51201
#> 
#> $lower
#> [1] NA
#> 
#> $upper
#> [1] NA
#> 
#> $note
#> [1] "NO3-N hardness blend (soft/mod/hard = 0.00/1.00/0.00)"
#> 
```
