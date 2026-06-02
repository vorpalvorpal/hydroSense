# Return HC50 for use in Concentration Addition (msPAF).

HC50 is the concentration at which 50% of species are predicted to be
affected. Used as the denominator when computing Toxic Units for the
Concentration Addition combination step of msPAF.

## Usage

``` r
ssd_hc50(
  analyte,
  method = c("multi", "anzecc"),
  hardness_mg_L = NULL,
  guideline_dir = getOption("leachatetools.guideline_dir")
)
```

## Arguments

- analyte:

  Character. Analyte name (key in .SSD_NAME_MAP). Supply "NO3-N"
  together with `hardness_mg_L` for automatic class selection, or supply
  the explicit class name ("NO3-N_soft" etc.).

- method:

  Character. "multi" (default) fits all 6 BCANZ distributions and
  model-averages; "anzecc" uses the per-analyte distribution that best
  matches the original ANZG derivation.

- hardness_mg_L:

  Numeric or NULL. Required for NO3-N analyte. Hardness in mg/L CaCO3 at
  the time of measurement.

- guideline_dir:

  Character. Path to the "guideline data" folder containing ANZG XLSX
  files. Falls back to getOption("leachatetools.guideline_dir").

## Value

Numeric scalar — HC50 in µg/L, or NA.

## Examples

``` r
# HC50 (50% of species affected) for copper, used as the
# Toxic-Unit denominator in the msPAF concentration-addition step:
ssd_hc50("Cu")
#> [1] 4.225759
```
