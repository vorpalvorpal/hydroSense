# Quick point-estimate only (no CI, no bootstrapping).

Quick point-estimate only (no CI, no bootstrapping).

## Usage

``` r
ssd_pct(
  analyte,
  conc_ug_L,
  method = "multi",
  hardness_mg_L = NULL,
  guideline_dir = getOption("leachatetools.guideline_dir")
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

- guideline_dir:

  Character. Path to the "guideline data" folder containing ANZG XLSX
  files. Falls back to getOption("leachatetools.guideline_dir").

## Value

Numeric — % species affected, or NA.

## Examples

``` r
ssd_pct("Zn", conc_ug_L = 30)
#> [1] 42.01064
```
