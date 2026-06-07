# Quick point-estimate only (no CI, no bootstrapping).

Quick point-estimate only (no CI, no bootstrapping).

## Usage

``` r
ssd_pct(
  analyte,
  conc,
  conc_units = NULL,
  method = "multi",
  hardness = NULL,
  hardness_units = NULL,
  guideline_dir = getOption("leachatetools.guideline_dir")
)
```

## Arguments

- analyte:

  Character. Analyte name (key in .SSD_NAME_MAP). Supply "NO3-N"
  together with `hardness` for automatic class selection, or supply the
  explicit class name ("NO3-N_soft" etc.).

- conc:

  Numeric or `units` object. Concentration to evaluate. Bare numeric
  requires `conc_units`.

- conc_units:

  Character. Unit of `conc` when it is bare numeric, e.g. `"ug/L"` or
  `"mg/L"`. Ignored when `conc` is a `units` object.

- method:

  Character. "multi" (default) fits all 6 BCANZ distributions and
  model-averages; "anzecc" uses the per-analyte distribution that best
  matches the original ANZG derivation.

- hardness:

  Numeric or `units` object, or `NULL`. Required for the NO3-N analyte.
  Hardness at the time of measurement. Bare numeric requires
  `hardness_units`.

- hardness_units:

  Character. Unit of `hardness` when it is bare numeric, e.g. `"mg/L"`.
  Ignored when `hardness` is a `units` object or `NULL`.

- guideline_dir:

  Character. Path to the "guideline data" folder containing ANZG XLSX
  files. Falls back to getOption("leachatetools.guideline_dir").

## Value

Numeric — % species affected, or NA.

## Examples

``` r
ssd_pct("Zn", conc = 30, conc_units = "ug/L")
#> [1] 42.01064
```
