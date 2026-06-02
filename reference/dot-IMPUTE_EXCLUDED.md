# Analytes that must never enter the imputation model as response variables

These are excluded from both the metals and organics groups. They are
typically non-concentration measurements (counts, qualitative, physical)
for which a log-normal concentration model is inappropriate.

## Usage

``` r
.IMPUTE_EXCLUDED
```

## Format

An object of class `character` of length 10.
