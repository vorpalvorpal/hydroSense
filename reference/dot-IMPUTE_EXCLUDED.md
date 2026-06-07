# Analytes that must never enter the imputation model as response variables

The default `exclude` set for
[`fit_imputation_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_imputation_model.md)
(a leachate-preset default; override via the `exclude` argument for
other domains). These are typically non-concentration measurements
(counts, qualitative, physical) for which a log-normal concentration
model is inappropriate, so they are excluded from every imputation
group.

## Usage

``` r
.IMPUTE_EXCLUDED
```

## Format

An object of class `character` of length 10.
