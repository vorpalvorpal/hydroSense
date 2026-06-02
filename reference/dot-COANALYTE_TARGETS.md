# Co-analytes required by ANZECC/ANZG metal normalisation formulas

These are imputed separately (after metals/organics imputation) via
[`impute_coanalytes()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_coanalytes.md)
so that
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md)
has values to normalise against. pH and EC are excluded — they are
always present (required vars).

## Usage

``` r
.COANALYTE_TARGETS
```

## Format

An object of class `character` of length 4.
