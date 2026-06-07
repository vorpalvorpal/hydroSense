# Leachate imputation-group preset

Returns the default
[`impute_group()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_group.md)
specification used by
[`fit_imputation_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_imputation_model.md)
for landfill-leachate monitoring chemistry: a **metals** group (hurdled
on metal presence) and a catch-all **organics** group (hurdled on
dissolved-organic-carbon presence). This is the leachate domain layer on
top of the otherwise domain-agnostic engine — pass your own list of
[`impute_group()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_group.md)
objects to model a different chemistry.

## Usage

``` r
leachate_impute_groups()
```

## Value

A list of two `"impute_group"` objects (`metals`, `organics`).

## Details

The metals set is
[.METAL_ANALYTES](https://vorpalvorpal.github.io/leachatetools/reference/dot-METAL_ANALYTES.md)
and the organic-carbon hurdle set is
[.DOC_LIKE_ANALYTES](https://vorpalvorpal.github.io/leachatetools/reference/dot-DOC_LIKE_ANALYTES.md).

## See also

[`impute_group()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_group.md),
[`fit_imputation_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_imputation_model.md)

## Examples

``` r
leachate_impute_groups()
#> [[1]]
#> <impute_group: metals>
#>   targets: Al, As, B, Ba, Be, Cd, Co, Cr, Cr-6, Cu, Fe, Hg, Mn, Mo, Ni, Pb, Sb, Se, Sn, Sr, V, Zn
#>   hurdle:  Al, As, B, Ba, Be, Cd, Co, Cr, Cr-6, Cu, Fe, Hg, Mn, Mo, Ni, Pb, Sb, Se, Sn, Sr, V, Zn
#> 
#> [[2]]
#> <impute_group: organics>
#>   targets: <catch-all>
#>   hurdle:  DOC, TOC, BOD, COD, cBOD
#> 
```
