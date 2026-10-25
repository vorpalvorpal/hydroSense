# Analyte-name aliases

The mapping from the analyte names used by sampleTidy (and other
upstream stores) to the names hydroSense computes with.
[`add_lmf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_lmf.md),
[`to_meq()`](https://vorpalvorpal.github.io/hydroSense/reference/to_meq.md),
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)
and
[`mspaf_pipeline()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_pipeline.md)
apply it to their inputs, so data stored under either spelling is
accepted, and return the caller's rows with their original names.

## Usage

``` r
analyte_aliases()
```

## Value

A data frame with character columns `alias` (an upstream name) and
`canonical` (the hydroSense name), one row per alias. Names that already
match map to themselves.

## Details

Two upstream quantities have no hydroSense counterpart of their own and
map to themselves: `"Alkalinity-total-CaCO3"` (the sum of `"CO3-CaCO3"`
and `"HCO3-CaCO3"`) and `"NO2+NO3-N"` (the sum of `"NO3-N"` and
`"NO2-N"`).
[`add_lmf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_lmf.md)
uses such a combined value in place of its two species, and ignores it
when the species themselves are present, so nothing is counted twice.

## Examples

``` r
analyte_aliases()
#> # A tibble: 21 × 2
#>    alias canonical
#>    <chr> <chr>    
#>  1 Na    Na       
#>  2 K     K        
#>  3 Ca    Ca       
#>  4 Mg    Mg       
#>  5 Cl    Cl       
#>  6 SO4²⁻ SO4²⁻    
#>  7 F     F        
#>  8 NH3-N NH3-N    
#>  9 NO3-N NO3-N    
#> 10 NO2-N NO2-N    
#> # ℹ 11 more rows
```
