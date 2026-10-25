# Convert analyte concentrations to milliequivalents per litre

Appends meq/L-converted rows to the input dataframe. Converted rows have
`analyte` suffixed with `"_"` to distinguish them from the original
rows, which are preserved unchanged.

## Usage

``` r
to_meq(df)
```

## Arguments

- df:

  Long-format dataframe. Required columns: `analyte`, `value`, and a
  unit column: `units.analyte` or sampleTidy's `units` (a unit string
  per row, e.g. `"mg/L"`, `"ug/L"`). Optional: `valence.analyte`,
  `atomic_mass.analyte`.

## Value

The input `df` with converted rows appended via `bind_rows`. Converted
rows have `analyte` suffixed with `"_"` and `value` in meq/L. Rows with
no known valence and molar mass (not a panel ion and none supplied) are
excluded from conversion but retained in the original rows.

## Details

Valence and molar mass for the LMF panel ions (Na, K, Ca, Mg, Cl, F,
SO4, NH3-N, NO3-N, NO2-N, CO3-CaCO3, HCO3-CaCO3) are built in, so a
frame with only `analyte`, `value` and a unit column converts. Nitrogen
species are taken as N and the carbonate species as CaCO3. A
`valence.analyte` / `atomic_mass.analyte` supplied by the caller takes
precedence for the rows it is non-`NA` on.

## Examples

``` r
df <- tibble::tibble(
  analyte             = "Ca",
  value               = 40.08,   # mg/L
  units.analyte       = "mg/L",
  valence.analyte     = 2,
  atomic_mass.analyte = 40.08    # g/mol
)
to_meq(df)   # appends a "Ca_" row of 2 meq/L (1 mmol/L x |valence| = 2)
#> # A tibble: 2 × 5
#>   analyte value units.analyte valence.analyte atomic_mass.analyte
#>   <chr>   <dbl> <chr>                   <dbl>               <dbl>
#> 1 Ca       40.1 mg/L                        2                40.1
#> 2 Ca_       2   mg/L                        2                40.1

# Panel ions need no valence or molar mass columns.
to_meq(tibble::tibble(analyte = "Ca", value = 40.08, units = "mg/L"))
#> # A tibble: 2 × 3
#>   analyte value units
#>   <chr>   <dbl> <chr>
#> 1 Ca      40.1  mg/L 
#> 2 Ca_      2.00 mg/L 
```
