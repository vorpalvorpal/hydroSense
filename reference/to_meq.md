# Convert analyte concentrations to milliequivalents per litre

Appends meq/L-converted rows to the input dataframe. Converted rows have
`name.analyte` suffixed with `"_"` to distinguish them from the original
rows, which are preserved unchanged.

## Usage

``` r
to_meq(df)
```

## Arguments

- df:

  Long-format dataframe. Required columns: `analyte`, `value`,
  `units.analyte`, `valence.analyte`, `atomic_mass.analyte`.

## Value

The input `df` with converted rows appended via `bind_rows`. Converted
rows have `analyte` suffixed with `"_"` and `value` in meq/L. Rows
lacking `valence.analyte` or `atomic_mass.analyte` are excluded from
conversion but retained in the original rows.

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
```
