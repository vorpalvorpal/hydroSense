# Parse a normalisation formula string into a quoted expression

Parse a normalisation formula string into a quoted expression

## Usage

``` r
.parse_normalisation_formula(formula_str)
```

## Arguments

- formula_str:

  Character string containing an R expression. The expression may
  reference `C` (the raw concentration) and any co-analyte names listed
  in `coanalytes_required`. Empty string or `NA` → identity (returns `C`
  unchanged).

## Value

A language object (quoted expression), or `NULL` for identity.
