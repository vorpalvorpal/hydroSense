# Apply a parsed normalisation formula to a concentration value

Apply a parsed normalisation formula to a concentration value

## Usage

``` r
.apply_normalisation(parsed_expr, C, coanalytes = numeric(0))
```

## Arguments

- parsed_expr:

  Parsed expression from
  [`.parse_normalisation_formula()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-parse_normalisation_formula.md),
  or `NULL` for identity.

- C:

  Numeric concentration (µg/L or relevant unit).

- coanalytes:

  Named numeric vector of co-analyte values available at this sample
  (e.g. `c(DOC = 0.5, pH = 7.2)`). May be empty.

## Value

Normalised concentration. Returns `C` when `parsed_expr` is `NULL`.
Returns `NA_real_` if evaluation fails (e.g. required co-analyte absent
from `coanalytes`).
