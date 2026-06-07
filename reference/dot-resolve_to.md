# Resolve a bare numeric or units object to a target unit

Resolve a bare numeric or units object to a target unit

## Usage

``` r
.resolve_to(x, target_unit, units_str = NULL, arg_name = "x")
```

## Arguments

- x:

  Numeric or units object.

- target_unit:

  Character udunits2 unit string (desired output).

- units_str:

  Companion unit string; required when `x` is bare numeric.

- arg_name:

  Name of the calling argument, used in error messages.

## Value

Bare numeric in `target_unit`.
