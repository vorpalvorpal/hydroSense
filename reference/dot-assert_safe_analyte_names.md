# Guard against make.names() collisions in a set of analyte names

The fitting/prediction path maps original analyte names to syntactically
valid R names via
[`make.names()`](https://rdrr.io/r/base/make.names.html) and relies on
that mapping being a bijection. If two distinct analytes collide under
[`make.names()`](https://rdrr.io/r/base/make.names.html) (e.g. `"Cr-6"`
and `"Cr.6"` both become `"Cr.6"`), a value-equality inverse lookup can
silently return the wrong original name. This guard errors early and
informatively instead.

## Usage

``` r
.assert_safe_analyte_names(targets)
```

## Arguments

- targets:

  Character vector of original analyte names.
