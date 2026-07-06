# Scale-aware log floor for a single numeric vector

Half the smallest observed positive value, so a genuine zero maps just
below the column's own real values rather than to a single absolute
constant shared across analytes of wildly different magnitude. Falls
back to `eps` when the vector has no observed positive values at all.

## Usage

``` r
.scale_aware_log_floor(x, eps = 1e-09)
```

## Arguments

- x:

  Numeric vector.

- eps:

  Absolute fallback floor (default `1e-9`).
