# Hydrology-weighted bracketing-anchor interpolation of the residual state

Pinches to the observed residual at each anchor; between two bracketing
anchors it leans toward the one whose antecedent hydrology more closely
matches the query day (and, secondarily, the nearer in time). Outside
the anchor span it carries the nearest anchor's residual
forward/backward (flat).

## Usage

``` r
.interp_residual(anchors, qdate, qshort, qlong)
```

## Arguments

- anchors:

  Tibble with `date`, `S`, `hydro_short`, `hydro_long` (≥1 row).

- qdate:

  Query Date scalar.

- qshort, qlong:

  Query-day hydro features.

## Value

Numeric interpolated residual `S`.
