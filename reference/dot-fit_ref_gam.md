# Fit a GAM on reference observations for one analyte

Tries [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) with cyclic
spline on day-of-year and thin-plate splines on the two hydro features.
Returns the fitted `gam` object, or `NULL` if fitting fails.

## Usage

``` r
.fit_ref_gam(df_model, eps = 1e-09)
```

## Arguments

- df_model:

  Tibble `(y, doy, hydro_short, hydro_long)`.

- eps:

  Guard already applied to `y` upstream (unused here, for docs).

## Value

A `gam` object or `NULL`.
