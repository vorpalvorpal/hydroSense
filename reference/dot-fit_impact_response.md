# Fit the season-blind impact response for one analyte

Selects the short/long reservoir recession constants (tau) by profiled
AIC, fits `I ~ s(hydro_short) + s(hydro_long)` (no seasonal term), and
keeps it only if it beats the intercept-only null. Stores the de-trended
residual `S = I - fitted` as the bridge-interpolation state.

## Usage

``` r
.fit_impact_response(
  obs,
  hydro,
  hydro_type,
  tau_bounds_short,
  tau_bounds_long,
  auto_select,
  min_obs_model,
  eps
)
```

## Value

List with `impact_fit`, `tau_short`, `tau_long`, `tier`, `n_obs`,
`anchors` (`date`, `I`, `S`, `hydro_short`, `hydro_long`).
