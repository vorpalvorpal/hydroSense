# Fit the season-blind impact response for one analyte

Selects short/long antecedent windows by AIC, fits
`I ~ s(hydro_short) + s(hydro_long)` (no seasonal term), and keeps it
only if it beats the intercept-only null. Stores the de-trended residual
`S = I - fitted` as the bridge-interpolation state.

## Usage

``` r
.fit_impact_response(
  obs,
  hydro,
  hydro_type,
  api_windows_short,
  api_windows_long,
  auto_select,
  min_obs_model,
  eps
)
```

## Value

List with `impact_fit`, `window_short`, `window_long`, `tier`, `n_obs`,
`anchors` (`date`, `I`, `S`, `hydro_short`, `hydro_long`).
