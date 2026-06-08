# Pooled (hierarchical) season-blind impact response across analytes

Fits a single factor-smooth GAM
`I ~ s(hydro_short, analyte, bs="fs") + s(hydro_long, analyte, bs="fs")`
at one AIC-selected common window pair, so each analyte's hydrological
response is shrunk toward a shared population response (partial
pooling). This regularises noisy, low-SNR analytes by borrowing a
response shape from co-varying ones — it does **not** add hydrological
coverage (co-sampled analytes share the same regimes). Analytes with
fewer than `min_obs_model` anchors get a per-analyte flat bridge; if
pooling fails or doesn't beat an analyte-intercept null, the poolable
analytes fall back to independent fits.

## Usage

``` r
.fit_pooled_impact_response(
  anchors_all,
  target_analytes,
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

Named list (per analyte) of
[`.fit_impact_response()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-fit_impact_response.md)-shaped
objects; pooled analytes additionally carry `pooled = TRUE`, `analyte`,
and `pool_levels` for prediction.
