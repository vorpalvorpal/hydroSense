# Pooled (hierarchical) season-blind impact response across analytes

Fits a single factor-smooth GAM
`z ~ s(hydro_short, analyte, bs="fs") + s(hydro_long, analyte, bs="fs")`
at one AIC-selected common window pair, where `z = (I - mu_a) / sd_a` is
the impact standardised **per analyte**. Pooling on this common unit
scale shares the response *shape* across analytes while leaving each
analyte's *magnitude* untouched (it is restored by de-standardising the
fitted shape, `fitted_I = mu_a + sd_a * z_hat`). This is essential: the
`bs = "fs"` penalty shrinks each analyte's level as well as its
wiggliness, so pooling the raw `I` would drag a large-signal analyte
(e.g. Cu) toward a population dominated by near-zero ones and inflate
the near-zero ones in turn.

## Usage

``` r
.fit_pooled_impact_response(
  anchors_all,
  target_analytes,
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

Named list (per analyte) of
[`.fit_impact_response()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-fit_impact_response.md)-shaped
objects; pooled analytes additionally carry `pooled = TRUE`, `analyte`,
`pool_levels`, and the de-standardisation pair `pool_center` (`mu_a`)
and `pool_scale` (`sd_a`) for prediction.

## Details

Pooling regularises noisy, low-SNR analytes by borrowing a response
shape from co-varying ones – it does **not** add hydrological coverage
(co-sampled analytes share the same regimes). Analytes with fewer than
`min_obs_model` anchors, or with ~no impact variance (no shape to
share), get a per-analyte flat bridge; if pooling fails or doesn't beat
the no-shared-shape null (`z ~ 1`), the poolable analytes fall back to
independent fits.
