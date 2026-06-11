# Fit the WQ-\>metal layer and its residual `d` for one analyte

A GAM of normalised concentration on the chemistry-PCA scores (the
non-Bayesian WQ prediction), kept only if it beats an intercept-only
null by AIC. The residual `d = value_norm - WQ-prediction` is
date-aggregated with hydro features for bracketing-bridge interpolation,
exactly like the impact state `S`.

## Usage

``` r
.fit_wq_layer(
  os_nm,
  pc_cols,
  hydro,
  hydro_type,
  window_short,
  window_long,
  min_obs_model
)
```

## Value

List with `wq_fit` (gam or `NULL`) and `d_anchors` (or `NULL`).
