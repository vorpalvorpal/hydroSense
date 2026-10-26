# Right-hand side of a marginal Stage-1 GAM, gated on the detected count

Uses the first `m = min(P, floor(n_det / 3) - 1)` PC scores, one smooth
each, so the unpenalised part of the model (intercept plus one linear
term per smooth) has at most `n_det / 3` coefficients. `m = 0` (fewer
than six detected values) gives the intercept-only model `lv ~ 1`. Each
smooth's basis dimension is `k = min(10, floor((n_det - 3) / m) + 1)`,
mgcv's default `k = 10` whenever there are enough points. The model then
has at most `1 + m * (k - 1) <= n_det - 2` coefficients, so its edf is
at most `n_det - 2` by construction (REML on a larger basis can still
interpolate a thin analyte), and `k >= 4` always, since
`m <= n_det / 3 - 1`.

## Usage

``` r
.marginal_gam_rhs(pc_cols, n_det)
```

## Arguments

- pc_cols:

  Character vector of PC score columns, in order (`PC1`, ...).

- n_det:

  Number of detected observations the GAM is fitted to.

## Value

A single string, e.g. `"s(PC1, k = 10) + s(PC2, k = 10)"` or `"1"`.
