# Fit the marginal (no-borrowing) group model

One left-censored GAM per target analyte (`lv ~ s(PC1) + ...` on
detected observations), with NO shared factor. Uses mgcv only – no brms,
no cmdstanr. The "borrowing" the factor model attempts is skipped
entirely, so this is the robust choice on panels where cross-analyte
correlation is weak or ragged (sparse analytes' spurious loadings can't
mis-condition anything).

## Usage

``` r
.fit_group_model_marginal(
  target_analytes,
  safe_analytes,
  base,
  pc_wide,
  pc_cols,
  log_floors,
  group_name = "group"
)
```
