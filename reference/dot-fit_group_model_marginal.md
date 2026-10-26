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

## Details

Each GAM is fitted by REML on the first `m` PC scores only (see
[`.marginal_gam_rhs()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-marginal_gam_rhs.md)):
the number of smooths is gated on the analyte's detected count so that
thin analytes are not fitted with more flexibility than their data
support. GCV on the full five-smooth model saturated such analytes (edf
equal to the number of points), which made their imputations sensitive
to ulp-level changes in the inputs.
