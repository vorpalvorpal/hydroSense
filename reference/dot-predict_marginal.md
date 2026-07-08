# Posterior-predictive prediction for the marginal method

Per analyte, per cell that needs imputing (BDL or entirely absent): draw
the GAM mean with parameter uncertainty (`beta ~ N(coef, Vp)`, `Vp`
unconditional so smoothing-parameter uncertainty is included) plus
residual noise `N(0, sig2)`; BDL cells draw the residual left-truncated
so the value never exceeds `log(DL)`. This is a proper posterior
predictive (mean + residual uncertainty), not the plug-in residual-only
draw. Returns the same `pm_long` shape as the other predictors so
[`.predict_and_merge()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-predict_and_merge.md)
is unchanged.

## Usage

``` r
.predict_marginal(group, pc_wide, df_eligible, return, ndraws, batch_size)
```
