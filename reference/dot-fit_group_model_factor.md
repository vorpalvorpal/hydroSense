# Fit the Route C two-stage censored factor group model

Stage 1: per-analyte `mgcv::gam(lv ~ s(PC1) + ...)` on detected
observations gives the spline mean `mu_j(X_i)`; detected residuals and
BDL censoring bounds (`log(DL) - mu_j`) are computed from it. Stage 2: a
Stan program (`inst/stan/factor_censored.stan`) fits a rank-`k` factor
model on those residuals, with BDL cells as upper-bounded latent
parameters (proper left-censoring, jointly with the factor and
covariance).

## Usage

``` r
.fit_group_model_factor(
  target_analytes,
  safe_analytes,
  base,
  pc_wide,
  pc_cols,
  log_floors,
  iter,
  warmup,
  chains,
  cores,
  group_name = "group",
  k = NULL,
  seed = NULL,
  ...
)
```
