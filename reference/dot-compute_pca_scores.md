# Project new data onto stored NIPALS PCA axes

Handles within-row missing values via NIPALS regression scoring: each
component score is estimated from observed variables only, then the
residual is deflated before the next component. Columns entirely absent
in `df` (not measured at all, not just BDL) are filled with training
medians.

## Usage

``` r
.compute_pca_scores(df, pca_obj)
```
