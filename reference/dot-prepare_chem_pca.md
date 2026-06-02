# Fit the unified chemistry PCA on training data

This PCA spans the full unified chemistry predictor set (`pca_vars` in
[`fit_imputation_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_imputation_model.md))
— major ions, pH, EC, NH3-N, DOC, nutrients and redox indicators.
Concentration-like variables are `log10`-transformed (see
[`.log_transform_pca()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-log_transform_pca.md))
before centring/scaling. PC score columns are named `PC1`, `PC2`, ….

## Usage

``` r
.prepare_chem_pca(df, wq_vars, min_var_explained = 0.75, max_pcs = 4L)
```
