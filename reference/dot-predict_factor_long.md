# Long-format prediction for the cens_factor shared-latent-factor model

Builds one row per (eligible sample × analyte), predicts the latent log
concentration from the univariate model (the per-sample
`(1 | sample_id)` factor couples analytes), and returns the same
`pm_long` shape the wide path produces so the merge step is identical.

## Usage

``` r
.predict_factor_long(group, pc_wide, return, ndraws, batch_size)
```
