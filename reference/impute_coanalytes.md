# Impute missing normalisation co-analytes from the fitted chemistry PCA

Fits a univariate log-Gaussian GAM
([`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html)) for each target
co-analyte using the PC scores already computed by
[`fit_imputation_model()`](https://www.kedumba.com.au/leachatetools/reference/fit_imputation_model.md).
Only samples where the co-analyte is entirely absent are filled; BDL
observations are left unchanged.

## Usage

``` r
impute_coanalytes(df, model, targets = NULL, min_obs = 10L)
```

## Arguments

- df:

  Long-format chemistry data frame (same schema as
  [`impute_chemistry()`](https://www.kedumba.com.au/leachatetools/reference/impute_chemistry.md),
  with `imputed`/`imputed_kind` columns if
  [`impute_chemistry()`](https://www.kedumba.com.au/leachatetools/reference/impute_chemistry.md)
  has already been called).

- model:

  Fitted model from
  [`fit_imputation_model()`](https://www.kedumba.com.au/leachatetools/reference/fit_imputation_model.md)
  (provides the PCA object and the list of `pca_vars`).

- targets:

  Co-analyte names to impute when missing. Default `.COANALYTE_TARGETS`
  (`"DOC"`, `"Ca"`, `"Mg"`, `"hardness"`). Only targets present in
  `model$pca_vars` are processed; others are skipped with a warning.

- min_obs:

  Minimum number of quantified observations required to fit a GAM for a
  target. Targets with fewer observations are skipped. Default `10L`.

## Value

`df` with missing co-analyte rows filled in, tagged with
`imputed = TRUE` and `imputed_kind = "missing"`. All other rows are
unchanged.

## Details

This step belongs **after**
[`impute_chemistry()`](https://www.kedumba.com.au/leachatetools/reference/impute_chemistry.md)
and **before**
[`time_weighted_aggregate()`](https://www.kedumba.com.au/leachatetools/reference/time_weighted_aggregate.md).
Imputed co-analyte values are never fed back into the metals/organics
model — the brms model ran on measured values only and is already done.

Using the chemistry PCA as the sole predictor set is appropriate
because: (a) the PCA already captures DOC/Ca/Mg variation in its axes;
(b) a univariate GAM on PC scores is unbiased and fast (no Stan
required); (c) the same PCA is used for the metals model so the
co-analyte predictions are conditioned on the same chemical environment
summary.

## See also

[`fit_imputation_model()`](https://www.kedumba.com.au/leachatetools/reference/fit_imputation_model.md),
[`impute_chemistry()`](https://www.kedumba.com.au/leachatetools/reference/impute_chemistry.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Deterministic GAM-based imputation of normalisation co-analytes
# (pH, DOC, hardness, ...) from the measured analyte suite.
impute_coanalytes(monitoring_long)
} # }
```
