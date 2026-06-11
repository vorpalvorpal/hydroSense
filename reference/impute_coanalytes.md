# Impute missing normalisation co-analytes from the fitted chemistry PCA

Fits a univariate log-Gaussian GAM
([`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html)) for each target
co-analyte using the PC scores already computed by
[`fit_imputation_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_imputation_model.md).
Only samples where the co-analyte is entirely absent are filled; BDL
observations are left unchanged.

## Usage

``` r
impute_coanalytes(
  df,
  model,
  targets = NULL,
  min_obs = 10L,
  return = c("point", "draws"),
  ndraws = NULL,
  seed = NULL
)
```

## Arguments

- df:

  Long-format chemistry data frame (same schema as
  [`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md),
  with `imputed`/`imputed_kind` columns if
  [`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md)
  has already been called).

- model:

  Fitted model from
  [`fit_imputation_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_imputation_model.md)
  (provides the PCA object and the list of `pca_vars`).

- targets:

  Co-analyte names to impute when missing. Default `.COANALYTE_TARGETS`
  (`"DOC"`, `"Ca"`, `"Mg"`, `"hardness"`). Only targets present in
  `model$pca_vars` are processed; others are skipped with a warning.

- min_obs:

  Minimum number of quantified observations required to fit a GAM for a
  target. Targets with fewer observations are skipped. Default `10L`.

- return:

  `"point"` (default) for the posterior mean of the GAM prediction —
  identical to pre-draws behaviour. `"draws"` for full
  posterior-predictive draws: each missing co-analyte cell emits `N`
  rows keyed by `draw_id 1..N`, reflecting both parameter uncertainty
  (`beta ~ N(coef(gam), Vp)`) and residual Gaussian noise (`gam$sig2`).
  Observed co-analyte cells keep `draw_id = NA`.

- ndraws:

  Number of draws to generate. Required when `return = "draws"` and `df`
  contains no existing draws. When `df` already carries draws (from
  [`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md)),
  `N` is inferred from the existing draw domain; `ndraws` must be `NULL`
  or equal to that count.

- seed:

  Optional integer seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before the sampling
  calls, for reproducibility.

## Value

`df` with missing co-analyte rows filled in, tagged with
`imputed = TRUE` and `imputed_kind = "missing"`. In `"draws"` mode each
imputed cell is replicated `N` times with `draw_id 1..N`; observed cells
keep `draw_id = NA`. In `"point"` mode the output schema is unchanged
from the pre-draws behaviour.

## Details

This step belongs **after**
[`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md)
and **before**
[`time_weighted_aggregate()`](https://vorpalvorpal.github.io/leachatetools/reference/time_weighted_aggregate.md).
Imputed co-analyte values are never fed back into the metals/organics
model — the brms model ran on measured values only and is already done.

Using the chemistry PCA as the sole predictor set is appropriate
because: (a) the PCA already captures DOC/Ca/Mg variation in its axes;
(b) a univariate GAM on PC scores is unbiased and fast (no Stan
required); (c) the same PCA is used for the metals model so the
co-analyte predictions are conditioned on the same chemical environment
summary.

## See also

[`fit_imputation_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_imputation_model.md),
[`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md),
[`summarise_draws()`](https://vorpalvorpal.github.io/leachatetools/reference/summarise_draws.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Deterministic GAM-based imputation (default, point mode)
impute_coanalytes(monitoring_long, model)

# Posterior-predictive draws when df already carries metals draws
impute_coanalytes(metals_draws, model, return = "draws")
} # }
```
