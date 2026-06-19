# Run the full daily msPAF pipeline end to end

Orchestrates the three steps a daily msPAF analysis usually chains by
hand: optionally fit a Bayesian imputation model
([`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md)),
fit the reference/impact model
([`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md)),
then compute the daily multi-substance PAF
([`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md))
for the target site. The fitted imputation model is threaded into
**both** downstream steps so that below-detection and entirely-absent
analytes are filled before the toxicity calculation.

## Usage

``` r
mspaf_pipeline(
  target,
  reference,
  hydro = NULL,
  impute = TRUE,
  impute_on = c("reference", "target"),
  required_vars = c("pH", "EC"),
  impute_groups = NULL,
  impute_seed = NULL,
  reference_args = list(),
  daily_args = list()
)
```

## Arguments

- target:

  Long-format target-site chemistry (same schema as
  [`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)'s
  `df`): `sample_id`, `site_id`, `datetime`, `analyte`, `value`,
  `detected`.

- reference:

  Long-format reference-site chemistry passed to
  [`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md).

- hydro:

  Optional hydrology series forwarded to
  [`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md).

- impute:

  Logical. Fit and apply a Bayesian imputation model? Default `TRUE`.
  When `FALSE`, no imputation model is fitted and `NULL` is passed
  downstream (identical to chaining the functions by hand without one).

- impute_on:

  Which chemistry to fit the imputation model on: `"reference"`
  (default, matching the documented design) or `"target"`.

- required_vars:

  Passed to
  [`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md).
  Default `c("pH", "EC")`.

- impute_groups:

  Optional `groups` for
  [`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md);
  `NULL` (default) uses
  [`leachate_impute_groups()`](https://vorpalvorpal.github.io/hydroSense/reference/leachate_impute_groups.md).

- impute_seed:

  Optional integer seed for the imputation fit (forwarded to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) via
  [`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md))
  so the default-on path is reproducible.

- reference_args, daily_args:

  Named lists of additional arguments forwarded to
  [`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md)
  and
  [`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)
  respectively. Do not include arguments the orchestrator sets itself
  (`reference`/`hydro`/ `imputation_model` for the reference fit;
  `df`/`reference_model`/ `imputation_model` for the daily call).
  `daily_args` may override the `interpolation = "model"` default.

## Value

The
[`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)
result, with the fitted models attached as attributes
`"reference_model"` and `"imputation_model"` (the latter `NULL` when
`impute = FALSE` or the fit was skipped). Note that most dplyr verbs
drop attributes, so read them off the returned frame directly.

## Details

**Imputation is on by default here.** The individual functions treat
imputation as opt-in (you must fit a model and pass it). This
orchestrator flips that default: with `impute = TRUE` it fits an
imputation model and applies it throughout. This is convenient but **not
free of consequences** — imputing below-detection cells (replacing
detection-limit values with modelled sub-limit estimates) and
fabricating entirely-absent analytes generally *moves* the msPAF result.
Whether that is appropriate is a missingness (MAR/MNAR) modelling
judgement, not a neutral default. Set `impute = FALSE` to reproduce the
non-imputed behaviour of calling
[`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md)
and
[`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)
directly.

## See also

[`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md),
[`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md),
[`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)

## Examples

``` r
if (FALSE) { # \dontrun{
demo <- leachate_demo()
out <- mspaf_pipeline(
  target = subset(demo, site_id == "downstream"),
  reference = subset(demo, site_id == "reference"),
  daily_args = list(require_temperature = FALSE, ndraws = 50L, seed = 1L)
)
} # }
```
