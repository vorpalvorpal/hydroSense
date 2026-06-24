# Run the full chronic msPAF pipeline end to end

Orchestrates the four steps a chronic msPAF analysis usually chains by
hand: optionally fit a Bayesian imputation model
([`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md)),
fit the reference/impact model
([`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md)),
compute the daily multi-substance PAF
([`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md))
for the target site, then aggregate the daily series into a **chronic**
time-weighted msPAF
([`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md)).
The fitted imputation model is threaded into **both** downstream
chemistry steps so that below-detection and entirely-absent analytes are
filled before the toxicity calculation.

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
  daily_args = list(),
  chronic = TRUE,
  focal_dates = NULL,
  focal_by = "day",
  tau = NULL,
  window = NULL,
  chronic_summary = "arith_mean"
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
  `daily_args` may override the `interpolation = "model"` default. When
  `chronic = TRUE` and `daily_args` requests draws (`ndraws`),
  `return = "draws"` is forced so per-draw rows can be propagated into
  the chronic summary; the `gap_uncertainty`, `interval` and `central`
  values from `daily_args` (or
  [`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)'s
  defaults `"bracket"`, `0.9`, `"median"`) shape the chronic envelope.

- chronic:

  Logical. Aggregate the daily series into chronic time-weighted msPAF
  as a final step? Default `TRUE`. When `FALSE`, the daily
  [`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)
  result is returned unchanged (with model attributes only).

- focal_dates:

  Focal dates for the chronic aggregation, passed to
  [`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md).
  `NULL` (default) derives a sequence spanning the daily date range via
  [`expand_focal_dates()`](https://vorpalvorpal.github.io/hydroSense/reference/expand_focal_dates.md)
  (spacing `focal_by`). May also be a `Date` vector or a
  `focal_date`/`site_id` data frame.

- focal_by:

  Spacing for the derived `focal_dates` when `focal_dates` is `NULL`;
  passed to
  [`expand_focal_dates()`](https://vorpalvorpal.github.io/hydroSense/reference/expand_focal_dates.md)'s
  `by` ([`base::seq.Date()`](https://rdrr.io/r/base/seq.Date.html)
  semantics): `"day"` (default), `7` (weekly), `"week"`, `"month"`, etc.
  Ignored when `focal_dates` is supplied.

- tau:

  Exponential-decay parameter in days for the chronic kernel, forwarded
  to
  [`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md).
  `NULL` (default) uses its default of 90 days.

- window:

  Look-back window in days for the chronic aggregation, forwarded to
  [`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md).
  `NULL` (default) uses its default of 365 days.

- chronic_summary:

  Aggregation method for the chronic step, passed to
  [`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md)'s
  `summary`. Default `"arith_mean"` (the recommended aggregation for
  bounded msPAF percentages).

## Value

When `chronic = TRUE` (default), the chronic msPAF frame, one row per
(focal date × site):

- **point mode** (`ndraws` absent): columns `focal_date`, `site_id`,
  `sample_id`, `analyte`, `value`, `detected`, `n_samples_in_window`,
  `n_imputed_in_window` (a single chronic `value`, no interval columns).

- **draws mode** (`ndraws` set): columns `focal_date`, `site_id` plus
  the bracket envelope summary (`median_*`/`lo_*`/`hi_*` per envelope,
  and `precautionary_lo`/`precautionary_hi` in `"bracket"` mode). The
  chronic frame carries the daily frame it was built from as attribute
  `"daily"`. When `chronic = FALSE`, the
  [`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)
  result is returned unchanged. In all cases the fitted models are
  attached as attributes `"reference_model"` and `"imputation_model"`
  (the latter `NULL` when `impute = FALSE` or the fit was skipped). Note
  that most dplyr verbs drop attributes, so read them off the returned
  frame directly.

## Details

The chronic step "smears" the daily msPAF line with an exponential-decay
memory kernel: because the daily series is equally spaced, forward-step
duration weighting is uniform (Δt = 1 day) and the chronic value is
shaped only by the `tau` decay. Set `chronic = FALSE` to stop at the
daily series and reproduce the previous (daily) return.

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
[`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md),
[`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md),
[`expand_focal_dates()`](https://vorpalvorpal.github.io/hydroSense/reference/expand_focal_dates.md)

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
