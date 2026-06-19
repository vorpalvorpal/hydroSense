# Impute missing and BDL chemistry using a fitted imputation model

Applies the models fitted by
[`fit_imputation_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_imputation_model.md)
to `df`, returning posterior estimates for below-detection-limit (BDL)
and missing observations in every fitted group.

## Usage

``` r
impute_chemistry(
  df,
  model,
  apply_hurdles = TRUE,
  bdl_cap = TRUE,
  return = c("point", "draws"),
  ndraws = NULL,
  batch_size = NULL
)
```

## Arguments

- df:

  Long-format chemistry data frame (same schema as used for fitting).

- model:

  Fitted model from
  [`fit_imputation_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_imputation_model.md).

- apply_hurdles:

  Logical. Apply each group's presence hurdle? Default `TRUE`. When
  `FALSE`, every sample is eligible for every group.

- bdl_cap:

  Logical. Cap imputed BDL values at the original detection limit?
  Default `TRUE`. Applied to all methods: an imputed below-detection
  value should never exceed its detection limit.

- return:

  `"point"` (default) for posterior mean per cell; `"draws"` for one row
  per (sample × analyte × draw).

- ndraws:

  Integer or `NULL`. Use only this many posterior draws for prediction
  (subsampled). `NULL` (default) uses all draws. Lowering it reduces
  memory/time, at some cost to interval precision.

- batch_size:

  Integer or `NULL`. Predict eligible samples in batches of this many
  rows to bound peak memory (important for `"rescor_mi"`, whose `mi()`
  prediction is memory-heavy). `NULL` (default) predicts all at once.

## Value

`df` with BDL and missing cells in every fitted group replaced by
posterior mean estimates, plus columns:

- `imputed` (logical) — `TRUE` for filled cells

- `imputed_kind` — `"observed"`, `"censored_left"`, or `"missing"`

When `bdl_cap = TRUE` and any imputed BDL cell exceeded its detection
limit, a per-cell audit of the cap activations is attached as the
`"bdl_cap_summary"` attribute; retrieve it with
[`bdl_cap_summary()`](https://vorpalvorpal.github.io/leachatetools/reference/bdl_cap_summary.md).

## Details

**Completing the panel**

For each eligible sample (see *Hurdles*), every target analyte the group
models is filled: BDL cells are replaced with their posterior estimate,
and analytes that are **entirely absent** for a sample gain a new
model-anchored row (`imputed_kind = "missing"`, `detected = TRUE`). This
is what lets a well-sampled analyte lift a sparsely-sampled one — e.g. a
sample with Zn but no Cu gains a Cu row predicted from the fitted Cu–Zn
relationship. Fabricated rows carry the originating sample's `site_id`,
`datetime`, and any other sample-level columns.

Fabricated rows are model predictions, not measurements. With
`return = "draws"` the imputation-model uncertainty propagates through
the draw carrier; with `return = "point"` the single anchor enters any
downstream model as if observed, which can understate that model's
uncertainty. Prefer `return = "draws"` when the imputed values feed a
further model (e.g. the reference/target GAMs).

**Hurdles**

Each fitted group may carry a presence hurdle (see
[`impute_group()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_group.md)).
When `apply_hurdles = TRUE`, imputed values for a group are only
returned for samples carrying at least one of that group's hurdle
analytes (detected or BDL) — e.g. under the leachate preset, a sample
with no metals recorded is not given imputed metals, because a leachate
metal pulse may simply not have arrived at that location yet. Samples
failing a hurdle pass through unchanged (non-imputed values preserved;
BDL values remain flagged as BDL).

## See also

[`fit_imputation_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_imputation_model.md),
[`bdl_cap_summary()`](https://vorpalvorpal.github.io/leachatetools/reference/bdl_cap_summary.md)

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_imputation_model(monitoring_long)
# Point estimates (default), or full posterior draws with return = "draws":
imputed <- impute_chemistry(monitoring_long, model, return = "point")
} # }
```
