# Impute missing and BDL chemistry using a fitted imputation model

Applies the models fitted by
[`fit_imputation_model()`](https://www.kedumba.com.au/leachatetools/reference/fit_imputation_model.md)
to `df`, returning posterior mean estimates for missing and
below-detection-limit (BDL) observations in the metals and organics
groups.

## Usage

``` r
impute_chemistry(
  df,
  model,
  metal_hurdle = TRUE,
  organic_hurdle = TRUE,
  bdl_cap = TRUE,
  return = c("point", "draws")
)
```

## Arguments

- df:

  Long-format chemistry data frame (same schema as used for fitting).

- model:

  Fitted model from
  [`fit_imputation_model()`](https://www.kedumba.com.au/leachatetools/reference/fit_imputation_model.md).

- metal_hurdle:

  Logical. Apply metal-presence hurdle? Default `TRUE`.

- organic_hurdle:

  Logical. Apply DOC-like-presence hurdle? Default `TRUE`.

- bdl_cap:

  Logical. Cap imputed BDL values at the original detection limit?
  Default `TRUE`.

- return:

  `"point"` (default) for posterior mean per cell; `"draws"` for one row
  per (sample × analyte × draw).

## Value

`df` with BDL and missing cells in the metals/organics groups replaced
by posterior mean estimates, plus columns:

- `imputed` (logical) — `TRUE` for filled cells

- `imputed_kind` — `"observed"`, `"censored_left"`, or `"missing"`

## Details

**Hurdles**

Imputed values are only returned for samples that meet the relevant
hurdle:

- *Metals*: at least one metal analyte present (detected or BDL) at the
  sample. Samples with no metals recorded are not imputed — a leachate
  metal pulse may simply not have arrived at this location yet.

- *Organics*: at least one of DOC, TOC, BOD, COD, cBOD present at the
  sample.

Samples failing a hurdle pass through unchanged (non-imputed values are
preserved; BDL values remain flagged as BDL).

## See also

[`fit_imputation_model()`](https://www.kedumba.com.au/leachatetools/reference/fit_imputation_model.md)

## Examples

``` r
if (FALSE) { # \dontrun{
model <- fit_imputation_model(monitoring_long)
# Point estimates (default), or full posterior draws with return = "draws":
imputed <- impute_chemistry(monitoring_long, model, return = "point")
} # }
```
