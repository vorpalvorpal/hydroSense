# Declare an imputation group

The imputation engine
([`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md))
is domain-agnostic: it imputes one or more **groups** of target analytes
from a shared PCA-compressed chemistry context, with cross-target
residual correlation within each group. `impute_group()` describes a
single group — which analytes it models and which (if any) presence
hurdle gates it. Pass a list of these as the `groups` argument of
[`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md).

## Usage

``` r
impute_group(name, targets = NULL, hurdle = NULL)
```

## Arguments

- name:

  Group label (a non-empty string). Used as the name of the group's slot
  in the fitted model's `$groups` list and in console messages.

- targets:

  Character vector of analyte names to model jointly in this group, or
  `NULL` to mark this as the **catch-all** group, which claims every
  remaining modellable analyte (those not excluded, not used as PCA
  predictors, and not already claimed by an earlier group). At most one
  catch-all group is allowed in a `groups` list.

- hurdle:

  Character vector of analyte names defining a *presence hurdle*, or
  `NULL` for no hurdle. When a hurdle is set, a sample is only imputed
  for this group if it carries at least one of these analytes (detected
  or below-detection). This stops silence being mistaken for absence —
  e.g. a sample with no metals recorded is left alone rather than given
  invented metal values.

## Value

An object of class `"impute_group"`: a list with elements `name`,
`targets`, and `hurdle`.

## Details

Each group is fitted as its own joint `brms` model, so analytes in
different groups do not borrow residual correlation from one another.
Group together analytes you expect to co-vary (e.g. metals that move
together in a plume).

## See also

[`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md),
[`leachate_impute_groups()`](https://vorpalvorpal.github.io/hydroSense/reference/leachate_impute_groups.md)

## Examples

``` r
# A metals group hurdled on metal presence, plus an everything-else group
# hurdled on dissolved-organic-carbon presence:
groups <- list(
  impute_group("metals", targets = c("Cu", "Zn", "Ni", "Pb"),
               hurdle = c("Cu", "Zn", "Ni", "Pb")),
  impute_group("organics", targets = NULL, hurdle = c("DOC", "TOC"))
)
```
