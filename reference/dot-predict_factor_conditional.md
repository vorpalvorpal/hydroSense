# Closed-form conditional prediction for the Route C factor model

For each eligible sample, builds the residual vector `y` (observed minus
the Stage-1 GAM mean for detected target cells; `NA` for BDL/missing
cells), then applies
[`.factor_condition()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-factor_condition.md)
/
[`.factor_condition_draw()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-factor_condition_draw.md)
per posterior draw of `(Lambda, Psi)` to predict the missing cells –
conditioned on this sample's *own* observed analytes (finding 3), with
any BDL target truncated at `log(DL)` (findings 1-2). Mirrors
[`.predict_factor_long()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-predict_factor_long.md)'s
`pm_long` return shape so
[`.predict_and_merge()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-predict_and_merge.md)'s
merge/fabricate logic is unchanged; only rows for cells that actually
need prediction (BDL or entirely absent) are emitted.

## Usage

``` r
.predict_factor_conditional(
  group,
  pc_wide,
  df_eligible,
  return,
  ndraws,
  batch_size
)
```

## Arguments

- group:

  A fitted `factor`-method group.

- pc_wide:

  Per-sample PC scores for the eligible samples (`sample_id`

  - `PC*` columns).

- df_eligible:

  The long-format input rows for the eligible samples (supplies each
  sample's own observed target values).
