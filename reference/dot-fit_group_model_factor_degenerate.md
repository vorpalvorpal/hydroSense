# Stage-1-only fit for a single-analyte factor group (no cross-analyte coupling is possible with J = 1)

Fits the Stage-1 GAM as usual and estimates a plug-in residual variance
from the detected residuals. No Stan model is fit; `Lambda` is a J x 0
matrix (`k = 0`) so
[`.factor_condition()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-factor_condition.md)'s
existing `length(obs) == 0` branch (always true here – a single-analyte
group's own target is either fully observed, in which case there is
nothing to predict, or missing, in which case there is nothing else in
the group to condition on) reduces exactly to the marginal
`N(mu_j, Psi_j)`, letting
[`.predict_factor_conditional()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-predict_factor_conditional.md)
and
[`.route_c_draw_params()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-route_c_draw_params.md)
handle this group with no special-casing.

## Usage

``` r
.fit_group_model_factor_degenerate(
  target_analytes,
  safe_analytes,
  base,
  pc_wide,
  pc_cols,
  log_floors,
  group_name = "group",
  post_ndraws = 1000L
)
```
