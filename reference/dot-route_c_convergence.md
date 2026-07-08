# Rotation-invariant convergence diagnostics for a Route C factor fit

A low-rank factor model is only identified up to rotation of `Lambda`
(`Lambda R` gives the same `Lambda Lambda'`). With the positive-lower-
triangular constraint this is mostly pinned, but at small `N` / weak
signal the sampler can still explore near-rotationally-equivalent modes,
inflating the per-element Rhat of `Lambda` **without any consequence for
imputation** – prediction uses only `Sigma = Lambda Lambda' + Psi`,
which is rotation- invariant (see
[`.factor_condition()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-factor_condition.md)).
So the meaningful convergence check monitors the invariant functionals:
the implied residual covariance `Sigma`, the idiosyncratic variances
`psi`, and the log-density `lp__` – not raw `Lambda`. Empirically
`Sigma` converges cleanly (Rhat ~ 1.0) even when a `Lambda` element
trips a naive `max(rhat(fit))` gate.

## Usage

``` r
.route_c_convergence(group)
```

## Arguments

- group:

  A fitted `factor`-method group (from
  [`.fit_group_model_factor()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-fit_group_model_factor.md)).

## Value

A list with `sigma_rhat`, `psi_rhat`, `lp_rhat` (max Rhat over the
`Sigma` entries, over `psi`, and for `lp__`); all `NA` for a degenerate
(single-analyte, Stan-free) group.
