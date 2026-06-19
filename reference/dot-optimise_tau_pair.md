# Minimise an AIC objective over a (tau_short, tau_long) pair

Deterministic golden-section coordinate descent
([`stats::optimize()`](https://rdrr.io/r/stats/optimize.html)) over the
two reservoir recession constants, holding one fixed while the other is
optimised (two passes). The fast/slow separation
`tau_long >= 1.5*tau_short` is imposed by capping tau_short's upper
bound at `tau_long/1.5` and lifting tau_long's lower bound to
`1.5*tau_short`. A degenerate bound `lo == hi` pins that store. Shared
by the reference
([`.select_api_tau()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-select_api_tau.md))
and impact
([`.fit_impact_response()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-fit_impact_response.md))
hydrology fits. No RNG.

## Usage

``` r
.optimise_tau_pair(
  aic_fn,
  tau_bounds_short,
  tau_bounds_long,
  default_short = .REF_TAU_DEFAULT_SHORT,
  default_long = .REF_TAU_DEFAULT_LONG
)
```

## Arguments

- aic_fn:

  Function `(tau_short, tau_long) -> AIC` (`Inf` on fit failure).

- tau_bounds_short, tau_bounds_long:

  Length-2 numeric `c(lo, hi)` (days).

- default_short, default_long:

  Starting tau (days).

## Value

List `(tau_short, tau_long)`.
