# Select reservoir recession constants (tau) per analyte by profiled AIC

Chooses the short- and long-store recession constants `tau_short`,
`tau_long` (days) for the rainfall reservoir by minimising the fitted
GAM's AIC. Given tau the GAM is conditionally linear, so tau is selected
by an *outer* 1-D optimisation of the profiled AIC — profile likelihood
for a nonlinear-in-parameter feature (Wood 2017). The two stores are
optimised by deterministic golden-section search
([`stats::optimize()`](https://rdrr.io/r/stats/optimize.html)) in a
short coordinate descent (tau_short with tau_long held, then tau_long
with tau_short held); no RNG is used.

## Usage

``` r
.select_api_tau(
  obs,
  hydro,
  hydro_type,
  tau_bounds_short,
  tau_bounds_long,
  default_short = 7,
  default_long = 60,
  eps = 1e-09
)
```

## Arguments

- obs:

  Tibble `(date, value_norm)` for one analyte.

- hydro:

  Daily hydro series.

- hydro_type:

  Hydro type string.

- tau_bounds_short, tau_bounds_long:

  Length-2 numeric `c(lo, hi)` search ranges (days). A degenerate
  `lo == hi` fixes that store's tau.

- default_short, default_long:

  Parsimonious fallback tau (days).

- eps:

  Log guard.

## Value

List `(tau_short, tau_long, best_aic, null_aic)`.

## Details

Guards:

- **Separation** `tau_long >= 1.5 * tau_short` keeps the fast and slow
  stores distinct (otherwise the two smooths become collinear).

- **Overfitting gate** — the optimised tau is adopted only when it
  improves AIC by `>= 2` over a parsimonious default
  `(default_short, default_long)`; otherwise the default is returned.
  ΔAIC ≥ 2 is the conventional threshold for a distinguishable model
  (Burnham & Anderson 2002), so a flat or uninformative AIC surface
  falls back to the simpler model rather than chasing noise.

## References

Wood, S.N. (2017) *Generalized Additive Models: An Introduction with R*,
2nd ed. CRC Press. Burnham, K.P. & Anderson, D.R. (2002) *Model
Selection and Multimodel Inference*, 2nd ed. Springer.
