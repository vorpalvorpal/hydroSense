# Antecedent Precipitation Index (API) — exact recursive linear reservoir

Implements the exact discrete solution of the single linear store
`dS/dt = -S/tau + P(t)`:

## Usage

``` r
.compute_api(hydro_values, hydro_dates, target_dates, tau)
```

## Arguments

- hydro_values:

  Numeric vector of daily rainfall values (`NA` is 0 input).

- hydro_dates:

  Date vector matching `hydro_values`.

- target_dates:

  Date vector; API is evaluated at each.

- tau:

  Positive numeric; reservoir memory constant (days). The weight for a
  gap of `dt` days is `exp(-dt / tau)`.

## Value

Numeric vector, same length as `target_dates`.

## Details

    S_t = k^{dt} * S_{t-1} + P_t,   k = exp(-1 / tau)

This is the convergent infinite-horizon form of the classical Antecedent
Precipitation Index (Kohler & Linsley 1951), whose theoretical basis is
the Maillet (1905) exponential recession. Unlike a windowed sum it has
no truncation horizon: all prior rainfall contributes with exponentially
decaying weight `k^{dt}` set by the actual day gap `dt`, so irregular
spacing is handled exactly.

Algorithm: lay the rainfall onto a complete daily grid (`NA` and gap
days contribute 0 input, so the store simply decays across them). On a
daily grid every step has `dt = 1`, so the reservoir reduces to the
first-order linear recursion `S_t = k * S_{t-1} + P_t`, evaluated in C
by [`stats::filter()`](https://rdrr.io/r/stats/filter.html)
(`method = "recursive"`). This is an exact reorganisation of the
per-event form `S_t = k^{dt} * S_{t-1} + P_t` (`k^{dt}` is just `dt`
unit steps with zero input between events), with none of the per-target
re-summation the old windowed form required. Target dates index straight
into the grid; any falling before the first hydro day return 0 (empty
reservoir).

## References

Maillet, E. (1905) *Essais d'hydraulique souterraine et fluviale.*
Hermann, Paris. Kohler, M.A. & Linsley, R.K. (1951) *Predicting the
runoff from storm rainfall.* US Weather Bureau Research Paper 34.
