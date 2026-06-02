# Un-ionised ammonia fraction (NH3) of total ammonia-N

Un-ionised ammonia fraction (NH3) of total ammonia-N

## Usage

``` r
.nh3_unionised_fraction(pH, temperature_C)
```

## Arguments

- pH:

  Sample pH.

- temperature_C:

  Sample temperature (°C).

## Value

Fraction (0–1) of total ammonia-N present as un-ionised NH3, using
`f = 1 / (1 + 10^(pKa - pH))` with `pKa = 0.09018 + 2729.92 / T(K)`
(Emerson et al. 1975), `T(K) = temperature_C + 273.15`.
