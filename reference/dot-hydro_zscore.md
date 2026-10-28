# Standardised short-window hydro feature over query dates (q-modulation input)

Drives the hydrology-modulated process variance of the residual smoother
(`q_mult = exp(kappa * z)`). `api_cache` (from
[`.api_cache_new()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-api_cache_new.md)
over `target_model$hydro`) shares the API across analytes at the same
tau.

## Usage

``` r
.hydro_zscore(target_model, qdates, m, api_cache = NULL)
```
