# Standardised short-window hydro feature over query dates (q-modulation input)

Drives the hydrology-modulated process variance of the residual smoother
(`q_mult = exp(kappa * z)`).

## Usage

``` r
.hydro_zscore(target_model, qdates, m)
```
