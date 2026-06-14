# Variance-stabilising transform of the ARA impact

Maps the signed impact `impact` (`I` = `C_norm - ref_norm`, normalised
concentration units) onto the `asinh(I / c)` scale, on which the daily
residual smoother's process variance is proportional to concentration
rather than constant. See
[`.g_inverse()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-g_inverse.md)
for the inverse and
[`.analyte_c()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-analyte_c.md)
for the per-analyte scale `c`.

## Usage

``` r
.g_transform(impact, scale_c)
```

## Arguments

- impact:

  Numeric vector of impacts `C_norm - ref_norm` (may be negative, `NA`,
  or `Inf`).

- scale_c:

  Single positive scale `c` (the additive-\>proportional crossover; the
  per-analyte SSD HC5, see
  [`.analyte_c()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-analyte_c.md)).

## Value

Numeric vector `asinh(impact / scale_c)`, same length as `impact`.
`NA`/`Inf` propagate; `+/-Inf` map to `+/-Inf`.
