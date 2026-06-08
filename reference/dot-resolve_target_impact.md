# Predict the site impact (and implied normalised concentration) at query dates

For each query (date × analyte): `I_hat = beta·f(hydro) + S_interp`,
with `beta·f(hydro)` from the fitted response (0 for bridge-tier
analytes) and `S_interp` from
[`.interp_residual()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-interp_residual.md).
Also returns `ref_norm` (from the reference model) and the implied
`C_norm = max(ref_norm + I_hat, 0)`.

## Usage

``` r
.resolve_target_impact(target_model, query, analytes = NULL, wq = NULL)
```

## Arguments

- target_model:

  A `target_model`.

- query:

  Tibble with `date` (Date) — the days to predict.

- analytes:

  Character; analytes to predict (default: all modelled).

- wq:

  Optional long-format water-quality data frame (`sample_id`, `analyte`,
  `value`) giving each query day's WQ, with `sample_id` equal to the
  query date as a character string. When supplied and the analyte has a
  fitted WQ layer, the day's concentration is predicted as
  `WQ-prediction + d_interp` (tier `"wq"`) instead of `ref + impact`.

## Value

Tibble `(date, analyte, ref_norm, impact, C_norm, impact_tier)`.
