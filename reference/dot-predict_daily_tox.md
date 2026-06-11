# Per-draw daily toxicant prediction from a pre-fitted scaffold (issue \#16)

Predicts normalised and raw concentrations for one draw (or for the
deterministic centre line when `residual_paths = NULL`). Always uses the
precomputed static scaffolding from `fdm`; accepts an
optionally-perturbed target model (`tm_p`) and a per-analyte residual
path (`residual_paths`) for draw mode. Pass overrides to `co_split` and
`wq_long` (from
[`.perturb_co_split()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-perturb_co_split.md))
for S7 co-analyte measurement-error draws.

## Usage

``` r
.predict_daily_tox(
  fdm,
  tm_p = fdm$tm,
  residual_paths = NULL,
  co_split = fdm$co_split,
  wq_long = fdm$wq_long
)
```

## Arguments

- fdm:

  Output of
  [`.fit_daily_target()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-fit_daily_target.md).

- tm_p:

  Target model to predict with (default: `fdm$tm`). For draw mode, pass
  a GAM-perturbed copy from
  [`.perturb_target_model()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-perturb_target_model.md).

- residual_paths:

  Named list of residual paths (impact `S` or WQ `d`), one numeric
  vector of length `length(fdm$qdates)` per analyte (`NA` outside the
  analyte's clipped grab span). `NULL` (the default) uses the smoother
  posterior means — the deterministic centre line.

- co_split:

  Per-date co-analyte lookup (default: `fdm$co_split`). Chunk D supplies
  a perturbed version for S7 draws.

- wq_long:

  WQ layer data (default: `fdm$wq_long`).

## Value

Tibble of model rows (`.date`, `value`, `detected`, `.measured`,
`analyte`, `units.analyte`) with `attr("impact_tiers")` attached, or
`NULL` if prediction produced no finite rows.
