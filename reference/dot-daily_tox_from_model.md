# Model-based daily toxicant interpolation (season-blind impact model)

Thin wrapper: calls
[`.fit_daily_target()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-fit_daily_target.md)
once then
[`.predict_daily_tox()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-predict_daily_tox.md)
once (point mode, no ε). Modelled-toxicant rows in `daily_long` are
replaced; co-analytes and non-modelled toxicants are left untouched. On
any failure the input is returned unchanged.

## Usage

``` r
.daily_tox_from_model(
  daily_long,
  site_rows,
  reference_model,
  imputation_model,
  conc_units,
  meta,
  tox_analytes,
  method = "multi",
  guideline_dir = NULL
)
```

## Arguments

- daily_long:

  Output of
  [`.build_daily_chem()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-build_daily_chem.md)
  (+ temperature fill).

- site_rows:

  Grab chemistry for the site (passed to the target model).

- reference_model:

  A `reference_model`.

- imputation_model:

  Optional `imputation_model` (tier-2 enrichment).

- conc_units:

  Units for `site_rows` toxicants when no `units.analyte`.

- meta:

  Analyte metadata (normalisation formulas).

- tox_analytes:

  SSD-eligible analyte names.

## Value

`daily_long` with modelled-toxicant rows replaced.

## Details

For draw-mode orchestration (Chunk E), call the two underlying helpers
directly: fit once with
[`.fit_daily_target()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-fit_daily_target.md),
then loop N times over
`.predict_daily_tox(fdm, perturbed_tm, eps_paths)`.
