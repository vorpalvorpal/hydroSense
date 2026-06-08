# Model-based daily toxicant interpolation (season-blind impact model)

Fits a
[`fit_target_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_target_model.md)
on the site's grab chemistry and the supplied `reference_model`,
predicts each toxicant's normalised concentration
`C_norm = ref_norm + impact` on every daily date, and reconstructs the
raw µg/L concentration `C_raw = C_norm / factor` (the normalisation is
multiplicative, so the factor is
`normalise(1, co-analytes_of_the_day)`). Modelled-toxicant rows in
`daily_long` are replaced with these estimates (in µg/L,
`units.analyte = "ug/L"`); co-analytes and non-modelled toxicants are
left untouched. On any failure the input is returned unchanged.

## Usage

``` r
.daily_tox_from_model(
  daily_long,
  site_rows,
  reference_model,
  imputation_model,
  conc_units,
  meta,
  tox_analytes
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
