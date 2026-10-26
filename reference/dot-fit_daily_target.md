# Fit-once scaffold for the season-blind daily target model (issue \#16)

Calls
[`fit_target_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_target_model.md)
once, builds the static scaffolding needed for per-draw prediction, and
precomputes the OU bridge Cholesky factors for every analyte over the
full daily date grid. The result ("fitted daily model", `fdm`) is the
input to
[`.predict_daily_tox()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-predict_daily_tox.md);
call `.predict_daily_tox(fdm)` once for point mode or N times (with a
perturbed target model and OU-sampled `eps_paths`) for draw mode.

## Usage

``` r
.fit_daily_target(
  site_rows,
  reference_model,
  imputation_model,
  conc_units,
  meta,
  tox_analytes,
  daily_long,
  ou_scale = 1,
  grab_cv = NULL,
  kappa = 0.5,
  method = "multi",
  guideline_dir = NULL,
  transform = "pseudo_log",
  ref_q = NULL
)
```

## Arguments

- site_rows:

  Grab chemistry for the site.

- reference_model, imputation_model, conc_units, meta, tox_analytes:

  Forwarded to
  [`fit_target_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_target_model.md).

- daily_long:

  Forward-filled daily chemistry for the site; supplies the date grid
  (`qdates`) and co-analyte values.

- ref_q:

  Optional pre-resolved reference norm
  `(sample_id, analyte, ref_norm, ref_tier, date)` covering every date
  of `qdates` (as resolved once per call by
  [`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)).
  `NULL` resolves it here.

## Value

A named list (`fdm`) or `NULL` on model failure.
