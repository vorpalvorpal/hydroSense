# Continuous daily AmsPAF time series from interpolated grab chemistry

Interpolates per-analyte grab chemistry onto a daily date grid and
computes AmsPAF for every day within the requested date range. The
result is a tidy tibble with one row per (site \\\times\\ day), suitable
for trend analysis, visualisation, and downstream
[`time_weighted_aggregate()`](https://vorpalvorpal.github.io/leachatetools/reference/time_weighted_aggregate.md)
calls.

## Usage

``` r
amspaf_daily(
  df,
  temperature = NULL,
  reference = NULL,
  reference_model = NULL,
  imputation_model = NULL,
  start = NULL,
  end = NULL,
  by = "day",
  interpolation = c("forward_fill", "linear", "model"),
  leading_edge = c("drop", "backfill"),
  analyte_metadata = NULL,
  method = c("multi", "anzecc"),
  guideline_dir = getOption("leachatetools.guideline_dir"),
  min_analytes = 3L,
  conc_units = NULL,
  require_temperature = TRUE
)
```

## Arguments

- df:

  Long-format grab chemistry data frame. Required columns: `sample_id`,
  `site_id`, `datetime` (Date or POSIXct), `analyte`, `value`,
  `detected`. Optional but propagated: `units.analyte`, `imputed`.
  Chemistry for multiple sites may be stacked; interpolation and AmsPAF
  are computed per site.

- temperature:

  Optional daily water temperature data frame for days without a grab
  temperature measurement. Required columns: `datetime` (Date or
  POSIXct) and `value` (temperature in \\{}^\circ\\C). The output of
  [`estimate_water_temp()`](https://vorpalvorpal.github.io/leachatetools/reference/estimate_water_temp.md)
  is accepted directly (extra columns are ignored). When both this
  argument and grab-sample temperature rows in `df` are present for the
  same day, the grab measurement takes priority. `NULL` (default) means
  temperature must come from `df` rows alone.

- reference:

  Background reference chemistry for ARA adjustment. Accepts the same
  four forms as
  [`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md).
  Controls **only** whether background is subtracted; it is independent
  of `interpolation`. With `interpolation = "model"`, pass the same
  `reference_model` here to assess the leachate-attributable increment,
  or `NULL` to assess total concentration.

- reference_model:

  A `reference_model` from
  [`fit_reference_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_reference_model.md).
  **Required** when `interpolation = "model"` — it supplies the
  background and catchment hydrology used by the season-blind target
  impact model
  ([`fit_target_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_target_model.md))
  that interpolates toxicants between grabs. Ignored for the
  `"forward_fill"` and `"linear"` paths.

- imputation_model:

  Optional `imputation_model` from
  [`fit_imputation_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_imputation_model.md),
  passed to
  [`fit_target_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_target_model.md)
  for tier-2 enrichment under `interpolation = "model"`. Requires
  **brms**.

- start, end:

  Date boundaries for the daily grid. Default: earliest and latest
  `datetime` values in `df`. Coerced to `Date`.

- by:

  Temporal resolution string passed to
  [`seq.Date()`](https://rdrr.io/r/base/seq.Date.html). Default `"day"`.

- interpolation:

  How to fill gaps between grab samples. One of `"forward_fill"`
  (default) or `"linear"`. See the *Interpolation* section.

- leading_edge:

  What to do with days before the first grab sample for an analyte.
  `"drop"` (default) or `"backfill"`. See the *Leading edge* section.

- analyte_metadata:

  Data frame of analyte metadata, or `NULL` to use the bundled metadata.
  Passed to
  [`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md).

- method:

  SSD method. One of `"multi"` (default, model-averaged) or `"anzecc"`.
  Passed to
  [`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md).

- guideline_dir:

  Path to the ANZG guideline data folder. Falls back to
  `getOption("leachatetools.guideline_dir")`.

- min_analytes:

  Minimum number of SSD-eligible analytes per day for AmsPAF to be
  computed. Default `3L`.

- conc_units:

  Character. Concentration units for SSD-eligible rows when `df` lacks a
  `units.analyte` column. Passed to
  [`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md).

- require_temperature:

  Logical (default `TRUE`). When `TRUE`, any daily sample with `NH3-N`
  must also carry a `temperature` value. Passed to
  [`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md).
  Set `FALSE` for datasets without ammonia.

## Value

A tibble with one row per (site \\\times\\ day) for days with sufficient
analyte coverage:

- `date`:

  Date of this daily estimate.

- `site_id`:

  Site identifier.

- `amspaf`:

  Daily AmsPAF (percentage, 0–100+).

- `n_analytes_used`:

  SSD-eligible analytes contributing to AmsPAF.

- `dominant_analyte`:

  Analyte with the highest individual PAF.

- `max_paf`:

  PAF of the dominant analyte (proportion 0–1).

- `n_measured_analytes`:

  SSD-eligible analytes with a direct grab sample on this day (not
  interpolated).

- `days_since_last_sample`:

  Days since the most recent grab sample for any SSD-eligible analyte.
  Helps identify heavily interpolated regions.

- `analyte_pafs`:

  List column of per-analyte diagnostic tibbles (same structure as from
  [`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md)).

An `"ara_summary"` attribute is attached; retrieve it with
[`ara_summary()`](https://vorpalvorpal.github.io/leachatetools/reference/ara_summary.md).

## Interpolation

Each analyte is interpolated independently:

- `"forward_fill"` carries the last directly observed value forward
  until the next grab. Produces a step function.

- `"linear"` linearly interpolates between consecutive observations. For
  SSD-eligible toxicants (metals, ammonia), interpolation is performed
  in log-concentration space, which is more appropriate for log-normally
  distributed data and avoids negative intermediate values. Co-analytes
  (pH, temperature, DOC, hardness) are interpolated linearly.

- `"model"` fits a season-blind
  [`fit_target_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_target_model.md)
  on the grab chemistry and the supplied `reference_model`, and predicts
  each toxicant's concentration between grabs as `reference + impact`,
  where the impact (the leachate-attributable increment) is modelled
  from hydrology and a persistent latent state but **never** from
  day-of-year. Co-analytes are forward-filled. Requires
  `reference_model`; see
  [`fit_target_model()`](https://vorpalvorpal.github.io/leachatetools/reference/fit_target_model.md)
  for the method.

Below-detection values are treated as their detection-limit value for
interpolation purposes, matching the treatment in
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md).

## Leading edge

Days before the first grab sample for an analyte are outside the
observation record:

- `"drop"` (default) excludes such days from the output.

- `"backfill"` assigns the first observed value to all earlier days. Use
  cautiously – it assumes the analyte was at its first observed level
  before sampling commenced.

## Temperature for ammonia

When `NH3-N` is in `df`, water temperature is required per sample for
the un-ionised fraction normalisation. Two sources are accepted:

- Temperature rows already in `df` (measured on grab-sample days) are
  interpolated along with other analytes.

- An external daily temperature series supplied via the `temperature`
  argument (e.g. from
  [`estimate_water_temp()`](https://vorpalvorpal.github.io/leachatetools/reference/estimate_water_temp.md))
  fills in days without a grab temperature. Grab-sample-day measurements
  take priority.

Set `require_temperature = FALSE` only for datasets that do not contain
ammonia.

## See also

[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md),
[`time_weighted_aggregate()`](https://vorpalvorpal.github.io/leachatetools/reference/time_weighted_aggregate.md),
[`estimate_water_temp()`](https://vorpalvorpal.github.io/leachatetools/reference/estimate_water_temp.md),
[`get_silo_air_temp()`](https://vorpalvorpal.github.io/leachatetools/reference/get_silo_air_temp.md)

## Examples

``` r
# \donttest{
demo <- leachate_demo()
ds  <- subset(demo, site_id == "downstream")
out <- amspaf_daily(ds, require_temperature = FALSE)
head(out[, c("date", "site_id", "amspaf", "n_measured_analytes",
             "days_since_last_sample")])
#> # A tibble: 6 × 5
#>   date       site_id    amspaf n_measured_analytes days_since_last_sample
#>   <date>     <chr>       <dbl>               <int>                  <int>
#> 1 2024-01-15 downstream   93.3                   3                      0
#> 2 2024-01-16 downstream   93.3                   0                      1
#> 3 2024-01-17 downstream   93.3                   0                      2
#> 4 2024-01-18 downstream   93.3                   0                      3
#> 5 2024-01-19 downstream   93.3                   0                      4
#> 6 2024-01-20 downstream   93.3                   0                      5
# }
```
