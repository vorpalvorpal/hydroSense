# Fit a temporal reference model for contemporaneous ARA background subtraction

Fits per-analyte temporal models on reference-site chemistry so that
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)
can subtract a *contemporaneous* background (what the reference site
would have shown at the same moment as the target sample) rather than a
static time-average.

## Usage

``` r
fit_reference_model(
  reference,
  hydro = NULL,
  hydro_type = "rainfall",
  latitude = NULL,
  longitude = NULL,
  conc_units = NULL,
  analyte_metadata = NULL,
  imputation_model = NULL,
  match_window_days = 5L,
  match_hydro_tol = NULL,
  api_tau_bounds_short = c(1, 30),
  api_tau_bounds_long = c(20, 365),
  auto_select = TRUE,
  min_obs_model = 20L,
  summary = "geom_mean",
  silo_start = NULL,
  silo_end = NULL,
  silo_api_key = NULL,
  eps = 1e-09
)
```

## Arguments

- reference:

  Long-format reference chemistry data frame. Required columns:
  `sample_id`, `site_id`, `datetime`, `analyte`, `value`, `detected`.
  Toxicant concentrations must be in µg/L; supply them via a
  `units.analyte` column or via `conc_units`.

- hydro:

  Daily hydrology data frame (`date`, `value`), or `NULL` to fetch SILO
  rainfall automatically. See *Hydrology input* above.

- hydro_type:

  Character; one of `"rainfall"`, `"stage"`, `"discharge"`. Ignored if
  `hydro` already has a `type` column. Default `"rainfall"`.

- latitude, longitude:

  Decimal-degree coordinates; required only when `hydro = NULL` (SILO
  auto-fetch path).

- conc_units:

  Unit string (e.g. `"mg/L"`, `"ug/L"`) for reference chemistry when no
  `units.analyte` column is present.

- analyte_metadata:

  Analyte metadata, or `NULL` for the bundled CSV.

- imputation_model:

  Optional `imputation_model` from
  [`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md)
  (fit on the reference site's own chemistry). When supplied, missing
  analytes are imputed in raw concentration space *before* the
  per-analyte models are fitted, so a well-sampled analyte lifts a
  sparsely-sampled one into a richer spread of hydrological regimes.
  Imputed rows (`detected = TRUE`) are used as model anchors alongside
  measured observations. Requires **brms**. Default `NULL` (measured
  observations only).

- match_window_days:

  Integer; tier-1 time tolerance in days (default 5).

- match_hydro_tol:

  Numeric; tier-1 hydro tolerance (default `NULL` → `0.5 × IQR` of the
  reference event-API series).

- api_tau_bounds_short:

  Length-2 numeric `c(lo, hi)` search range (days) for the short-store
  recession constant `tau_short` (default `c(1, 30)`). A degenerate
  `c(x, x)` fixes `tau_short = x`.

- api_tau_bounds_long:

  Length-2 numeric `c(lo, hi)` search range (days) for the long-store
  recession constant `tau_long` (default `c(20, 365)`).

- auto_select:

  Logical; if `TRUE` (default), select `tau_short`, `tau_long` per
  analyte by profiled AIC. If `FALSE`, use the parsimonious defaults
  (`tau_short = 7`, `tau_long = 60`) for all analytes.

- min_obs_model:

  Integer; minimum detected observations required to attempt tier-2
  modelling (default `20L`). Analytes below this fall directly to tier
  3.

- summary:

  Static-fallback summary statistic (tier 3): one of `"geom_mean"`
  (default), `"median"`, `"arith_mean"`, `"p80"`, `"p90"`, `"p95"`.

- silo_start, silo_end:

  Start/end dates for SILO auto-fetch. Default `NULL`: derived from the
  reference chemistry date range, padded on the left by
  `5 × max(api_tau_bounds_long)` days so the recursive reservoir has
  enough burn-in (about 5 tau) to converge before the first observation.

- silo_api_key:

  Passed to
  [`get_silo_rainfall()`](https://vorpalvorpal.github.io/hydroSense/reference/get_silo_rainfall.md)
  when `hydro = NULL`.

- eps:

  Small positive guard for log transform (default `1e-9`).

## Value

An object of class `reference_model`:

- `$models`:

  Named list (one per analyte) carrying `gamm_fit`, `tau_short`,
  `tau_long`, `best_aic`, `null_aic`, `tier`, `n_obs`, `static_ref`, and
  `obs` (normalised observations with hydro features — used for tier-1
  matching).

- `$hydro`:

  Daily hydro series used for index computation.

- `$hydro_type`:

  `"rainfall"`, `"stage"`, or `"discharge"`.

- `$match_window_days`:

  Tier-1 time window.

- `$match_hydro_tol`:

  Tier-1 hydro tolerance (computed or supplied).

- `$static_ref`:

  Tibble `(analyte, ref_norm)` — tier-3 fallback for all analytes.

- `$fit_date`:

  Date the model was fitted.

## Details

**Package invariant:** reference and target are assumed to share the
same catchment. The hydrology series supplied here (or fetched
automatically from SILO) applies to both sites. This invariant should
hold for near-field leachate assessments where a headwater reference is
paired with a downstream target in the same sub-catchment.

**Three-tier resolver** (evaluated per analyte × target date):

1.  **Tier 1 — direct match.** A reference observation within
    `±match_window_days` days whose event-API is within
    `match_hydro_tol` (default: `0.5 × IQR` of the reference API
    series). This gate rejects time-close but hydrologically-mismatched
    grabs (e.g. a dry-weather reference next to a wet-event target).

2.  **Tier 2 — GAM prediction.** A per-analyte
    [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) with
    `s(doy, bs="cc") + s(hydro_short) + s(hydro_long)`. The reservoir
    recession constants `tau_short`, `tau_long` are selected per analyte
    by profiled AIC over the
    `api_tau_bounds_short`/`api_tau_bounds_long` ranges (continuous,
    with a `tau_long >= 1.5*tau_short` separation and a ΔAIC ≥ 2
    adoption gate over a parsimonious default). The analyte falls back
    to tier 3 if the best model has higher AIC than the null
    (intercept-only) model, or if fewer than `min_obs_model` detected
    observations are available.

3.  **Tier 3 — static fallback.** Geometric mean of all normalised
    reference observations (identical to
    [`prepare_reference()`](https://vorpalvorpal.github.io/hydroSense/reference/prepare_reference.md)
    with `summary = "geom_mean"`).

**Hydrology input** — supply exactly one of:

- `hydro`: a data frame with columns `date` (Date) and `value`
  (numeric), plus a `type` column OR supply `hydro_type =` separately.
  Supported types: `"rainfall"` (daily mm; → API), `"stage"` (m; →
  antecedent mean), `"discharge"` (m³/s; → antecedent mean).

- `latitude` + `longitude`: if `hydro = NULL`, SILO daily rainfall is
  fetched automatically via
  [`get_silo_rainfall()`](https://vorpalvorpal.github.io/hydroSense/reference/get_silo_rainfall.md).
  Requires the **weatherOz** package and a valid SILO API key.

## See also

[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md),
[`get_silo_rainfall()`](https://vorpalvorpal.github.io/hydroSense/reference/get_silo_rainfall.md),
[`prepare_reference()`](https://vorpalvorpal.github.io/hydroSense/reference/prepare_reference.md),
[`ara_summary()`](https://vorpalvorpal.github.io/hydroSense/reference/ara_summary.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Fit using SILO auto-fetch
ref <- subset(leachate_demo(), site_id == "reference")
ref_model <- fit_reference_model(
  reference  = ref,
  latitude   = -33.87,
  longitude  = 151.21,
  conc_units = "ug/L"
)
ref_model

# Or supply your own gauge record
ref_model2 <- fit_reference_model(
  reference  = ref,
  hydro      = my_stage_df, # data.frame(date, value)
  hydro_type = "stage",
  conc_units = "ug/L"
)
} # }
```
