# Continuous daily msPAF time series from interpolated grab chemistry

Interpolates per-analyte grab chemistry onto a daily date grid and
computes msPAF for every day within the requested date range. The result
is a tidy tibble with one row per (site \\\times\\ day), suitable for
trend analysis, visualisation, and downstream
[`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md)
calls.

## Usage

``` r
mspaf_daily(
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
  guideline_dir = getOption("hydroSense.guideline_dir"),
  min_analytes = 3L,
  conc_units = NULL,
  require_temperature = TRUE,
  ndraws = NULL,
  seed = NULL,
  return = c("summary", "draws"),
  interval = 0.9,
  central = c("median", "mean"),
  grab_cv = NULL,
  ou_scale = 1,
  kappa = 0.5,
  parallel = FALSE,
  couple_residuals = TRUE,
  gap_uncertainty = c("bracket", "ignorable", "informative"),
  transform = c("pseudo_log", "additive")
)
```

## Arguments

- df:

  Long-format grab chemistry data frame. Required columns: `sample_id`,
  `site_id`, `datetime` (Date or POSIXct), `analyte`, `value`,
  `detected`. Optional but propagated: `units.analyte`, `imputed`.
  Chemistry for multiple sites may be stacked; interpolation and msPAF
  are computed per site.

- temperature:

  Optional daily water temperature data frame for days without a grab
  temperature measurement. Required columns: `datetime` (Date or
  POSIXct) and `value` (temperature in \\{}^\circ\\C). The output of
  [`estimate_water_temp()`](https://vorpalvorpal.github.io/hydroSense/reference/estimate_water_temp.md)
  is accepted directly (extra columns are ignored). When both this
  argument and grab-sample temperature rows in `df` are present for the
  same day, the grab measurement takes priority. `NULL` (default) means
  temperature must come from `df` rows alone.

- reference:

  Background reference chemistry for ARA adjustment. Accepts the same
  four forms as
  [`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md).
  Controls **only** whether background is subtracted; it is independent
  of `interpolation`. With `interpolation = "model"`, pass the same
  `reference_model` here to assess the leachate-attributable increment,
  or `NULL` to assess total concentration.

- reference_model:

  A `reference_model` from
  [`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md).
  **Required** when `interpolation = "model"` — it supplies the
  background and catchment hydrology used by the season-blind target
  impact model
  ([`fit_target_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_target_model.md))
  that interpolates toxicants between grabs. Ignored for the
  `"forward_fill"` and `"linear"` paths.

- imputation_model:

  Optional `imputation_model` from
  [`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md),
  passed to
  [`fit_target_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_target_model.md)
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
  [`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md).

- method:

  SSD method. One of `"multi"` (default, model-averaged) or `"anzecc"`.
  Passed to
  [`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md).

- guideline_dir:

  Path to the ANZG guideline data folder. Falls back to
  `getOption("hydroSense.guideline_dir")`.

- min_analytes:

  Minimum number of SSD-eligible analytes per day for msPAF to be
  computed. Default `3L`.

- conc_units:

  Character. Concentration units for SSD-eligible rows when `df` lacks a
  `units.analyte` column. Passed to
  [`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md).

- require_temperature:

  Logical (default `TRUE`). When `TRUE`, any daily sample with `NH3-N`
  must also carry a `temperature` value. Passed to
  [`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md).
  Set `FALSE` for datasets without ammonia.

- ndraws:

  Positive integer or `NULL` (default). `NULL` returns the
  **deterministic** daily msPAF: the fast, grabs-exact point estimate
  (the default and recommended best guess). When non-`NULL`, runs the
  full OU/GAM uncertainty propagation for `ndraws` draws and returns a
  **draws** product instead (see Details for the distinction). Requires
  `interpolation = "model"`.

- seed:

  Integer or `NULL`. RNG seed for reproducibility of draws.

- return:

  `"summary"` (default) or `"draws"`; relevant only when `ndraws` is
  supplied. `"summary"` collapses the draws to per-envelope central
  estimates (`median_*`, the draws' own central tendency — see
  `central`) plus credible-interval bounds (`lo_*`, `hi_*`) for the
  chosen `gap_uncertainty`, alongside the `deterministic` centre line.
  `"draws"` returns one row per (site \\\times\\ day \\\times\\ draw)
  with a `draw_id` column and the per-draw envelope value(s). See
  *Value*. Ignored in point mode (`ndraws = NULL`).

- interval:

  Credible interval width for `return = "summary"` (default `0.9`). The
  lower bound is the `(1 - interval)/2` quantile, the upper is the
  `1 - (1 - interval)/2` quantile.

- central:

  Central tendency for the `median_*` envelope columns when
  `return = "summary"`: `"median"` (default) or `"mean"` of the per-day
  draws. Because it is a summary *of the draws*, each `median_*` is
  coherent with its band (it lies inside `[lo_*, hi_*]` by construction)
  and depends on `ndraws`/`seed`. It is NOT the deterministic point
  estimate, which is reported separately in the `deterministic` column
  (and is the sole estimate when called with `ndraws = NULL`).

- grab_cv:

  Numeric scalar or named numeric vector of coefficients of variation
  for grab-sample measurement error. A scalar applies the same CV to all
  analytes; a named vector (e.g. `c(Cu = 0.1, pH = 0.02)`) applies
  analyte-specific CVs. Controls two uncertainty sources: S6 (anchor
  S-value spread at measured dates) and S7 (co-analyte normalisation
  spread). `NULL` (default) disables S6 and S7.

- ou_scale:

  Scale factor for the OU bridge envelope (default `1`). Multiplies
  \\\sigma^2\\ and \\\gamma\\ (marginal variance) without changing
  \\\theta\\ (correlation length); values \> 1 widen the between-grab
  uncertainty bands. Either a single number applied to every analyte,
  **or** a numeric vector named by impact tier to scale the tiers
  independently, e.g. `c(model = 1, bridge = 2.5)` (a tier not named
  falls back to `1`). Measure each tier's calibration with
  [`loo_anchor_coverage()`](https://vorpalvorpal.github.io/hydroSense/reference/loo_anchor_coverage.md)
  and raise the scale of an under-covered tier — the sparse-anchor
  **bridge** tier is typically the under-covered one, and a per-tier
  scale widens it without inflating the well-calibrated **model** tier.
  See
  [`vignette("daily-uncertainty")`](https://vorpalvorpal.github.io/hydroSense/articles/daily-uncertainty.md).

- kappa:

  Non-negative numeric (default `0.5`). Hydrological modulation
  exponent: the smoother's process variance is multiplied by
  \\\exp(\kappa \cdot z_h)\\ where \\z_h\\ is the standardised daily
  flow (capped at ±4 SD). `kappa = 0` disables hydrological modulation;
  larger values widen the credible band more aggressively across
  high-flow gaps. Requires a hydrology column in `df`.

- parallel:

  Logical (default `FALSE`). When `TRUE`, the per-**site** loop
  (interpolation + impact-model fit + residual-draw generation for each
  site) is parallelised via
  [`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html),
  which honours whatever
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  the caller has established. This is the embarrassingly-parallel level;
  the cross-site toxicity combine downstream is a single vectorised pass
  and runs once. Only speeds up multi-site `df` (one site = one task).
  Requires the **future.apply** package and a parallel plan, e.g.
  `future::plan(future::multisession, workers = 4)`; note workers need
  the **installed** package (use
  [`future::multicore`](https://future.futureverse.org/reference/multicore.html)
  to parallelise a `devtools::load_all()` session). In parallel mode
  `future.apply` manages a reproducible per-site L'Ecuyer-CMRG stream
  from `seed`, so draws differ from the sequential stream (point mode is
  deterministic and identical either way), but are themselves
  reproducible.

- couple_residuals:

  Logical (default `TRUE`). When `TRUE` and \>= 2 analytes have fitted
  residual smoothers, daily residual draws are correlated across
  analytes using the empirical anchor-residual correlation (see
  [`.anchor_residual_cor()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-anchor_residual_cor.md)).
  This widens the combined msPAF interval to reflect co-movement of
  co-toxicants on breach events while leaving per-analyte marginals
  unchanged. Set to `FALSE` to reproduce the pre-#32 independent-draw
  path exactly.

- gap_uncertainty:

  One of `"bracket"` (default), `"ignorable"`, or `"informative"`;
  relevant only in draws mode. In observation gaps the latent residual
  reverts to its marginal variance, widening the band — the honest
  posterior **only under ignorable (MAR) missingness**. Field sampling
  is often **informative (MNAR)**: gaps exist *because* the system was
  judged quiescent, so that band over-states gap uncertainty. The model
  cannot identify the mechanism from grabs, so the default **brackets**
  both: the *ignorable* (upper) envelope keeps the simulation-smoother
  draw; the *informative* (lower) envelope freezes the residual at its
  posterior mean on in-gap days (the fully-informative extreme — gaps
  perfectly predictable). The two are nested and coincide at observation
  days. `"ignorable"` / `"informative"` return only that envelope. See
  *Interpreting gap uncertainty* below. Reference: Rubin (1976)
  [doi:10.1093/biomet/63.3.581](https://doi.org/10.1093/biomet/63.3.581)
  .

- transform:

  `"pseudo_log"` (default) or `"additive"`. Controls the
  variance-stabilising transform for the daily impact residual smoother.
  `"pseudo_log"` applies `g = asinh(I / c)` with per-analyte scale
  `c = HC5` (issue \#15); `"additive"` keeps `g = I` (pre-#15 behaviour,
  no HC5 scaling). Forwarded to
  [`fit_target_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_target_model.md);
  ignored unless `interpolation = "model"`.

## Value

**Point mode** (`ndraws = NULL`): a tibble with one row per (site
\\\times\\ day) for days with sufficient analyte coverage; `mspaf` is
the deterministic daily estimate.

**Draws mode** (`ndraws > 0`, `return = "summary"`): one row per (site
\\\times\\ day) with the `deterministic` centre line and the envelope
columns for the chosen `gap_uncertainty`: `median_*`, `lo_*`, `hi_*` per
envelope (`informative` and/or `ignorable`, the central tendency per
`central` and the `interval` credible bounds), plus
`precautionary_lo`/`precautionary_hi` in `"bracket"` mode.

**Draws mode** (`ndraws > 0`, `return = "draws"`): one row per (site
\\\times\\ day \\\times\\ draw), with a `draw_id` integer column and the
per-draw msPAF value(s) `mspaf_ignorable` and/or `mspaf_informative` for
the chosen `gap_uncertainty`.

Common columns:

- `date`:

  Date of this daily estimate.

- `site_id`:

  Site identifier.

- `n_analytes_used`:

  SSD-eligible analytes contributing to msPAF.

- `dominant_analyte`:

  Analyte with the highest individual PAF.

- `max_paf`:

  PAF of the dominant analyte (proportion 0–1).

- `n_measured_analytes`:

  SSD-eligible analytes with a direct grab sample on this day (not
  interpolated).

- `days_since_last_sample`:

  Days since the most recent grab sample for any SSD-eligible analyte.

`"analyte_pafs"` (per-analyte PAF breakdown, re-keyed by date) and
`"ara_summary"` attributes are attached; retrieve them with
[`analyte_pafs()`](https://vorpalvorpal.github.io/hydroSense/reference/analyte_pafs.md)
and
[`ara_summary()`](https://vorpalvorpal.github.io/hydroSense/reference/ara_summary.md).
(`analyte_pafs` is now a flat attribute, not a list-column — issue
\#30.)

## Details

`mspaf_daily()` returns one of **two distinct products**, chosen by
`ndraws`; they answer different questions and should not be expected to
coincide:

- **Deterministic** (`ndraws = NULL`, default). A single grabs-exact
  daily line: the residual smoother is pinned to the measured grab
  samples, so the curve threads the reported lab values. *Pro:* fast,
  reproducible, the best single-number estimate; ideal as a sanity
  check. *Con:* carries no uncertainty and over-fits noisy grabs.

- **Draws** (`ndraws > 0`). A Monte-Carlo posterior propagating trend
  (GAM) and between-grab (OU bridge) uncertainty, plus grab measurement
  error when `grab_cv` is set. `return = "summary"` reports the central
  tendency plus a credible band; `return = "draws"` returns per-draw
  paths. *Pro:* honest uncertainty quantification. *Con:* the central
  estimate is seed/`ndraws`-dependent and, by Jensen's inequality on the
  bounded msPAF index, generally differs from the deterministic line.

The summary centre is a summary *of the draws* (issue \#42), so it
always lies within its own credible band; it is not the deterministic
line overlaid on a band built from a different posterior.

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
  [`fit_target_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_target_model.md)
  on the grab chemistry and the supplied `reference_model`, and predicts
  each toxicant's concentration between grabs as `reference + impact`,
  where the impact (the leachate-attributable increment) is modelled
  from hydrology and a persistent latent state but **never** from
  day-of-year. Co-analytes are forward-filled. Requires
  `reference_model`; see
  [`fit_target_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_target_model.md)
  for the method.

Below-detection values are treated as their detection-limit value for
interpolation purposes, matching the treatment in
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md).

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
  [`estimate_water_temp()`](https://vorpalvorpal.github.io/hydroSense/reference/estimate_water_temp.md))
  fills in days without a grab temperature. Grab-sample-day measurements
  take priority.

Set `require_temperature = FALSE` only for datasets that do not contain
ammonia.

## Interpreting gap uncertainty

In draws mode the width of the band across an observation gap is the
honest posterior **only under ignorable (MAR) missingness**. Where you
have external grounds that a gap was quiescent (informative/MNAR
sampling), the ignorable band over-states uncertainty there. The
`"bracket"` output gives both extremes: read the **informative** (lower)
envelope where you vouch the gap was quiet, the **ignorable** (upper)
envelope otherwise. The `precautionary_lo`/`precautionary_hi` columns
are the composite `[lo_informative, hi_ignorable]` — a **decision bound,
not a calibrated credible interval** (its coverage exceeds nominal and
is undefined). Applied blanket, the informative envelope under-covers
genuinely eventful gaps, so use it per-gap; automatic per-gap
conditioning on a continuous proxy is future work (#18). The
`deterministic` line is **not** a safe blanket alternative — it
under-states risk by Jensen's inequality (#39/#42).

## See also

[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md),
[`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md),
[`estimate_water_temp()`](https://vorpalvorpal.github.io/hydroSense/reference/estimate_water_temp.md),
[`get_silo_air_temp()`](https://vorpalvorpal.github.io/hydroSense/reference/get_silo_air_temp.md)

## Examples

``` r
# \donttest{
demo <- leachate_demo()
ds <- subset(demo, site_id == "downstream")
out <- mspaf_daily(ds, require_temperature = FALSE)
head(out[, c(
  "date", "site_id", "mspaf", "n_measured_analytes",
  "days_since_last_sample"
)])
#> # A tibble: 6 × 5
#>   date       site_id    mspaf n_measured_analytes days_since_last_sample
#>   <date>     <chr>      <dbl>               <int>                  <int>
#> 1 2024-01-15 downstream  93.3                   3                      0
#> 2 2024-01-16 downstream  93.3                   0                      1
#> 3 2024-01-17 downstream  93.3                   0                      2
#> 4 2024-01-18 downstream  93.3                   0                      3
#> 5 2024-01-19 downstream  93.3                   0                      4
#> 6 2024-01-20 downstream  93.3                   0                      5
# }
```
