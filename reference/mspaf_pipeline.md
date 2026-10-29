# Run the full chronic msPAF pipeline end to end

Orchestrates the four steps a chronic msPAF analysis usually chains by
hand: optionally fit a Bayesian imputation model
([`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md)),
fit the reference/impact model
([`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md)),
compute the daily multi-substance PAF
([`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md))
for the target site, then aggregate the daily series into a **chronic**
time-weighted msPAF
([`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md)).
The fitted imputation model is threaded into **both** downstream
chemistry steps so that below-detection and entirely-absent analytes are
filled before the toxicity calculation.

## Usage

``` r
mspaf_pipeline(
  target,
  reference = NULL,
  hydro = NULL,
  impute = TRUE,
  impute_on = c("reference", "target"),
  required_vars = c("pH", "EC"),
  impute_groups = NULL,
  impute_seed = NULL,
  reference_args = list(),
  daily_args = list(),
  chronic = TRUE,
  focal_dates = NULL,
  focal_by = "day",
  tau = NULL,
  window = NULL,
  chronic_summary = "arith_mean",
  also_total = FALSE,
  ara = c("static", "contemporaneous", "none"),
  contributions = FALSE,
  ara_form = c("excess", "added_hazard"),
  min_coverage = 0.2
)
```

## Arguments

- target:

  Long-format target-site chemistry (same schema as
  [`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)'s
  `df`): `sample_id`, `site_id`, `datetime`, `analyte`, `value`,
  `detected`, and a unit column (`units.analyte`, or sampleTidy's
  `units`). Analyte names may be hydroSense's or sampleTidy's (see
  [`analyte_aliases()`](https://vorpalvorpal.github.io/hydroSense/reference/analyte_aliases.md)).

- reference:

  Long-format reference-site chemistry passed to
  [`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md).
  `NULL` (default) runs the pipeline without a reference site: no
  reference model is fitted (the `reference_model` attribute is `NULL`),
  `ara` defaults to `"none"`, daily interpolation defaults to `"linear"`
  and `impute_on` defaults to `"target"`. See Details.

- hydro:

  Optional hydrology series forwarded to
  [`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md).

- impute:

  Logical. Fit and apply a Bayesian imputation model? Default `TRUE`.
  When `FALSE`, no imputation model is fitted and `NULL` is passed
  downstream (identical to chaining the functions by hand without one).

- impute_on:

  Which chemistry to fit the imputation model on: `"reference"`
  (default, matching the documented design) or `"target"`. With
  `reference = NULL` the default becomes `"target"`; an explicit
  `"reference"` is an error of class `hydroSense_impute_needs_reference`
  when `impute = TRUE`.

- required_vars:

  Passed to
  [`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md).
  Default `c("pH", "EC")`.

- impute_groups:

  Optional `groups` for
  [`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md);
  `NULL` (default) uses
  [`leachate_impute_groups()`](https://vorpalvorpal.github.io/hydroSense/reference/leachate_impute_groups.md).

- impute_seed:

  Optional integer seed for the imputation fit (forwarded to
  [`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md))
  so a stochastic fit is reproducible. The default `"marginal"` method's
  GAM fit is deterministic; the seed matters for the Stan-based methods
  (`"factor"` and the brms methods).

- reference_args, daily_args:

  Named lists of additional arguments forwarded to
  [`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md)
  and
  [`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)
  respectively. Do not include arguments the orchestrator sets itself
  (`reference`/`hydro`/ `imputation_model` for the reference fit;
  `df`/`reference_model`/ `imputation_model` for the daily call).
  `daily_args` may override the `interpolation = "model"` default. When
  `chronic = TRUE` and `daily_args` requests draws (`ndraws`),
  `return = "draws"` is forced so per-draw rows can be propagated into
  the chronic summary; the `gap_uncertainty`, `interval` and `central`
  values from `daily_args` (or
  [`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)'s
  defaults `"bracket"`, `0.9`, `"median"`) shape the chronic envelope.

- chronic:

  Logical. Aggregate the daily series into chronic time-weighted msPAF
  as a final step? Default `TRUE`. When `FALSE`, the daily
  [`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)
  result is returned unchanged (with model attributes only).

- focal_dates:

  Focal dates for the chronic aggregation, passed to
  [`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md).
  `NULL` (default) derives a sequence spanning the daily date range via
  [`expand_focal_dates()`](https://vorpalvorpal.github.io/hydroSense/reference/expand_focal_dates.md)
  (spacing `focal_by`). May also be a `Date` vector or a
  `focal_date`/`site_id` data frame.

- focal_by:

  Spacing for the derived `focal_dates` when `focal_dates` is `NULL`;
  passed to
  [`expand_focal_dates()`](https://vorpalvorpal.github.io/hydroSense/reference/expand_focal_dates.md)'s
  `by` ([`base::seq.Date()`](https://rdrr.io/r/base/seq.Date.html)
  semantics): `"day"` (default), `7` (weekly), `"week"`, `"month"`, etc.
  Ignored when `focal_dates` is supplied.

- tau:

  Exponential-decay parameter in days for the chronic kernel, forwarded
  to
  [`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md).
  `NULL` (default) uses its default of 90 days.

- window:

  Look-back window in days for the chronic aggregation, forwarded to
  [`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md).
  `NULL` (default) uses its default of 365 days.

- chronic_summary:

  Aggregation method for the chronic step, passed to
  [`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md)'s
  `summary`. Default `"arith_mean"` (the recommended aggregation for
  bounded msPAF percentages).

- also_total:

  Logical (default `FALSE`). Point mode only. Passed to
  [`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md),
  which adds a total-concentration (no ARA) `mspaf_total` daily series
  alongside the ARA one; the chronic step aggregates it exactly as it
  aggregates `mspaf` and returns it as `value_total`. The static ARA
  reference is `daily_args$reference` (the `reference` argument here is
  the reference-model fit data); without it `value_total` equals
  `value`. An error in draws mode (`daily_args$ndraws` set). Composes
  with every `ara` mode: `value_total` is always the no-ARA value.

- ara:

  Area-reference adjustment (ARA) applied by
  [`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md):

  - `"static"`: static ARA against the per-analyte summary of
    `daily_args$reference` (a raw reference chemistry frame, see
    [`prepare_reference()`](https://vorpalvorpal.github.io/hydroSense/reference/prepare_reference.md)),
    which is then required.

  - `"contemporaneous"`: the reference model this pipeline fits is
    passed to
    [`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)
    as `reference`, so each day is adjusted by the reference model's
    resolution for that day (tier-1 match, GAM prediction or static
    fallback; see
    [`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md)).
    The model is fitted once and reused. `daily_args$reference` must not
    be given.

  - `"none"`: no ARA. `daily_args$reference` must not be given.

  With `reference = NULL`, `"contemporaneous"` (and `"static"` without
  `daily_args$reference`) is an error of class
  `hydroSense_ara_needs_reference`.

  When `ara` is not supplied the behaviour is unchanged from earlier
  versions: static ARA if `daily_args$reference` is given, otherwise
  none. Conflicting combinations are errors raised before any model is
  fitted.

- contributions:

  Logical (default `FALSE`). Point mode only. When `TRUE`, the result
  carries attribute `"contributions"`: the chronic msPAF apportioned to
  analytes and mode-of-action groups, with columns `site_id`,
  `focal_date`, `analyte`, `group` and `contribution` (%). Each
  analyte's daily contribution series (from
  [`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md);
  zero on days it did not enter the msPAF) is aggregated with the same
  time weights, window and rows as the chronic `value`, so the
  contributions sum to `value` for every (site, focal date). With
  `also_total = TRUE` a `variant` column (`"ara"` / `"total"`)
  distinguishes the contributions to `value` and `value_total`. Requires
  `chronic_summary = "arith_mean"` (the only additive aggregation) when
  `chronic = TRUE`. An error of class `hydroSense_contributions_draws`
  in draws mode (`daily_args$ndraws` set). With `chronic = FALSE` the
  daily breakdown from
  [`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)
  is attached to the daily result instead. For
  `ara_form = "added_hazard"` the definition differs; see the section
  *Added-hazard form*.

- ara_form:

  The form of the area-reference adjustment:

  - `"excess"` (default): the existing behaviour described above, where
    the reference is subtracted from each analyte's concentration before
    the msPAF is formed.

  - `"added_hazard"`: the reference is compared at the msPAF level, as
    the added hazard of the target over the reference (see the section
    *Added-hazard form*). It has its own engine: no model is fitted and
    nothing is imputed.

- min_coverage:

  Used only by `ara_form = "added_hazard"`: the smallest share (0 to 1)
  of the time weights that covered days must carry for a chronic value
  to be reported. Default `0.2`.

## Value

When `chronic = TRUE` (default), the chronic msPAF frame, one row per
(focal date × site):

- **point mode** (`ndraws` absent): columns `focal_date`, `site_id`,
  `sample_id`, `analyte`, `value`, `detected`, `n_samples_in_window`,
  `n_imputed_in_window` (the number of daily rows in the window with no
  grab sample, i.e. model-predicted or forward-filled; a single chronic
  `value`, no interval columns), plus `value_total` when
  `also_total = TRUE`.

- **draws mode** (`ndraws` set): columns `focal_date`, `site_id` plus
  the bracket envelope summary (`median_*`/`lo_*`/`hi_*` per envelope,
  and `precautionary_lo`/`precautionary_hi` in `"bracket"` mode). The
  chronic frame carries the daily frame it was built from as attribute
  `"daily"`. When `chronic = FALSE`, the
  [`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)
  result is returned unchanged. In all cases the fitted models are
  attached as attributes `"reference_model"` and `"imputation_model"`
  (the latter `NULL` when `impute = FALSE` or the fit was skipped). Note
  that most dplyr verbs drop attributes, so read them off the returned
  frame directly. For `ara_form = "added_hazard"` see the section
  *Added-hazard form*.

## Details

The chronic step "smears" the daily msPAF line with an exponential-decay
memory kernel: because the daily series is equally spaced, forward-step
duration weighting is uniform (Δt = 1 day) and the chronic value is
shaped only by the `tau` decay. Set `chronic = FALSE` to stop at the
daily series and reproduce the previous (daily) return.

**Imputation is on by default here.** The individual functions treat
imputation as opt-in (you must fit a model and pass it). This
orchestrator flips that default: with `impute = TRUE` it fits an
imputation model and applies it throughout. This is convenient but **not
free of consequences** — imputing below-detection cells (replacing
detection-limit values with modelled sub-limit estimates) and
fabricating entirely-absent analytes generally *moves* the msPAF result.
Whether that is appropriate is a missingness (MAR/MNAR) modelling
judgement, not a neutral default. Set `impute = FALSE` to reproduce the
non-imputed behaviour of calling
[`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md)
and
[`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)
directly.

With `reference = NULL` (no reference site of the right water type) no
reference model is fitted. `daily_args$interpolation = "model"` is then
an error (class `hydroSense_model_needs_reference`); the default is
`"linear"`. All such argument errors are raised before any model is
fitted.

## Added-hazard form

With `ara_form = "added_hazard"` each point (every target site, and the
reference) is built the same way from its own samples, on a daily grid.

**Input.** `target` and `reference` are raw long frames in the same
format (`sample_id`, `site_id`, `datetime`, `analyte`, `value`,
`detected`, a unit column). A below-detection row has `detected = FALSE`
and its detection limit (LOR) as `value`. Analyte names go through
[`analyte_aliases()`](https://vorpalvorpal.github.io/hydroSense/reference/analyte_aliases.md).
Toxicants are the analytes with SSD parameters for the given
`daily_args$analyte_metadata` / `method` / `guideline_dir`; the rest are
co-analytes. Rows are reduced to one value per (site, analyte, day): the
mean of the detected results when there is one, otherwise the smallest
LOR (below detection). `reference` must hold exactly one `site_id`
(error class `hydroSense_reference_sites`).

**Grid.** Daily, from the earliest sample date of the targets and the
reference to the latest target sample date plus the span cap (the latest
date itself when the cap is infinite), unless `daily_args$start` / `end`
are given.

**Point build.** Toxicants are interpolated with
[`interpolate_loglinear()`](https://vorpalvorpal.github.io/hydroSense/reference/interpolate_loglinear.md)
and the span cap `daily_args$span_cap` (default 120 days): days in a gap
between samples longer than the cap, more than the cap after the last
sample, or before the first sample are uncovered. Below-detection
results are anchored at the LOR for a target and at LOR/2 for the
reference. The co-analytes are interpolated without a cap and held
constant beyond the first and last sample (detected positive results
only; one result is a constant): pH and temperature linearly, DOC,
hardness, Ca and Mg linearly in log. Hardness is the measured value
where it exists, else `2.497 Ca + 4.118 Mg` on days with both.
`daily_args$temperature` (a `datetime`/`value` series, per site only
when it has a `site_id` column) replaces every point's own temperature
results when given; it is interpolated linearly between its days and
held beyond its ends. Each toxicant is normalised with the point's own
co-analytes using the metadata normalisation formulas, as
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)
does; a day lacking a required co-analyte leaves that analyte out of
that point on that day.

**Detection-limit rules and the gap-filler**, per target, analyte and
day: where the target is covered and the reference is not, the reference
takes its static summary
`s = max(geometric mean of its results with below-detection results at LOR/2, median(LOR of its below-detection results) / 2)`
(the geometric mean alone when it has none), normalised with the
reference's co-analytes that day. An analyte the reference never
measured is dropped. When the target is below detection that day (both
bracketing samples below detection) its value is set to the reference's,
so the analyte adds no hazard. An analyte enters a day only when both
sides have a value.

**Daily added hazard.** `msPAF` combines the normalised concentrations
by concentration addition within mode-of-action group
(`TU = sum(C / HC50)`, mean slope) and independent action across groups.
With `msT` and `msR` the target's and reference's msPAF on the same
analytes, `D = max((msT - msR) / (1 - msR), 0)` (0 when `msR` is 1).
Without a reference (`ara = "none"`) below-detection target days enter
at 0 and `D = msT`. On this path `ara` defaults to `"contemporaneous"`
with a reference and `"none"` without; `"static"`, `"contemporaneous"`
without a reference, `impute = TRUE`, a `daily_args$interpolation` other
than `"loglinear"` and a `chronic_summary` other than `"arith_mean"` are
errors of class `hydroSense_ara_conflict`, raised before anything is
computed (`daily_args$interpolation = "loglinear"` with
`ara_form = "excess"` is the same error). An unspecified `impute` means
`FALSE` here. `hydro`, `reference_args` and the `daily_args` entries
other than `interpolation`, `span_cap`, `temperature`, `start`, `end`,
`analyte_metadata`, `method`, `guideline_dir` and `conc_units` are
ignored.

**Time-weighting and coverage.** For a focal date the days of age
`a = 0..window` carry weight `exp(-a / tau)` (`tau` 90 and `window` 365
days by default). A day is covered when `D` exists. The chronic `value`
is `100 * sum(w D) / sum(w)` over the covered days; `coverage` is the
covered share of all the weights. `value` is `NA` when
`coverage < min_coverage` or nothing is covered. `value_total` (with
`also_total = TRUE`) is the same time-weighting (same weights, covered
days and `NA` gate) of `msT0`, the target's own msPAF with no reference:
its normalised values with below-detection analyte-days at concentration
0, on the same analytes as that day's `D`. Non-detects are therefore not
counted at the reference's level in the total.

**Contributions.** Per analyte and day the added single-substance PAF is
`max(PAF_i(T) - PAF_i(R), 0)`, with
`PAF_i(C) = pnorm(log10(C / HC50_i) / sigma_i)`. The day's `D` is
apportioned in proportion to it (all zero when no analyte adds any), and
the shares are time-weighted with exactly the weights and denominator of
`value`, so the contributions of a focal date sum to `value` and are
`NA` where it is. A target below detection gives that analyte exactly
zero. The `"total"` variant apportions `msT0` by `PAF_i(T0)` (`T0` the
target's values with below-detection days at 0, so a below-detection
analyte-day contributes 0) the same way and sums to `value_total`. A
group's contribution is the sum over its analytes.

**Output.** One row per (target site, focal date): `focal_date`,
`site_id`, `sample_id` (`"chronic_<focal_date>_<site_id>"`), `analyte`
(`"msPAF"`), `value`, `value_total` (with `also_total`), `detected`,
`n_samples_in_window` (target sample dates with a toxicant result in the
window), `n_imputed_in_window` (covered days in the window that are not
sample dates) and `coverage`. Attribute `"daily"` holds `date`,
`site_id`, `mspaf` (`100 D`), `mspaf_total` (`100 msT0`),
`mspaf_reference` and `n_analytes`; attribute `"sample_days"` has one
row per target sample date (`site_id`, `date`, `paired`, `TRUE` when the
reference has a toxicant result that day). With `chronic = FALSE` the
daily frame is returned, carrying `"sample_days"` (and the daily
`"contributions"`, by `date`, when requested).

## See also

[`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md),
[`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md),
[`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md),
[`time_weighted_aggregate()`](https://vorpalvorpal.github.io/hydroSense/reference/time_weighted_aggregate.md),
[`expand_focal_dates()`](https://vorpalvorpal.github.io/hydroSense/reference/expand_focal_dates.md),
[`interpolate_loglinear()`](https://vorpalvorpal.github.io/hydroSense/reference/interpolate_loglinear.md)

## Examples

``` r
if (FALSE) { # \dontrun{
demo <- leachate_demo()
out <- mspaf_pipeline(
  target = subset(demo, site_id == "downstream"),
  reference = subset(demo, site_id == "reference"),
  daily_args = list(require_temperature = FALSE, ndraws = 50L, seed = 1L)
)
} # }
```
