# Time-weighted chronic aggregation for any long-format value column

For each (focal date × monitoring feature × analyte) combination,
aggregates values from the preceding `window` (in days) using
exponential-decay temporal weighting and forward-step duration
weighting. This is the chronic exposure / response predictor used for
downstream calibration against biological community state.

## Usage

``` r
time_weighted_aggregate(
  df,
  focal_dates,
  tau = NULL,
  tau_units = NULL,
  window = NULL,
  window_units = NULL,
  summary = c("geom_mean", "arith_mean", "p90"),
  anchor_outside_window = TRUE,
  eps = 1e-09,
  return = c("summary", "draws"),
  interval = 0.9,
  central = c("median", "mean")
)
```

## Arguments

- df:

  Long-format data frame. Required columns:

  - `sample_id` (character) — unique sample identifier

  - `site_id` (character) — monitoring feature identifier

  - `datetime` (Date or POSIXct) — sample collection date/time

  - `analyte` (character) — analyte / index name

  - `value` (numeric) — concentration or index value

  - `detected` (logical) — whether value is a detected observation
    Optional: `imputed` (logical) — propagated to `n_imputed_in_window`.

- focal_dates:

  Either a `Date` vector (applied to all features in `df`) or a data
  frame with columns `focal_date` (Date) and `site_id` (character)
  specifying per-feature focal dates.

- tau:

  Exponential-decay rate parameter. Numeric or `units` object; bare
  numeric requires `tau_units`. Default `NULL` → 90 days. The effective
  half-life is `tau * log(2)` ≈ 62 days at the default. Choose to match
  the response timescale of the downstream biology (algae: days–weeks;
  macroinvertebrates: weeks–months; fish: months–years).

- tau_units:

  Character unit string for `tau` when it is bare numeric, e.g. `"d"`.
  Ignored when `tau` is a `units` object or `NULL`.

- window:

  Look-back window. Numeric or `units` object; bare numeric requires
  `window_units`. Default `NULL` → 365 days.

- window_units:

  Character unit string for `window` when it is bare numeric, e.g.
  `"d"`. Ignored when `window` is a `units` object or `NULL`.

- summary:

  Aggregation method: `"geom_mean"` (default), `"arith_mean"`, `"p90"`.
  See *Choosing `summary`* above.

- anchor_outside_window:

  Logical (default `TRUE`). If `TRUE`, the most recent sample *before*
  the window start is included as an anchor to provide duration
  weighting at the leading edge of the window. Its duration is clipped
  to start at `window_start` so it only contributes the in-window
  portion of its interval.

- eps:

  Small positive guard added inside the log for geometric mean to avoid
  `log(0)`. Default `1e-9`.

- return:

  Output mode for draw-carrier input (see
  [`summarise_draws()`](https://vorpalvorpal.github.io/hydroSense/reference/summarise_draws.md)).
  `"summary"` (default) collapses posterior draws to a central estimate
  plus a credible interval (`value`, `value_lower`, `value_upper`,
  `n_draws`); `"draws"` returns the raw per-draw chronic rows
  (`draw_id 1..N`). For point (non-draw) input both modes return
  byte-identical output with no interval columns. Draws are paired
  across time by `draw_id` (an index-pairing approximation that assumes
  temporal independence; see the OU/Kalman smoothing roadmap item).

- interval:

  Width of the credible interval when `return = "summary"`. Default
  `0.90` (5th/95th percentile bounds).

- central:

  Central-tendency statistic when `return = "summary"`: `"median"`
  (default) or `"mean"`.

## Value

A long-format tibble with columns:

- `focal_date` (Date)

- `site_id` (character)

- `sample_id` (character) — synthetic key
  `"chronic_<focal_date>_<site>"`

- `analyte` (character)

- `value` (numeric) — chronic time-weighted value

- `detected` (logical) — always `TRUE`

- `n_samples_in_window` (integer)

- `n_imputed_in_window` (integer) — count of `imputed == TRUE` samples
  contributing; 0 if `imputed` column absent from `df`

## Details

This function is value-agnostic: pass any long-format data frame with
`analyte` and `value` columns and it will compute one time-weighted
value per (focal_date × site_id × analyte). Use it to aggregate raw
chemistry (`analyte` = chemical species, `summary = "geom_mean"`) or
per-sample msPAF values (`analyte = "msPAF"`, `summary = "arith_mean"`)
into a chronic predictor.

**Forward-step duration weighting** treats each sample as representing
the period from its collection date to the next sample's date (or to
`focal_date` for the most recent sample). This corrects for pulse-biased
sampling where storm events are sampled more frequently than base-flow
periods, which would otherwise over-weight episodic concentrations.

**Exponential-decay temporal weighting** assigns higher weight to recent
samples, with a half-life of approximately `tau * log(2)` days (at the
default of 90 d).

These two components are the only weighting scheme offered, by design:
forward-step duration weighting is the minimal correction for irregular
/ pulse-biased sampling, and exponential decay is the standard memory
kernel for a community integrating a fluctuating exposure (one
interpretable parameter, `tau`, tied to the target biology's response
time). Richer kernels would add parameters without a defensible way to
fit them from routine monitoring data, so they are left out — tune `tau`
rather than swapping the kernel.

**Combined weight** for sample *i* is
`w_i = Δt_i × exp(-(focal_date - midpoint_i) / tau_d)`, where
`midpoint_i` is the midpoint of sample *i*'s representative interval.

**Choosing `summary`:**

- `"geom_mean"`: max-likelihood central tendency for log-normal data.
  Use for chemistry concentrations. Caveat: with highly pulsed exposure
  and a non-linear SSD response, this underestimates the time-averaged
  risk metric (Jensen's inequality on the upper tail of the SSD).

- `"arith_mean"`: weighted arithmetic mean. Use for bounded indices like
  msPAF percentages, or for chemistry when comparison against an
  arithmetic-mean compliance trigger is wanted.

- `"p90"`: 90th percentile (duration-weighted empirical CDF). Diagnostic
  for "what's the upper-end exposure the community sees most of the
  time."

For chronic msPAF specifically, the recommended pipeline is
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)
on per-sample chemistry, then `time_weighted_aggregate()` on the
resulting msPAF rows with `summary = "arith_mean"`. This computes the
time-averaged msPAF, which integrates the toxic response over time
(consistent with how biological communities respond to fluctuating
exposure) — rather than the msPAF computed at a single time-averaged
chemistry, which can substantially under-state risk for pulsed
exposures.

## Examples

``` r
# Chronic (time-weighted geometric-mean) chemistry for one downstream
# analyte at two focal dates, from the bundled demo data.
cu <- subset(leachate_demo(), site_id == "downstream" & analyte == "Cu")
time_weighted_aggregate(
  cu,
  focal_dates = as.Date(c("2024-06-01", "2024-12-01")),
  tau = 90, tau_units = "d", window = 365, window_units = "d",
  summary = "geom_mean"
)
#> # A tibble: 2 × 8
#>   focal_date site_id    sample_id     analyte value detected n_samples_in_window
#>   <date>     <chr>      <chr>         <chr>   <dbl> <lgl>                  <int>
#> 1 2024-06-01 downstream chronic_2024… Cu       7.56 TRUE                       3
#> 2 2024-12-01 downstream chronic_2024… Cu       7.50 TRUE                       6
#> # ℹ 1 more variable: n_imputed_in_window <int>
```
