# Leave-one-anchor-out coverage of the daily Kalman bands

Diagnoses whether the between-grab uncertainty bands produced by
[`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)
(`interpolation = "model"`) are calibrated on **your** data. For each
modelled analyte it holds out each grab anchor (or a run of `block`
consecutive anchors), refits the OU/Kalman smoother on the rest, and
checks whether the held-out value lands inside the nominal `interval`
predictive band. The result is reported **per analyte and per impact
tier**, because the two tiers are usually calibrated differently:

## Usage

``` r
loo_anchor_coverage(target_model, interval = 0.9, block = 1L, scale = 1)
```

## Arguments

- target_model:

  A fitted
  [`fit_target_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_target_model.md)
  object (class `"target_model"`).

- interval:

  Nominal coverage to test, in `(0, 1)`. Default `0.9`.

- block:

  Number of consecutive anchors held out per fold (default `1`).
  `block >= 2` probes mid-gap width rather than just the near-anchor
  pinch.

- scale:

  OU process-variance scale used while refitting, matching the
  `ou_scale` argument of
  [`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)
  (a single number here). Sweep it to find the value that brings a tier
  to nominal coverage.

## Value

A tibble with one row per analyte (`analyte`, `tier`, `n` = number of
anchors, `coverage`, `mean_width`) plus a final `"(pooled)"` row
(anchor-weighted). `coverage`/`mean_width` are `NA` for analytes with
too few anchors to refit.

## Details

- **`"model"`** analytes (enough detected anchors to fit a hydrology
  impact model) are normally close to nominal.

- **`"bridge"`** analytes (too sparse for their own fit, interpolated
  from a shared shape) tend to be **under-covered** — their bands are
  too narrow.

Use this to *measure*, then act with
[`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)'s
tier-aware `ou_scale`: if the pooled `bridge` coverage is below your
`interval`, raise the bridge tier's scale (e.g.
`ou_scale = c(model = 1, bridge = 2)`) and re-measure. See
[`vignette("daily-uncertainty")`](https://vorpalvorpal.github.io/hydroSense/articles/daily-uncertainty.md)
for the full workflow.

Coverage is a finite-sample estimate; with few anchors (especially in
the bridge tier) it is noisy, so prefer **pooling across multiple
sites** before choosing a scale.

## See also

[`fit_target_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_target_model.md),
[`mspaf_daily()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_daily.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# 1. Fit the daily target model (the same object mspaf_daily() builds).
tm <- fit_target_model(target_chem, reference_model = ref)

# 2. Measure calibration by tier at the default scale.
loo_anchor_coverage(tm, interval = 0.9)
#>   analyte  tier      n coverage mean_width
#>   <chr>    <chr> <int>    <dbl>      <dbl>
#>   Cu       model    38    0.895      0.71
#>   Pb       bridge   10    0.500      1.62   # under-covered
#>   (pooled) NA      194    0.85       ...

# 3. Find a bridge scale that reaches nominal, then act in the pipeline.
loo_anchor_coverage(tm, scale = 2)$coverage          # re-measure
daily <- mspaf_daily(target_chem, reference_model = ref,
                     interpolation = "model", ndraws = 200,
                     ou_scale = c(model = 1, bridge = 2))
} # }
```
