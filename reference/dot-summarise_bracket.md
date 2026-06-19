# Summarise per-draw bracket AmsPAF into a tidy per-day frame

Summarise per-draw bracket AmsPAF into a tidy per-day frame

## Usage

``` r
.summarise_bracket(
  draws_df,
  interval = 0.9,
  central = "median",
  gap_uncertainty = "bracket"
)
```

## Arguments

- draws_df:

  Long per-draw frame with `date`, `site_id`, `draw_id` and the envelope
  value columns needed for `gap_uncertainty`: `amspaf_ignorable` and/or
  `amspaf_informative`.

- interval:

  Credible-interval width (default 0.9). The lower bound is the
  `(1 - interval)/2` quantile, the upper the `1 - (1 - interval)/2`
  quantile (type-7, matching the package's other draw summaries).

- central:

  Per-day central tendency: `"median"` (default) or `"mean"` of the
  draws (issue \#42 – the centre is a summary of the draws themselves).

- gap_uncertainty:

  `"bracket"` (both envelopes + precautionary), `"ignorable"` (ignorable
  columns only), or `"informative"` (informative columns only).

## Value

A tibble keyed by (`date`, `site_id`) with the envelope columns for the
requested mode: `median_*`, `lo_*`, `hi_*` per envelope, and
`precautionary_lo`/`precautionary_hi` in `"bracket"` mode.
