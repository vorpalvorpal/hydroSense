# Log-linear daily interpolation of grab samples, with a span cap

Fills a grid of days from a point's own grab samples ("anchors").
Between two anchors the value is interpolated linearly in `log(value)`
(so the midpoint of a span is the geometric mean of its ends), or
linearly in `value` when `log = FALSE`. This is the toxicant
interpolation used by
[`mspaf_pipeline()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_pipeline.md)
with `ara_form = "added_hazard"`.

## Usage

``` r
interpolate_loglinear(
  dates,
  values,
  detected,
  grid,
  span_cap = 120,
  log = TRUE
)
```

## Arguments

- dates:

  `Date` vector of anchor (sample) dates. Must be unique; any order.

- values:

  Numeric anchor values, one per date.

- detected:

  Logical, one per date. `FALSE` marks a below-detection anchor; `NA` is
  treated as `TRUE`.

- grid:

  `Date` vector of the days to fill, in any order.

- span_cap:

  Maximum anchor-to-anchor span (and hold after the last anchor) in
  days. Default 120. `Inf` disables the cap.

- log:

  Logical. Interpolate in `log(value)` (default) or in `value`. With
  `log = TRUE` anchors with an `NA` or non-positive value are dropped;
  with `log = FALSE` only `NA` anchors are dropped.

## Value

A data frame with one row per element of `grid`, in grid order, and
columns `date` (= `grid`), `value` (numeric) and `detected` (logical).
With no usable anchor every row is `NA`.

## Span cap

A day is only filled when the anchors that bracket it are close enough
together. For a grid day `g`, with `x_i` the last anchor on or before
`g`:

- there is no anchor on or before `g`: the day is uncovered (leading
  days are never back-filled);

- `g` is an anchor date: the anchor's own value and detection flag;

- `g` is after the last anchor: covered only while
  `g - x_last <= span_cap`, holding the last anchor's value and flag;

- otherwise, with `x_j` the next anchor: covered only when
  `x_j - x_i <= span_cap`, and then interpolated between the two
  anchors.

Uncovered days have `value = NA` and `detected = NA`. `span_cap = Inf`
disables the cap, so every day after the first anchor is covered.

## Below-detection anchors

A below-detection anchor (`detected = FALSE`) enters the interpolation
at the value supplied for it (the caller decides: the detection limit,
half of it, ...). An interpolated day is flagged `detected = FALSE` only
when **both** bracketing anchors are below detection; one detected
anchor is enough to make the day detected. An anchor date, or a held day
after the last anchor, carries that anchor's own flag.

## See also

[`mspaf_pipeline()`](https://vorpalvorpal.github.io/hydroSense/reference/mspaf_pipeline.md)
(`ara_form = "added_hazard"`)

## Examples

``` r
d0 <- as.Date("2024-01-01")
interpolate_loglinear(
  dates = d0 + c(0, 10, 200), values = c(2, 8, 4),
  detected = c(TRUE, FALSE, TRUE), grid = d0 + c(-1, 0, 5, 10, 100, 330)
)
#>         date value detected
#> 1 2023-12-31    NA       NA
#> 2 2024-01-01     2     TRUE
#> 3 2024-01-06     4     TRUE
#> 4 2024-01-11     8    FALSE
#> 5 2024-04-10    NA       NA
#> 6 2024-11-26    NA       NA
```
