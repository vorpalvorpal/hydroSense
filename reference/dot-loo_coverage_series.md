# Leave-one(-block)-out coverage of one residual series

Holds out each anchor (or each run of `block` consecutive anchors),
refits the smoother on the rest, and checks whether the held-out value
falls inside the predictive interval (state posterior + observation
noise). Plain LOO mainly tests the pinch zone; `block >= 2` probes
mid-gap width.

## Usage

``` r
.loo_coverage_series(
  anchor_dates,
  anchor_S,
  interval = 0.9,
  block = 1L,
  n_fit_min = 8L,
  scale = 1
)
```

## Arguments

- anchor_dates, anchor_S:

  The series.

- interval:

  Nominal coverage (default 0.9).

- block:

  Number of consecutive anchors held out per fold (default 1).

- n_fit_min, scale:

  Passed through to the smoother.

## Value

list(coverage, mean_width, n, n_held). `coverage`/`mean_width` are `NA`
when too few anchors remain to fit.
