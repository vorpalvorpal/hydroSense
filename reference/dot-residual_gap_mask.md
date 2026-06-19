# Flag grid days where the residual posterior variance has ballooned in a gap

The simulation-smoother posterior variance is small at (and near) the
grab anchors and inflates between them, peaking in the middle of long
gaps (Brownian-bridge shape). A grid day counts as "in-gap" when its
variance has risen meaningfully above the at-anchor floor — i.e. more
than a fraction `tol` of the way from that floor to the series' peak
variance. This is the set of days on which the informative envelope
(residual frozen at its posterior mean) may differ from the ignorable
envelope (residual drawn); everywhere else the two coincide. Observation
days, a densely-sampled series with no real gaps, and a degenerate
(zero-variance) series all return `FALSE` throughout.

## Usage

``` r
.residual_gap_mask(sm, anchor_dates, tol = 0.05)
```

## Arguments

- sm:

  A smoother list with `grid_dates` and `var` (e.g. from
  [`.residual_smoother()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-residual_smoother.md)
  or the per-analyte scaffold).

- anchor_dates:

  Grab/anchor dates; used to locate the at-anchor floor.

- tol:

  Fraction of the floor-to-peak variance range above the floor at which
  a day is declared in-gap (default 0.05).

## Value

Logical vector along `sm$grid_dates`; `logical(0)` for an empty grid.

## Details

The floor is taken at the anchor positions (the grab days carry only the
tiny observation-noise variance); using a *relative* threshold against
the floor-to-peak range makes the test scale-free, so it collapses to
"no gaps" automatically when sampling is dense.
