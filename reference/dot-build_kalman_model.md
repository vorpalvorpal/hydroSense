# Build a KFAS state-space model for one residual series on a grid

**Invariant:** `grid_dates` must be *consecutive days* (step = 1 day).
The transition `phi = exp(-theta)` assumes a unit time step between
adjacent grid points, so a sparse/irregular grid would mis-state the
correlation length.

## Usage

``` r
.build_kalman_model(
  grid_dates,
  anchor_dates,
  anchor_S,
  theta,
  gamma,
  r_vec = NULL,
  q_mult = NULL
)
```

## Arguments

- grid_dates:

  Sorted, consecutive-day Date vector: the grid to predict on.

- anchor_dates, anchor_S:

  Observed anchors (must fall within `grid_dates`).

- theta, gamma:

  OU parameters from
  [`.estimate_ou_kalman_params()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-estimate_ou_kalman_params.md).

- r_vec:

  Per-anchor observation variance (length = \#anchors), or `NULL` (-\> a
  tiny value relative to gamma, so the smoother pins to the anchors).

- q_mult:

  Per-grid-day positive multiplier on the process variance (hydrology
  modulation), length `length(grid_dates)`, or `NULL` (-\> 1).

## Value

A KFAS `SSModel`, or `NULL` if `gamma <= 0` or no anchors on the grid.
The grid is stored in `attr(., "grid_dates")`.
