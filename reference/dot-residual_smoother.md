# Fit + smooth the residual for one analyte over a (clipped) daily grid

Clips `target_dates` to the analyte's grab span
`[first_anchor, last_anchor]` (no extrapolation beyond the grabs), fits
the OU parameters, builds the KFAS model with optional
hydrology-modulated process variance, and returns the smoothed
posterior. For degenerate series returns a flat mean with zero variance.

## Usage

``` r
.residual_smoother(
  anchor_dates,
  anchor_S,
  target_dates,
  z_hydro = NULL,
  kappa = 0.5,
  n_fit_min = 8L,
  scale = 1,
  r_vec = NULL
)
```

## Arguments

- anchor_dates, anchor_S:

  Observed anchors.

- target_dates:

  Candidate daily grid (will be clipped to the grab span).

- z_hydro:

  Optional standardised hydro feature aligned to `target_dates` (same
  length); drives `q_mult = exp(kappa * z)`. `NULL` -\> no modulation.

- kappa:

  Hydrology process-variance sensitivity (default 0.5).

- n_fit_min, scale:

  Passed to
  [`.estimate_ou_kalman_params()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-estimate_ou_kalman_params.md).

- r_vec:

  Per-anchor observation variance (default tiny).

## Value

list(grid_dates, params, model, mean, var). `model` is `NULL` for the
degenerate path.
