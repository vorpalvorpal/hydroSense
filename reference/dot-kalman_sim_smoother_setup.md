# Precompute the DK simulation smoother gain matrix for a KFAS model

Builds the gain matrix `L` (n_grid × m, standardised units) such that
`L %*% y_anch` recovers the KFS-smoothed state when `y_anch` is placed
at the anchor positions. This relies on the linearity of the KFS
smoother in `y` when `a1 = 0`: the smoother output is a linear function
of the observations, so it can be precomputed column-by-column by
running KFS on unit vectors (Durbin & Koopman 2002, Biometrika 89(3),
§4).

## Usage

``` r
.kalman_sim_smoother_setup(model)
```

## Arguments

- model:

  A KFAS `SSModel` from
  [`.build_kalman_model()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-build_kalman_model.md).

## Value

A list with components:

- `L`:

  `[n_grid × m]` gain matrix (standardised units).

- `pos`:

  Integer indices of anchor positions in the grid.

- `x_hat`:

  KFS posterior mean in un-standardised units (length n_grid).

- `phi`:

  Scalar AR(1) coefficient.

- `q_sd_vec`:

  Standardised process noise SD per grid day (length n_grid).

- `h_sd_vec`:

  Standardised observation noise SD at anchor positions (length m).

- `resid_scale`:

  Un-standardising factor (`sqrt(gamma)`).

- `n_grid`:

  Number of grid days.

## Details

The returned `x_hat` is in **original (un-standardised) units** to match
the output of
[`.kalman_smooth()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-kalman_smooth.md).

## References

Durbin J, Koopman SJ (2002). A simple and efficient simulation smoother
for state space time series analysis. *Biometrika* **89**(3), 603–616.
<https://doi.org/10.1093/biomet/89.3.603>
