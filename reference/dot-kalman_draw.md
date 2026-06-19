# Draw coherent residual trajectories via the simulation smoother

Draw coherent residual trajectories via the simulation smoother

## Usage

``` r
.kalman_draw(model, nsim = 1L)
```

## Arguments

- model:

  A KFAS model from
  [`.build_kalman_model()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-build_kalman_model.md).

- nsim:

  Number of draws.

## Value

Numeric matrix `[n_grid x nsim]`. Call
[`set.seed()`](https://rdrr.io/r/base/Random.html) beforehand for
reproducibility.
