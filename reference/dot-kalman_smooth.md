# Posterior (smoothed) mean and variance of the residual on the grid

Posterior (smoothed) mean and variance of the residual on the grid

## Usage

``` r
.kalman_smooth(model)
```

## Arguments

- model:

  A KFAS model from
  [`.build_kalman_model()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-build_kalman_model.md).

## Value

list(mean, var), each length `length(grid_dates)`.
