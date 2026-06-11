# Build the state-space residual smoother for one analyte

Uses the WQ residual `d` (`m$d_anchors`, WQ tier) or the impact residual
`S` (`m$anchors`, impact/bridge tier). The daily grid is clipped to the
analyte's grab span; hydrology modulates the process variance. See
[`.residual_smoother()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-residual_smoother.md).

## Usage

``` r
.analyte_residual_smoother(
  m,
  target_model,
  qdates,
  kappa = 0.5,
  scale = 1,
  r_vec = NULL
)
```
