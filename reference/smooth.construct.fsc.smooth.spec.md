# Factor smooth with a canonical null-space basis (`bs = "fsc"`)

A drop-in variant of mgcv's factor-smooth interaction (`bs = "fs"`, see
[mgcv::factor.smooth](https://rdrr.io/pkg/mgcv/man/factor.smooth.html))
whose model matrix is stable under tiny (ulp-level) changes in the
covariate. Use it as `s(x, fac, bs = "fsc", k = ...)`; arguments `xt`,
`m` and `k` behave as for `"fs"`.

## Usage

``` r
# S3 method for class 'fsc.smooth.spec'
smooth.construct(object, data, knots)
```

## Arguments

- object:

  A smooth specification object of class `"fsc.smooth.spec"`, as
  produced by `s(..., bs = "fsc")`.

- data:

  A list or data frame containing the covariate and the factor.

- knots:

  Optional list of knots, passed to the base smooth.

## Value

A smooth object of class `"fs.interaction"`, as returned by the `"fs"`
constructor.

## Details

The `"fs"` constructor reparameterises the base smooth with
`mgcv:::nat.param()` and gives each of the null-space columns its own
ridge penalty. For the usual second-order base smooth the null space is
two-dimensional (constant and linear) and the two basis vectors that
`nat.param()` returns for it span a degenerate eigenspace, so their
orientation within that span is arbitrary and can rotate under a 1e-12
relative change in `x`. Because each column carries its own penalty, the
rotation changes the fitted model, which makes criteria such as AIC
discontinuous in the data.

`"fsc"` replaces those two columns with `r * [1, (x - mean(x)) / sd(x)]`
– the same span in a fixed orientation, with `r` the root-mean-square
size of the original null-space columns – and updates the
reparameterisation matrix `P` to match, so the `"fs"` prediction method
([`mgcv::Predict.matrix()`](https://rdrr.io/pkg/mgcv/man/Predict.matrix.html)
for class `"fs.interaction"`) reproduces the fit on new data. Column
space, penalty structure and number of smoothing parameters are those of
`"fs"`. When the null space of the base smooth is not two-dimensional,
the smooth is not one-dimensional, or the canonical columns are not in
the null-space span, the result is plain `"fs"`.

## Examples

``` r
d <- data.frame(x = runif(90, 40, 80), f = factor(rep(c("a", "b", "c"), 30)))
d$y <- sin(d$x / 10) + as.numeric(d$f) + rnorm(90, 0, 0.2)
fit <- mgcv::gam(y ~ s(x, f, bs = "fsc", k = 5), data = d, method = "REML")
```
