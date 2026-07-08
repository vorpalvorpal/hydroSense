# Rhat for a Route C Stage-2 Stan fit

[`.fit_group_model_factor()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-fit_group_model_factor.md)
tags its Stage-2 `CmdStanMCMC` fit with class `route_c_stanfit` so
generic diagnostics like
[`brms::rhat()`](https://mc-stan.org/posterior/reference/rhat.html)
(itself just
[`posterior::rhat()`](https://mc-stan.org/posterior/reference/rhat.html),
a plain S3 generic) work on it directly. The raw cmdstanr object has no
`rhat()` method of its own to dispatch to; this delegates to its
`$summary()` method, which already computes per-parameter Rhat correctly
(unlike calling
[`posterior::rhat()`](https://mc-stan.org/posterior/reference/rhat.html)
on the whole multi- parameter draws array, which silently blends every
parameter into one meaningless number).

## Usage

``` r
# S3 method for class 'route_c_stanfit'
rhat(x, ...)
```

## Arguments

- x:

  A `route_c_stanfit` object (a `CmdStanMCMC` with this class
  prepended).

- ...:

  Unused; present for S3 signature compatibility.

## Value

Named numeric vector of per-parameter Rhat.
