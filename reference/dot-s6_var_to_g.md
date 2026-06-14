# Map an observation variance onto the g scale (delta method)

The grab measurement error (S6) is specified as a variance in
concentration/impact space; the residual smoother now works on the
`g = asinh(x / c)` scale (issue \#15), so the observation noise must be
transformed. By the delta method `Var(g) = Var(x) * g'(x)^2`, and with
`g'(x) = 1 / sqrt(x^2 + c^2)` this is `Var(x) / (x^2 + c^2)`. With
multiplicative grab error `Var(x) = (cv*x)^2` this plateaus at `cv^2`
for `|x| >> c` (proportional error) and vanishes at baseline.

## Usage

``` r
.s6_var_to_g(var_x, x, scale_c)
```

## Arguments

- var_x:

  Observation variance in concentration/impact space (numeric).

- x:

  The level at which the error is evaluated (impact `I` for the impact
  tier, concentration for the WQ tier).

- scale_c:

  Single positive transform scale `c`.

## Value

The observation variance on the `g` scale, `var_x / (x^2 + c^2)`.
