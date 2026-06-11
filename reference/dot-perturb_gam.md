# Perturb a single mgcv GAM by drawing from its posterior coefficient distribution

Draws one set of coefficients from `N(coef(g), Vp)` and assigns them to
`g$coefficients` so that subsequent `predict.gam()` calls use the draw.
`predict.gam()` dispatches through `g$coefficients`, so this gives an
honest posterior-predictive draw without refitting.

## Usage

``` r
.perturb_gam(g)
```

## Arguments

- g:

  A fitted [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) object,
  or `NULL`.

## Value

`g` with `$coefficients` replaced by a single posterior draw, or `g`
unchanged on failure.

## Details

Returns `g` unperturbed when `$Vp` is absent, NULL, or when the draw
produces non-finite values (near-singular posterior). A symmetric PD
jitter is attempted once before giving up.
