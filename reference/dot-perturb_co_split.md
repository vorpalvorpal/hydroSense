# Draw coherent co-analyte perturbations for one draw (S7)

For each co-analyte grab in `fdm`, draws one lognormal multiplier
(geo-mean = 1, per-analyte or global CV from `grab_cv_co`). Each day in
`co_split` inherits the multiplier of its source grab (from
`fdm$co_grab_map`), preserving temporal coherence: consecutive forward-
filled days that came from the same grab receive the same perturbation.

## Usage

``` r
.perturb_co_split(fdm, grab_cv_co)
```

## Arguments

- fdm:

  Fitted daily scaffold from
  [`.fit_daily_target()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-fit_daily_target.md).

- grab_cv_co:

  Named numeric vector of CVs per co-analyte, or a single scalar. `NULL`
  or `NA` → return `fdm$co_split` and `fdm$wq_long` unchanged.

## Value

Named list `list(co_split, wq_long)` with perturbed values.
