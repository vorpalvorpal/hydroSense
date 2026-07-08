# Extract per-draw (Lambda, Psi) pairs from a fitted Route C group

Extract per-draw (Lambda, Psi) pairs from a fitted Route C group

## Usage

``` r
.route_c_draw_params(group, ndraws = NULL, return = "draws")
```

## Arguments

- group:

  A fitted `factor`-method group (from
  [`.fit_group_model_factor()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-fit_group_model_factor.md)).

- ndraws:

  Optional cap on the number of posterior draws to use (the first
  `ndraws`); `NULL` uses every draw.

## Value

`list(Lambda = list of J x k matrices, Psi = list of length-J vectors)`,
one entry per used draw, plus `n` (the number used).
