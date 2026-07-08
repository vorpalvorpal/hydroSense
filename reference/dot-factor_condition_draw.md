# Draws from the conditional missing-analyte distribution, truncating BDL targets at their detection limit

Draws from the conditional missing-analyte distribution, truncating BDL
targets at their detection limit

## Usage

``` r
.factor_condition_draw(y, mu, Lambda, psi, ndraws, upper = NULL)
```

## Arguments

- y:

  Named or positional length-J vector: observed residual (or raw value,
  on whatever scale `mu`/`Lambda`/`psi` share) for observed analytes,
  `NA` for missing/BDL ones.

- mu:

  Length-J mean vector.

- Lambda:

  J x k loading matrix.

- psi:

  Length-J idiosyncratic variances (not sd).

- ndraws:

  Number of draws.

- upper:

  Optional length-`length(which(is.na(y)))` vector aligned with the
  missing block (same order
  [`.factor_condition()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-factor_condition.md)
  returns): a finite value truncates that column's draws to
  `(-inf, upper]` (a BDL target's `log(DL)`); `NA` leaves the column
  unbounded.

## Value

`ndraws` x `length(missing)` matrix.
