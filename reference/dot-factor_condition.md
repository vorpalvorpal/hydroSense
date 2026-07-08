# Gaussian conditional of missing analytes given observed ones under a low-rank factor covariance

Uses the k x k Woodbury form (posterior over the latent factor `f`)
rather than inverting the full J x J residual covariance \`Sigma =
Lambda %\*% t(Lambda)

- diag(psi)\` – the whole reason for the low-rank structure. Stable even
  from a single observed analyte.

## Usage

``` r
.factor_condition(y, mu, Lambda, psi)
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

## Value

`list(mean, cov)` for the missing block, in the order `which(is.na(y))`.
