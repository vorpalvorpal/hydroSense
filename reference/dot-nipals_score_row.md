# NIPALS regression scoring for one centred/scaled observation

Computes PC scores for a single row by projecting onto each loading
vector using only observed (non-NA) elements, then deflating the
residual before moving to the next component. For a fully observed row
this is identical to the standard `x %*% loadings` projection. For a row
with missing values it correctly down-weights the loading vectors to the
observed subspace — without the bias that zero/median imputation
introduces.

## Usage

``` r
.nipals_score_row(x, loadings, n_pcs)
```

## Arguments

- x:

  Numeric vector (centred + scaled); `NA` marks missing variables.

- loadings:

  p × K loading matrix (unit-normalised columns from nipals).

- n_pcs:

  Number of components to score.
