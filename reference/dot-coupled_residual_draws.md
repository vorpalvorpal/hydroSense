# Draw cross-analyte correlated Kalman residual trajectories

Generates `ndraws` draws from each analyte's Durbin–Koopman simulation
smoother, with innovations correlated **across analytes** according to
the `[p × p]` correlation sub-matrix `cor_R`. Each analyte's per-analyte
marginal distribution is unchanged by the correlation (DK identity);
only the *joint* distribution across analytes is affected.

## Usage

``` r
.coupled_residual_draws(
  smoothers,
  modelled,
  ndraws,
  cor_R,
  cor_analytes = rownames(cor_R),
  seed = NULL
)
```

## Arguments

- smoothers:

  Named list by analyte. Each element is a list with `grid_dates` (Date
  vector), `mean` (numeric vector, KFS posterior mean), and `draw_model`
  (KFAS model from
  [`.build_kalman_model()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-build_kalman_model.md),
  or `NULL` for a degenerate / non-modelled analyte).

- modelled:

  Character vector of analyte names to include in the output (subset of
  `names(smoothers)`).

- ndraws:

  Positive integer; number of draws to generate.

- cor_R:

  `[p × p]` numeric positive-definite correlation matrix (no dimnames
  required). Typically the `$R` component from
  [`.anchor_residual_cor()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-anchor_residual_cor.md).

- cor_analytes:

  Character vector giving the analyte order for the rows and columns of
  `cor_R`. Defaults to `rownames(cor_R)` when `cor_R` has dimnames; must
  be supplied explicitly otherwise.

- seed:

  Integer or `NULL`. When non-`NULL`, `set.seed(seed)` is called at the
  start of the function so results are reproducible. In
  [`amspaf_daily()`](https://vorpalvorpal.github.io/leachatetools/reference/amspaf_daily.md),
  the seed is set by the caller (via its own `seed` argument) before
  this function is invoked, so `seed = NULL` is used there.

## Value

Named list by analyte (same order as `modelled`). Each element is
`list(grid_dates = <Date>, draws = <matrix [n_grid × ndraws]>)`.
Analytes with `draw_model = NULL` receive flat draws equal to their KFS
posterior mean (the degenerate path). Analytes with empty `grid_dates`
receive `draws = NULL`.

## Details

### Statistical justification

**DK identity** (Durbin & Koopman 2002, §4): the draw \$\$\tilde{\alpha}
= \hat{\alpha} + (\alpha^+ - L y^+)\$\$ guarantees each analyte's
marginal distribution is unchanged when the process innovations
\\\eta^+\\ are correlated: the smoother gain `L` is precomputed on the
actual data, and the correction \\\alpha^+ - L y^+\\ has expectation
zero. Correlating \\\eta^+\\ across analytes widens only the *joint*
distribution (i.e. the combined AmsPAF interval).

**Cholesky coupling**: `Z_raw %*% U_chol` (where `U_chol = chol(R)`,
upper triangular) gives rows with covariance `R` when rows of `Z_raw`
are iid \\N(0, I_p)\\. The same `U_chol` is applied to both the process
innovations and the initial-state draws so that the coupling is coherent
along the entire trajectory.

**Independent observation noise**: `eps_std` is drawn independently per
analyte because measurement errors from separate grab analyses are
uncorrelated even when the true concentrations co-move.

**Ragged grids**: `match(grid_a, U)` slices the union-grid field to each
analyte's own clipped grid. Days where an analyte has no entry (outside
its grab span) have no field element — no coupling on those days, which
is correct.

## References

Durbin J, Koopman SJ (2002). A simple and efficient simulation smoother
for state space time series analysis. *Biometrika* **89**(3), 603–616.
<https://doi.org/10.1093/biomet/89.3.603>
