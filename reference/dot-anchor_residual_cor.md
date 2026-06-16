# Estimate a regularised cross-analyte correlation matrix from anchor residuals

For each analyte in `analytes`, the anchor residual column `S` (from
`d_anchors` if available and \>= 2 rows, otherwise `anchors`) is
extracted and joined into a wide `[date × analyte]` matrix. Pairwise
Pearson correlation is computed and shrunk toward 0 in Fisher-z space by
a reliability weight driven by each pair's co-observation count, and the
result is projected onto the positive-definite cone via
[`Matrix::nearPD()`](https://rdrr.io/pkg/Matrix/man/nearPD.html) only if
still needed.

## Usage

``` r
.anchor_residual_cor(tm, analytes)
```

## Arguments

- tm:

  A target-model list with element `models`, a named list of per-analyte
  model objects. Each model must have an `anchors` data frame with
  `date` and `S` columns, and optionally a `d_anchors` data frame with
  the same structure.

- analytes:

  Character vector of analyte names (must be names of `tm$models`).

## Value

A named list with five elements:

- `R` – a `[p × p]` positive-definite correlation matrix (numeric
  matrix, no dimnames).

- `analytes` – the `analytes` argument, preserved for downstream use.

- `tau2` – numeric scalar; the estimated prior variance of the true
  Fisher-z correlation across pairs (the EB shrinkage target; 0 when no
  pair is reliable).

- `n_co` – `[p × p]` integer matrix of co-observation counts;
  `n_co[i, j]` is the number of dates with finite `S` for both analyte i
  and j (diagonal = per-analyte finite-S count).

- `nearpd_applied` – logical scalar; `TRUE` iff the
  [`Matrix::nearPD()`](https://rdrr.io/pkg/Matrix/man/nearPD.html)
  fallback fired.

## Details

### Statistical choices

**Pairwise complete observations** (`use = "pairwise.complete.obs"`):
The ARA analyte panel is ragged across sampling cadences (e.g. B.S01:
NH3-N on 121 grabs, metals on 78–79, Al on 7; complete-case would be 0
rows). Pairwise estimation uses every available co-measured date for
each pair, maximising the information extracted from the sparse overlap.

**Reliability-weighted Fisher-z shrinkage**: a fixed ridge toward the
identity damps every pair equally regardless of how well it is
estimated, so a pair with n = 3 co-observed dates and a pair with n =
100 are shrunk by the same fraction. Instead, each pairwise correlation
`r` (from `n` co-observed dates) is Fisher-z transformed,
`z = atanh(r)`, which has approximate sampling variance
`v = 1 / (n - 3)` (Fisher 1915). An empirical-Bayes (James–Stein-type)
weight `w = tau2 / (tau2 + v)` shrinks `z` toward 0: `z_shrunk = w * z`,
`r_shrunk = tanh(z_shrunk)`. Low-*n* pairs have large `v`, hence small
`w`, hence near-total shrinkage to independence; well-sampled pairs have
small `v`, `w` near 1, and survive close to their raw estimate. Pairs
with `n < 4` have undefined/negative `v` and are set to `r_shrunk = 0`
directly, killing spurious near-±1 correlations from a handful of
coincidentally aligned dates *before* the projection step below. `tau2`,
the prior variance of the true z across pairs, is estimated by
reliability-weighted method of moments over the reliable pairs
(`n >= 4`, both analytes non-degenerate):
`tau2 = max(0, weighted.mean(z^2 - v, w = n - 3))`, following the
data-driven shrinkage-intensity idea of Schäfer & Strimmer (2005) (their
covariance-shrinkage estimator is adapted here to a per-pair,
count-aware weight rather than a single global intensity, because the
ragged anchor panel gives every pair a different effective sample size).
When no pair is reliable, `tau2 = 0`, every weight is 0, and the whole
matrix collapses to the identity — the conservative, independence
assumption.

**[`Matrix::nearPD`](https://rdrr.io/pkg/Matrix/man/nearPD.html)**
(Higham 2002, alternating projection algorithm): shrinking toward 0 in
z-space already removes the spurious low-*n* off-diagonals that
previously made the pairwise matrix indefinite, so this safety net
should rarely fire. It is retained as a fallback: applied only when the
minimum eigenvalue of the shrunk matrix is \<= 1e-8, in which case it
projects onto the nearest (Frobenius norm) positive-definite correlation
matrix; otherwise the shrunk matrix is used directly.

**Degenerate analytes** (constant `S` or \< 2 finite observations):
forced independent (identity row/col in `R_hat`) rather than propagating
NA.

## References

Fisher, R. A. (1915). Frequency distribution of the values of the
correlation coefficient in samples from an indefinitely large
population. *Biometrika*, **10**(4), 507–521.
<https://doi.org/10.2307/2331838>

Schäfer, J., & Strimmer, K. (2005). A shrinkage approach to large-scale
covariance matrix estimation and implications for functional genomics.
*Statistical Applications in Genetics and Molecular Biology*, **4**(1),
Article 32. <https://doi.org/10.2202/1544-6115.1175>

Higham, N. J. (2002). Computing the nearest correlation matrix — a
problem from finance. *IMA Journal of Numerical Analysis*, **22**(3),
329–343. <https://doi.org/10.1093/imanum/22.3.329>
