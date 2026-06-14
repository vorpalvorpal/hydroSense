# Estimate a regularised cross-analyte correlation matrix from anchor residuals

For each analyte in `analytes`, the anchor residual column `S` (from
`d_anchors` if available and \>= 2 rows, otherwise `anchors`) is
extracted and joined into a wide `[date × analyte]` matrix. Pairwise
Pearson correlation is computed, ridge-shrunk toward the identity, and
the result is projected onto the positive-definite cone via
[`Matrix::nearPD()`](https://rdrr.io/pkg/Matrix/man/nearPD.html).

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

A named list with three elements:

- `R` – a `[p × p]` positive-definite correlation matrix (numeric
  matrix, dimnames = `list(analytes, analytes)`).

- `analytes` – the `analytes` argument, preserved for downstream use.

- `lambda` – the ridge hyperparameter (numeric scalar, currently 0.1).

## Details

### Statistical choices

**Pairwise complete observations** (`use = "pairwise.complete.obs"`):
The ARA analyte panel is ragged across sampling cadences (e.g. B.S01:
NH3-N on 121 grabs, metals on 78–79, Al on 7; complete-case would be 0
rows). Pairwise estimation uses every available co-measured date for
each pair, maximising the information extracted from the sparse overlap.

**Ridge shrinkage (fixed λ = 0.1)**: pulls low-*n* pairs toward
independence (conservative when co-observed dates are sparse). This is
more appropriate here than analytic Schäfer–Strimmer shrinkage, which
requires a complete-case matrix. The shrinkage formula is
`R_ridge = (1 - λ) * R_hat + λ * I`.

**[`Matrix::nearPD`](https://rdrr.io/pkg/Matrix/man/nearPD.html)**
(Higham 2002, alternating projection algorithm): Pairwise estimates with
unequal *n* per pair can produce indefinite matrices. nearPD projects
onto the nearest (in Frobenius norm) positive-definite correlation
matrix. Only applied when the minimum eigenvalue of `R_ridge` is \<=
1e-8; in the typical case (ridge alone suffices) the result is returned
directly without the projection step.

**Degenerate analytes** (constant `S` or \< 2 finite observations):
forced independent (identity row/col in `R_hat`) rather than propagating
NA.

## References

Higham, N. J. (2002). Computing the nearest correlation matrix — a
problem from finance. *IMA Journal of Numerical Analysis*, **22**(3),
329–343. <https://doi.org/10.1093/imanum/22.3.329>
