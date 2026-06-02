# Log10-transform the concentration-like columns of a PCA matrix

Applies [`log10()`](https://rdrr.io/r/base/Log.html) to every column
whose name is not in
[.PCA_NO_LOG_VARS](https://www.kedumba.com.au/leachatetools/reference/dot-PCA_NO_LOG_VARS.md),
leaving pH / temperature / ORP / DO on their native scale. `NA` cells
are preserved (so NIPALS can still handle within-sample missingness),
and a small floor guards against `log10(0)` for genuine zeros. Both the
training PCA
([`.prepare_chem_pca()`](https://www.kedumba.com.au/leachatetools/reference/dot-prepare_chem_pca.md))
and the scoring projection
([`.compute_pca_scores()`](https://www.kedumba.com.au/leachatetools/reference/dot-compute_pca_scores.md))
call this so the transform is identical on both paths.

## Usage

``` r
.log_transform_pca(mat, eps = 1e-09)
```

## Arguments

- mat:

  Numeric matrix with named columns (samples × variables).

- eps:

  Positive floor applied before the log (default `1e-9`).
