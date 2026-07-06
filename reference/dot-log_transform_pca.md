# Log10-transform the concentration-like columns of a PCA matrix

Applies [`log10()`](https://rdrr.io/r/base/Log.html) to every column
whose name is not in
[.PCA_NO_LOG_VARS](https://vorpalvorpal.github.io/hydroSense/reference/dot-PCA_NO_LOG_VARS.md),
leaving pH / temperature / ORP / DO on their native scale. `NA` cells
are preserved (so NIPALS can still handle within-sample missingness). A
genuine zero is floored at half the column's own smallest observed
positive value
([`.scale_aware_log_floor()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-scale_aware_log_floor.md))
rather than one absolute constant shared across analytes of wildly
different scale. Both the training PCA
([`.prepare_chem_pca()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-prepare_chem_pca.md))
and the scoring projection
([`.compute_pca_scores()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-compute_pca_scores.md))
call this so the transform is identical on both paths; scoring passes
the training-derived `floors` so a cell transforms identically at fit
and predict regardless of how many rows the scoring call sees.

## Usage

``` r
.log_transform_pca(
  mat,
  no_log_vars = .PCA_NO_LOG_VARS,
  eps = 1e-09,
  floors = NULL
)
```

## Arguments

- mat:

  Numeric matrix with named columns (samples × variables).

- no_log_vars:

  Column names to leave on their native scale. Default
  [.PCA_NO_LOG_VARS](https://vorpalvorpal.github.io/hydroSense/reference/dot-PCA_NO_LOG_VARS.md).

- eps:

  Absolute fallback floor for columns with no positive values (default
  `1e-9`).

- floors:

  Optional named numeric vector of per-column floors (as produced by a
  prior call's `attr(., "floors")`). When `NULL`, floors are computed
  from `mat` itself. Supply the training floors at scoring time.
