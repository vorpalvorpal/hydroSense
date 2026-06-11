# Expand observed measurements into lognormal posterior-predictive draws

Turns observed, detected, exact (`draw_id = NA`, `!imputed`) cells that
carry a per-measurement uncertainty column into N lognormal draws,
propagating analytical error through the pipeline alongside imputation
uncertainty from
[`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md)
/
[`impute_coanalytes()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_coanalytes.md).

## Usage

``` r
draw_measurement_error(
  df,
  error_col,
  error_type = c("cv", "sd"),
  ndraws = NULL,
  seed = NULL
)
```

## Arguments

- df:

  Long-format data frame following the draw-carrier contract.

- error_col:

  Name (string) of the column carrying per-measurement uncertainty. Must
  exist in `df`. Rows with NA or zero values are not expanded.

- error_type:

  `"cv"` (default): `error_col` contains coefficients of variation (e.g.
  `0.05` = 5 % RSD, matching the
  [`add_lmf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_lmf.md)
  convention). `"sd"`: `error_col` contains absolute standard deviations
  on the natural concentration scale; CV is derived as `sd / value`.

- ndraws:

  Number of draws. Required when `df` carries no existing draws. When
  draws are already present N is inferred from the draw domain; `ndraws`
  must match or be `NULL`.

- seed:

  Optional integer seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before sampling,
  for reproducibility.

## Value

`df` with eligible observed cells replaced by N lognormal draw rows
(`draw_id 1..N`). Ineligible cells are unchanged. Returns `df`
unmodified when no eligible cells exist (including when all `error_col`
values are `NA`).

## Details

Cells without reported error (NA or zero in `error_col`), BDL cells
(`detected = FALSE`), and already-drawn or imputed cells pass through
unchanged.

**Lognormal parameterisation — geometric-mean convention:** the reported
value \\v\\ is the geometric mean (median) of the draw distribution. For
coefficient of variation \\CV\\:

\$\$\sigma\_{\log} = \sqrt{\log(1 + CV^2)},\quad \mu = \log v,\quad
\text{draws} = \exp(\mathcal{N}(\mu,\\ \sigma\_{\log}^2))\$\$

The arithmetic mean of the draws is \\v \exp(\sigma\_{\log}^2/2)\\,
which exceeds \\v\\ by \< 0.5 % at 10 % RSD.

**Documented independence approximations:** measurement error is treated
as independent of imputation error (cross-source index-pairing) and
independent across analytes within a sample.

**Call order:**
[`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md)
→
[`impute_coanalytes()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_coanalytes.md)
→ `draw_measurement_error()` →
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md).
Imputation establishes the draw count N; this function reuses it.

## See also

[`impute_chemistry()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_chemistry.md),
[`impute_coanalytes()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_coanalytes.md),
[`summarise_draws()`](https://vorpalvorpal.github.io/leachatetools/reference/summarise_draws.md)
