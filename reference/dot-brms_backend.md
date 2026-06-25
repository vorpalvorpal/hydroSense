# Choose the brms/Stan backend for imputation fits

Prefers **cmdstanr** when both the package and a working CmdStan install
are present, otherwise falls back to **rstan**. cmdstanr caches compiled
model binaries (a warm refit skips compilation: ~3x faster in the
package benchmark) and samples faster, while being statistically
equivalent to rstan (posterior-mean correlation 0.99999, max
standardized difference ~0.12 SD on the imputation model; see
`dev/bench-backend.R`). Override with
`options(hydroSense.brms_backend = "rstan")` (or `"cmdstanr"`). An
explicit `backend` passed through `...` to
[`fit_imputation_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_imputation_model.md)
always wins.

## Usage

``` r
.brms_backend()
```
