# Stop with a friendly message if cmdstanr is unavailable

The Route C `"factor"` imputation method uses cmdstanr + mgcv (it does
**not** use brms), so this is the method-appropriate engine guard – the
analogue of
[`.require_brms()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-require_brms.md)
for the three brms-based methods.

## Usage

``` r
.require_cmdstanr()
```
