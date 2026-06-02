# Stop with a friendly, actionable message if brms is not installed

The Bayesian imputation step
([`fit_imputation_model()`](https://www.kedumba.com.au/leachatetools/reference/fit_imputation_model.md)
/
[`impute_chemistry()`](https://www.kedumba.com.au/leachatetools/reference/impute_chemistry.md))
is the only part of the package that needs brms, so brms is an optional
("Suggests") dependency rather than a hard requirement. This keeps the
package quick to install for users who only need the LMF or AmsPAF
tools. When someone actually calls an imputation function without brms
installed, this guard explains — in plain language — what to install and
why.

## Usage

``` r
.require_brms()
```
