# Construct a standardised LMF output row

Internal helper that builds the one-row tibble returned by
[`compute_lmf_for_sample`](https://vorpalvorpal.github.io/leachatetools/reference/compute_lmf_for_sample.md)
for both successful and failed computations. All paths call this
function so the output schema is guaranteed consistent.

## Usage

``` r
make_lmf_row(
  value,
  lmf_naive,
  reason,
  n_ions,
  n_downweighted,
  sigma,
  chi2,
  quantified,
  uuid_lmf
)
```

## Arguments

- value:

  Numeric robust LMF estimate, or `NA_real_` on failure.

- lmf_naive:

  Numeric naive (non-robust) LMF estimate for comparison.

- reason:

  Character reason code, or `NA_character_` on success.

- n_ions:

  Integer count of ions used in the calculation.

- n_downweighted:

  Integer count of ions whose weight was meaningfully reduced (\> 1%
  reduction) by robust reweighting.

- sigma:

  Numeric `sigma_lmf` uncertainty estimate (from robust weights).

- chi2:

  Numeric chi-squared per degree of freedom (on original weights;
  diagnostic only).

- quantified:

  Logical; `TRUE` on success, `FALSE` on failure.

- uuid_lmf:

  UUID of the LMF analyte entry in `analyteDF`.

## Value

A one-row tibble with columns: `value` (robust LMF), `lmf_naive`,
`lmf_reason`, `n_ions_used`, `n_ions_downweighted`, `sigma_lmf`,
`chi2_per_df`, `uuid.analyte`, `uuid`, `quantified`.
