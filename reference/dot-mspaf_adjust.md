# Normalise chemistry and apply the ARA shift for a block of sample rows

Filters to SSD-eligible analytes, joins SSD params + the reference
table, applies the per-analyte chemistry normalisation (BDL → 0),
records the `ref_source` diagnostic, and computes the
added-risk-adjusted concentration `C_adj = max(C_norm - ref_norm, 0)`.
Operates on one *or many* samples; the caller is responsible for any
per-sample grouping downstream.

## Usage

``` r
.mspaf_adjust(sample_rows, ref_table, ssd_params, ara_enabled = TRUE)
```

## Arguments

- sample_rows:

  Long-format chemistry rows (one or more samples). Must carry
  `detected`.

- ref_table:

  Tibble `(analyte, ref_norm)`.

- ssd_params:

  Tibble from
  [`derive_ssd_params()`](https://vorpalvorpal.github.io/hydroSense/reference/derive_ssd_params.md).

- ara_enabled:

  Logical; whether the caller supplied a reference.

## Value

A list `list(tox, dropped)` where `tox` has the eligible,
normalisation-resolved rows (with `C_norm`, `ref_norm`, `ref_source`,
`C_adj`) and `dropped` is a `(analyte, reason)` tibble of rows dropped
for a missing required co-analyte.

## Details

`ref_source` distinguishes:

- `"disabled"` — ARA off (no reference supplied to
  [`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md));

- `"matched"` — ARA on and a reference value was found for the analyte;

- `"unmatched"` — ARA on but no reference match (assessed against raw
  normalised concentration, i.e. `ref_norm = 0`).
