# Retrieve per-cell ARA diagnostics from an `add_mspaf()` result

After calling
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)
with a `reference_model` (or any reference), this accessor returns a
tibble describing what happened in the ARA subtraction for every (sample
× analyte) that was assessed. This is the primary tool for auditing the
"reference higher than target" case (floored to zero) and for
understanding which tier was used per cell.

## Usage

``` r
ara_summary(x)
```

## Arguments

- x:

  A data frame returned by
  [`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md).

## Value

A tibble with columns:

- `sample_id`:

  Sample identifier.

- `analyte`:

  Analyte name.

- `ref_norm`:

  Normalised reference concentration subtracted.

- `C_norm`:

  Normalised target concentration (before ARA).

- `C_adj`:

  ARA-adjusted concentration (`max(C_norm - ref_norm, 0)`).

- `C_excess`:

  Unfloored difference `C_norm - ref_norm`; negative values indicate the
  reference exceeded the target — possibly a geogenic artefact (e.g.
  low-pH upstream metal mobilisation).

- `floor_fired`:

  Logical; `TRUE` when `C_norm < ref_norm`.

- `ref_source`:

  `"disabled"`, `"matched"`, or `"unmatched"`.

- `ref_tier`:

  `"matched"`, `"model"`, `"model_integrated"`, or `"static"`. `NA` for
  non-temporal reference.

Returns `NULL` (with a message) if the attribute is absent.

## Details

The attribute is stored by
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)
and is dropped by most dplyr verbs, so read the summary before further
wrangling.

## See also

[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md),
[`fit_reference_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_reference_model.md)
