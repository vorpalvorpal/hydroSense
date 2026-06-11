# Summarise analytes assessed without a reference match (ARA coverage)

When ARA is enabled but some analytes had no matching reference value,
those analytes were assessed against their raw normalised concentration
(`ref_norm = 0`). This emits a single end-of-call cli message tallying
how many samples each such analyte affected, so the caller knows where
the ARA adjustment did *not* apply. No-op when ARA is disabled.

## Usage

``` r
.summarise_ara_coverage(pafs_long, ara_enabled)
```

## Arguments

- pafs_long:

  Flat per-analyte PAF breakdown (the `analyte_pafs` attribute).

- ara_enabled:

  Logical; whether the caller supplied a reference.
