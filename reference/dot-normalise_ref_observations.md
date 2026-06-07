# Normalise reference observations to the SSD index condition

Applies per-analyte bioavailability / physicochemical normalisation
formulas (hardness, pH, DOC) to a long-format chemistry data frame,
returning the frame with a `value_norm` column added. BDL rows
(`detected == FALSE`) are normalised at their detection-limit value but
callers are responsible for setting `value = 0` for BDL rows *before*
calling this if they want BDL to contribute zero to downstream summaries
(as
[`prepare_reference()`](https://vorpalvorpal.github.io/leachatetools/reference/prepare_reference.md)
does).

## Usage

``` r
.normalise_ref_observations(ref_df, original_data, meta)
```

## Arguments

- ref_df:

  Long-format chemistry data frame with at minimum `analyte`, `value`,
  `detected`, and `sample_id`.

- original_data:

  The original reference chemistry frame used to look up co-analyte
  values (pH, DOC, hardness, Ca, Mg) per sample.

- meta:

  Analyte metadata tibble from `.load_analyte_metadata()`.

## Value

`ref_df` with a `value_norm` numeric column appended (NA where
normalisation fails due to a missing required co-analyte).
