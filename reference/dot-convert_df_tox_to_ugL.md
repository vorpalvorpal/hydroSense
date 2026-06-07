# Convert SSD-eligible analyte rows in a long-format data frame to µg/L

Uses `units.analyte` per-row when the column is present; otherwise
applies `conc_units` uniformly to all SSD-eligible rows. Both paths
error loudly rather than silently assuming a unit.

## Usage

``` r
.convert_df_tox_to_ugL(df, ssd_analytes, conc_units = NULL, call_arg = "df")
```

## Arguments

- df:

  Long-format data frame with `analyte` and `value` columns.

- ssd_analytes:

  Character vector of SSD-eligible analyte names.

- conc_units:

  Character unit string; required when `units.analyte` is absent from
  `df`.

- call_arg:

  Name of the data-frame argument in the caller, for errors.

## Value

`df` with SSD-eligible `value` rows converted to µg/L.
