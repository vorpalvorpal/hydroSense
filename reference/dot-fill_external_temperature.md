# Fill temperature from an external daily series on non-grab days

Augments `daily_long` with temperature values from `temperature_df` for
dates where there is no directly measured grab temperature. Dates with a
`.measured == TRUE` temperature row (actual grab measurement) are kept
unchanged.

## Usage

``` r
.fill_external_temperature(daily_long, temperature_df)
```

## Arguments

- daily_long:

  Output of
  [`.build_daily_chem()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-build_daily_chem.md).

- temperature_df:

  Data frame with `datetime` and `value` columns.

## Value

`daily_long` with temperature rows augmented.
