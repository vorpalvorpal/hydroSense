# Assert every ammonia-bearing sample carries a water temperature

The NH3-N un-ionised-fraction normalisation requires water temperature.
Rather than silently dropping ammonia when temperature is absent, we
fail loudly: any `sample_id` that has an `NH3-N` row but no non-missing
`temperature` row triggers an error listing the offending samples.

## Usage

``` r
.assert_temperature_present(df)
```

## Arguments

- df:

  Long-format chemistry data frame (`sample_id`, `analyte`, `value`).

## Value

Invisibly `TRUE`; called for its side effect (error on violation).
