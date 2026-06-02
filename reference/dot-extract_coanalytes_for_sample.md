# Extract co-analyte values for a given sample from a long-format df

Extract co-analyte values for a given sample from a long-format df

## Usage

``` r
.extract_coanalytes_for_sample(df, sample_id, coanalytes_str)
```

## Arguments

- df:

  Long-format chemistry df

- sample_id:

  Sample identifier (may be NA for non-sample data)

- coanalytes_str:

  Comma-separated string of required co-analyte names

## Value

Named numeric vector; may be empty
