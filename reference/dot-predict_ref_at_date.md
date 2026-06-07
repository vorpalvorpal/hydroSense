# Predict reference value_norm at a single date for one analyte model

Predict reference value_norm at a single date for one analyte model

## Usage

``` r
.predict_ref_at_date(m, target_date, hydro, hydro_type, eps = 1e-09)
```

## Arguments

- m:

  Analyte model object from `reference_model$models[[nm]]`.

- target_date:

  Date scalar.

- hydro:

  Daily hydro series.

- hydro_type:

  Hydro type string.

- eps:

  Log guard.

## Value

List `(ref_norm, ref_tier)`.
