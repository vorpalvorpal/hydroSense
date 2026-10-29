# Canonicalise, unit-convert and reduce a long chemistry frame

One value per (site, analyte, date): when any row that day is detected,
the mean of the detected values (`detected = TRUE`); otherwise the
smallest reported limit (`detected = FALSE`). Toxicants are converted to
ug/L.

## Usage

``` r
.ah_reduce_input(df, ssd_analytes, conc_units, arg)
```

## Value

`list(red, dates)`: the reduced data frame (`site_id`, `analyte`, `date`
(numeric day), `value`, `detected`) and the numeric dates of every input
row (for the grid).
