# Active draw domain from a draw-carrier frame

Returns the sorted integer vector of draw IDs present in `df`, or
`integer(0)` for a deterministic (all-exact or legacy) frame.

## Usage

``` r
.draw_domain(df)
```

## Arguments

- df:

  Long-format data frame possibly carrying a `draw_id` column.

## Value

Sorted integer vector `1..N`, or `integer(0)`.

## Details

Asserts the domain is contiguous `1..N` and that every draw-bearing
`(sample_id, analyte)` cell carries the full set (no ragged N).
