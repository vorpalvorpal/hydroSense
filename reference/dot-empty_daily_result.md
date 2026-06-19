# Empty tibble matching the mspaf_daily() return schema

Empty tibble matching the mspaf_daily() return schema

## Usage

``` r
.empty_daily_result(
  mode = c("point", "summary", "draws"),
  gap_uncertainty = "bracket"
)
```

## Arguments

- mode:

  One of `"point"`, `"summary"`, or `"draws"` — governs which extra
  columns are included.

- gap_uncertainty:

  Bracket mode (`"bracket"`/`"ignorable"`/ `"informative"`); governs the
  envelope columns in draws/summary mode. Ignored for point mode.
