# Inspect detection-limit cap activations from `impute_chemistry()`

An imputed below-detection (BDL) cell must not exceed its detection
limit (DL). The posterior prediction is not itself constrained below the
DL, so
[`impute_chemistry()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_chemistry.md)
caps any imputed BDL cell whose estimate came out above the limit
(`bdl_cap = TRUE`). Frequent capping signals tension between the
modelled chemistry and the reported limits, so the cells that triggered
the cap are worth auditing rather than trusting blindly.

## Usage

``` r
bdl_cap_summary(x)
```

## Arguments

- x:

  A data frame returned by
  [`impute_chemistry()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_chemistry.md).

## Value

A tibble with one row per (`sample_id`, `analyte`) cell that exceeded
its detection limit, with columns `detection_limit`, `n_rows` (rows over
the DL — one per draw when `return = "draws"`), `max_imputed`,
`max_ratio` (`max_imputed / detection_limit`) and `capped` (whether the
cap was applied). The tibble carries the overall **`fire_rate`**
(fraction of capable BDL cells that were clipped) and **`n_bdl_cells`**
as attributes, also reported via a message. Returns `NULL` invisibly
when no cell exceeded its DL.

## Details

[`impute_chemistry()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_chemistry.md)
attaches a per-cell audit summary to its result as the
`"bdl_cap_summary"` attribute; this accessor returns it. Because plain
attributes are dropped by most dplyr verbs, call this on the frame **as
returned by
[`impute_chemistry()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_chemistry.md)**,
before further wrangling.

## See also

[`impute_chemistry()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_chemistry.md)
