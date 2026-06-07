# Route a candidate analyte pool into imputation groups

Explicit-target groups claim their analytes first, in declaration order
(so an analyte listed in two groups goes to the earlier one); the single
catch-all group (`targets = NULL`) then takes whatever remains.

## Usage

``` r
.route_groups(candidate_pool, groups)
```

## Arguments

- candidate_pool:

  Character vector of modellable analyte names.

- groups:

  A validated list of
  [`impute_group()`](https://vorpalvorpal.github.io/leachatetools/reference/impute_group.md)
  objects.

## Value

A named list (one entry per group, named by `name`) of the analyte names
assigned to each group.
