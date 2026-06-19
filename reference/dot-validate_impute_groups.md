# Validate and normalise a list of imputation groups

Checks that `groups` is a non-empty list of
[`impute_group()`](https://vorpalvorpal.github.io/hydroSense/reference/impute_group.md)
objects with unique names and at most one catch-all (`targets = NULL`)
entry.

## Usage

``` r
.validate_impute_groups(groups)
```

## Arguments

- groups:

  A list of `impute_group` objects.

## Value

The validated `groups` list (invisibly the same object).
