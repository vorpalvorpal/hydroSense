# Load and lightly validate the analyte metadata table

Reads the bundled CSV when `meta` is NULL. Returns a tibble. Caches the
result in an internal environment to avoid re-reading on every
[`prepare_reference()`](https://vorpalvorpal.github.io/hydroSense/reference/prepare_reference.md)
/
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)
call within a session.

## Usage

``` r
.meta_cache_env
```

## Format

An object of class `environment` of length 0.
