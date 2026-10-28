# API at `target_dates` for one tau, through the per-fit cache

Targets outside the cache window (or a cache without one) are evaluated
directly by
[`.compute_api()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-compute_api.md)
and not stored.

## Usage

``` r
.api_cache_get(cache, target_dates, tau)
```
