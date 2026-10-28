# Per-fit cache of the antecedent precipitation index

A tau search evaluates the API of ONE rain series at dozens of tau
values, for every analyte, then again at the chosen tau for the anchors,
the WQ layer, the static daily features and the smoother's hydro
z-score. The cache sits in front of
[`.compute_api()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-compute_api.md):
the first request for a tau evaluates the API once over every day of
`[min(dates), max(dates)]` (the dates the fit will ask for); later
requests at that tau index into it. Because the reservoir recursion runs
forward from the first hydro day, a longer evaluation grid leaves every
earlier value unchanged, so a lookup is exactly what
[`.compute_api()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-compute_api.md)
returns for those targets.

## Usage

``` r
.api_cache_new(hydro, hydro_type, dates)
```

## Arguments

- hydro:

  Data frame with `date` and `value` (as for
  [`.compute_hydro_features()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-compute_hydro_features.md)).

- hydro_type:

  `"rainfall"`, `"stage"` or `"discharge"`.

- dates:

  Every date the fit will request features for (the window).

## Value

An environment read by
[`.api_cache_features()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-api_cache_features.md).

## Details

The cache belongs to one hydro series and lives for one fit: create it
per
[`fit_target_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_target_model.md)
/
[`.fit_daily_target()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-fit_daily_target.md)
call and pass it down. Stage / discharge hydrology is not cached
([`.compute_antecedent_mean()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-compute_antecedent_mean.md)).
