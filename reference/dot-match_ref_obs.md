# Tier-1 contemporaneous match of target dates to reference observations

For each target date, the candidate observations are those within
`match_window_days` of it whose `hydro_short` is within
`match_hydro_tol` of the target's event API. The closest candidate in
time is chosen, ties going to the first in `obs` order
([`which.min()`](https://rdrr.io/r/base/which.min.html) semantics).

## Usage

``` r
.match_ref_obs(
  obs_date,
  obs_hydro,
  target_dates,
  target_api,
  match_window_days,
  match_hydro_tol
)
```

## Arguments

- obs_date, obs_hydro:

  Reference observation dates and stored `hydro_short`.

- target_dates:

  Distinct target dates.

- target_api:

  Event API at `target_dates`.

- match_window_days, match_hydro_tol:

  Tier-1 tolerances.

## Value

List `(matched, best)`: a logical vector (a candidate exists) and the
selected observation index (`NA` where unmatched).

## Details

The common case (no missing dates or hydrology) is resolved for all
dates at once from the observations inside each date's time window.
Inputs with missing values take the per-date path, which reproduces the
per-date `which(mask)[which.min(dt[mask])]` selection exactly, including
its behaviour when the mask itself holds `NA`.
