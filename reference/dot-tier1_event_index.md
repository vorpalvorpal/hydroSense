# Tier-1 event index: the antecedent rain index behind the matching gate

The single definition of the hydrological index that tier-1 direct
matching compares. It is used for all three quantities in the gate – the
index stored on each reference observation (`obs$hydro_event`), the
target date's index, and the IQR-based default tolerance
`match_hydro_tol` – so that the comparison and its tolerance are on one
scale.
[`.compute_api()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-compute_api.md)
is an un-normalised store whose level grows roughly with `tau`, so
mixing taus (as the gate once did: a 30-d store against a 5-d one, gated
by the spread of a 1-d one) leaves same-day, same-rain pairs tens of mm
apart and tier 1 practically never fires.

## Usage

``` r
.tier1_event_index(hydro_values, hydro_dates, target_dates, match_window_days)
```

## Arguments

- hydro_values, hydro_dates:

  Daily hydrology series.

- target_dates:

  Date vector; the index is evaluated at each.

- match_window_days:

  Tier-1 time window (days).

## Value

Numeric vector, same length as `target_dates`.

## Details

The memory constant is `match_window_days`: the gate asks whether the
rain over roughly the matching window differs between the two dates,
which is the time scale on which a grab within that window is a valid
proxy. It is floored at 1 day so `match_window_days = 0` (same-day
matches only) still yields a valid reservoir.
