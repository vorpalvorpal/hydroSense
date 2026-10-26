# Resolve ref_norm for instantaneous target samples (one date per sample)

Vectorised over target dates: the antecedent precipitation index is
evaluated in one pass over all distinct dates (one event-API pass for
the tier-1 gate, one short- and one long-store pass per tier-2 analyte
model), and each tier-2 reference GAM is predicted once on a frame of
every date left unmatched by tier 1. Results are identical to resolving
each date on its own.

## Usage

``` r
.resolve_ref_norm_instant(ref_model, target_df)
```

## Arguments

- ref_model:

  A `reference_model` object.

- target_df:

  Tibble with columns `sample_id`, `datetime` (Date).

## Value

Tibble `(sample_id, analyte, ref_norm, ref_tier)`.

## Details

Tiers, per (sample, analyte):

1.  `"matched"`: a reference observation within `match_window_days`
    whose stored `hydro_short` lies within `match_hydro_tol` of the
    target's event API (computed at `tau = match_window_days`); the
    closest in time wins, ties going to the earlier observation in
    `m$obs`.

2.  `"model"`: the reference GAM's prediction at the date.

3.  `"static"`: the static reference value (also the fallback when the
    prediction is unavailable).
