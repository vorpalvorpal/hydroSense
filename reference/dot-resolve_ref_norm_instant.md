# Resolve ref_norm for instantaneous target samples (one date per sample)

Vectorised over target dates: the antecedent precipitation index is
evaluated in one pass over all distinct dates (one event-API pass for
the tier-1 gate, one short- and one long-store pass per tier-2 analyte
model), and each tier-2 reference GAM is predicted once on a frame of
every date left unmatched by tier 1 (a failing batch is bisected, so
only dates whose own prediction errors fall back to static). Results are
identical to resolving each date on its own, including the error
fallback; an NA date resolves to static (the per-date code errored on
it).

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
    whose stored event index (`hydro_event`) lies within
    `match_hydro_tol` of the target's; both are
    [`.tier1_event_index()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-tier1_event_index.md)
    (the API at `tau = match_window_days`), the scale `match_hydro_tol`
    is fitted on. The closest in time wins, ties going to the earlier
    observation in `m$obs`.

2.  `"model"`: the reference GAM's prediction at the date.

3.  `"static"`: the static reference value (also the fallback when the
    prediction is unavailable).
