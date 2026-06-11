# Broadcast a draw-carrier frame to uniform draw coverage

After this call every row has a concrete integer `draw_id`; downstream
code can `group_by(draw_id)` with no `NA`-special-casing.

## Usage

``` r
.broadcast_draws(df, draws = .draw_domain(df))
```

## Arguments

- df:

  Long-format data frame following the draw-carrier contract.

- draws:

  Integer vector of active draw IDs (default: `.draw_domain(df)`).

## Value

`df` with `draw_id` filled for all rows; exact rows replicated N times.

## Details

Exact cells (`draw_id = NA`) are replicated once per draw in `draws`
with `draw_id` filled. Draw cells pass through unchanged, order
preserved.

**Call inside a per-site or per-sample scope** — not on the whole
dataset — because replication multiplies exact rows by N.
