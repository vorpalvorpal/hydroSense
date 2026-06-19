# Build per-draw co-analyte rows for the synthetic frame (S7, G2)

Creates a copy of the exact co-analyte rows from `daily_long_exact` with
values replaced by the corresponding perturbed values from
`co_pert_split` (output of
[`.perturb_co_split()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-perturb_co_split.md)),
and tags the result with `draw_id`. Analytes absent from `co_pert_split`
(e.g. non-modelled toxicants) are left unchanged. The resulting rows go
into the synthetic frame so that
[`add_mspaf()`](https://vorpalvorpal.github.io/hydroSense/reference/add_mspaf.md)'s
normalisation sees the per-draw co-analyte perturbation.

## Usage

``` r
.co_draw_rows(daily_long_exact, co_pert_split, draw_id)
```

## Arguments

- daily_long_exact:

  Rows of the forward-filled daily chemistry that belong to the
  non-modelled subset (co-analytes + unmodelled toxicants).

- co_pert_split:

  Named list (date-string → data.frame(analyte, value)) from
  [`.perturb_co_split()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-perturb_co_split.md).

- draw_id:

  Integer draw index to assign to the output rows.

## Value

A copy of `daily_long_exact` with perturbed values and `draw_id` set to
`draw_id`.
