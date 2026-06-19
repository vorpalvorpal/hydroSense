# Perturb all GAM fits inside a target_model for one draw

Draws posterior-coefficient samples for the impact GAMs and WQ-layer
GAMs. Pooled analytes (all sharing one factor-smooth `impact_fit`)
receive **one shared draw** to preserve cross-analyte coherence — the
same perturbed fit object is assigned to every pooled analyte.
Non-pooled `impact_fit`s and all `wq_fit`s are perturbed independently.

## Usage

``` r
.perturb_target_model(tm, perturb_reference = FALSE)
```

## Arguments

- tm:

  A `target_model` from
  [`fit_target_model()`](https://vorpalvorpal.github.io/hydroSense/reference/fit_target_model.md).

- perturb_reference:

  Logical; when `TRUE`, also perturb the embedded `reference_model`'s
  GAMs (use only when `mspaf_daily(reference = NULL)` —
  total-concentration mode — so that `ref_norm` uncertainty is not
  inadvertently cancelled by a downstream subtraction).

## Value

A modified copy of `tm`. The original object is not mutated (R
copy-on-modify semantics).
