# Add per-analyte PAF to adjusted tox rows, batched per analyte

Evaluates each analyte's SSD once across **all** its `C_adj` values in a
single vectorised
[`ssdtools::ssd_hp()`](https://bcgov.github.io/ssdtools/reference/ssd_hp.html)
call (via
[`.ssd_paf_vec()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-ssd_paf_vec.md)),
rather than one
[`ssd_paf()`](https://vorpalvorpal.github.io/leachatetools/reference/ssd_paf.md)
call per row. The fitted SSD object is taken from the `fit` list-column
of `ssd_params` (loaded once in
[`derive_ssd_params()`](https://vorpalvorpal.github.io/leachatetools/reference/derive_ssd_params.md)).

## Usage

``` r
.amspaf_add_paf(tox, ssd_params, method, guideline_dir)
```

## Arguments

- tox:

  Adjusted tox rows from
  [`.amspaf_adjust()`](https://vorpalvorpal.github.io/leachatetools/reference/dot-amspaf_adjust.md)
  (needs `analyte`, `C_adj`); may span multiple samples.

- ssd_params:

  Tibble from
  [`derive_ssd_params()`](https://vorpalvorpal.github.io/leachatetools/reference/derive_ssd_params.md)
  (uses `analyte`, `fit`).

- method:

  SSD method (NULL-fit fallback only).

- guideline_dir:

  Path to ANZG XLSX folder (NULL-fit fallback only).

## Value

`tox` with a numeric `PAF` column added (proportion, 0–1).
