# Vectorised SSD PAF lookup for one analyte

Returns the proportion of species affected at each concentration in
`conc`. Concentrations that are `NA`, `<= 0`, or non-finite map to
`PAF = 0`. When a fitted SSD object is available it is evaluated via a
spline lookup (fast path) or a single
[`ssdtools::ssd_hp()`](https://bcgov.github.io/ssdtools/reference/ssd_hp.html)
call (exact fallback); otherwise it falls back to a per-value
[`ssd_paf()`](https://vorpalvorpal.github.io/leachatetools/reference/ssd_paf.md)
lookup.

## Usage

``` r
.ssd_paf_vec(fit, conc, analyte, method, guideline_dir)
```

## Arguments

- fit:

  Fitted SSD object (from the `fit` list-column) or `NULL`.

- conc:

  Numeric vector of adjusted concentrations.

- analyte:

  Analyte name (fallback
  [`ssd_paf()`](https://vorpalvorpal.github.io/leachatetools/reference/ssd_paf.md)
  lookup only).

- method:

  SSD method (fallback only).

- guideline_dir:

  Path to ANZG XLSX folder (fallback only).

## Value

Numeric vector the same length as `conc`, PAF as a proportion (0–1).
