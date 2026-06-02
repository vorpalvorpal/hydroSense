# Vectorised SSD PAF lookup for one analyte

Returns the proportion of species affected at each concentration in
`conc`. Concentrations that are `NA` or `<= 0` map to `PAF = 0`. When a
fitted SSD object is available it is evaluated in a single
[`ssdtools::ssd_hp()`](https://bcgov.github.io/ssdtools/reference/ssd_hp.html)
call over the positive concentrations; otherwise it falls back to a
per-value
[`ssd_paf()`](https://www.kedumba.com.au/leachatetools/reference/ssd_paf.md)
lookup (which re-resolves the model internally).

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
  [`ssd_paf()`](https://www.kedumba.com.au/leachatetools/reference/ssd_paf.md)
  lookup only).

- method:

  SSD method (fallback only).

- guideline_dir:

  Path to ANZG XLSX folder (fallback only).

## Value

Numeric vector the same length as `conc`, PAF as a proportion (0–1).
