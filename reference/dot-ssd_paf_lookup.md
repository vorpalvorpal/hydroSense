# Resolve or build a spline-based PAF lookup closure for one analyte

Returns a [`stats::splinefun`](https://rdrr.io/r/stats/splinefun.html)
closure that maps `log10(conc)` to PAF. Shipped tables
(`NULL guideline_dir` + known method/analyte) are loaded from
`inst/extdata/ssd_paf_lookup.qs2` and cached in the session. Runtime
tables are built adaptively and likewise cached.

## Usage

``` r
.ssd_paf_lookup(analyte, method, fit, guideline_dir)
```

## Arguments

- analyte:

  Analyte name string.

- method:

  SSD method string (e.g. `"multi"`, `"anzecc"`).

- fit:

  `fitdists` object or `NULL`. Required for runtime builds.

- guideline_dir:

  Path to ANZG XLSX folder, or `NULL` for shipped tables.

## Value

A `splinefun` closure, or `NULL` if the lookup cannot be built.
