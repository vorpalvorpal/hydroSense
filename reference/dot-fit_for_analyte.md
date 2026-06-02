# Internal: load data and fit an SSD for one analyte.

Called by `.load_or_fit()` in paf.R on a cache miss.

## Usage

``` r
.fit_for_analyte(analyte, stem, meta, dists, guideline_dir)
```

## Arguments

- analyte:

  Character. Canonical analyte name.

- stem:

  Character. Safe file stem (from .SSD_NAME_MAP).

- meta:

  One-row data.frame from anzecc_analyte_metadata.csv.

- dists:

  Character vector. Distribution names for ssd_fit_dists().

- guideline_dir:

  Character. Path to the "guideline data/" folder. Only required for
  ANZG_XLSX analytes; may be NULL for Warne2000 analytes.

## Value

A fitted ssdtools object with provenance attributes, or NULL on error.
