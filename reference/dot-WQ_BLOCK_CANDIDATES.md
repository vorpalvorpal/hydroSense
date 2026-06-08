# WQ analytes used in the PCA pre-processing step

All candidates; the actual variables used at any given site are the
intersection of this list with analytes that pass prescreen in the
training data. ORP and DO are included because they are valuable at
sites where they are measured, but at most sites they will be filtered
out by prescreen.

## Usage

``` r
.WQ_BLOCK_CANDIDATES
```

## Format

An object of class `character` of length 34.

## Details

Note: the sulfate entry encodes its superscript charge (the SO4
two-minus symbol) with Unicode escape sequences in the source string, so
that matching works regardless of how the source file is parsed.
