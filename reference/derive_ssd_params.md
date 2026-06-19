# Derive SSD parameters for msPAF computation

Reads analyte eligibility and mode-of-action groups from the bundled
metadata CSV, then calls
[`ssd_hc50()`](https://vorpalvorpal.github.io/hydroSense/reference/ssd_hc50.md)
and `.ssd_sigma()` to populate the HC50 and effective sigma needed for
Concentration Addition. Chemistry normalisation formulas (currently all
stubs) are parsed once and stored as a list column for use in
[`compute_mspaf_per_sample()`](https://vorpalvorpal.github.io/hydroSense/reference/compute_mspaf_per_sample.md).

## Usage

``` r
derive_ssd_params(meta, method, guideline_dir)
```

## Arguments

- meta:

  Analyte metadata tibble from `.load_analyte_metadata()`.

- method:

  SSD method: `"multi"` or `"anzecc"`.

- guideline_dir:

  Path to ANZG guideline data folder.

## Value

Tibble with columns `analyte`, `hc50`, `sigma`, `moa_group`,
`parsed_formula` (list of language objects or NULLs), `coanalytes_req`
(character), and `fit` (list column of fitted SSD objects, one per
analyte). The fitted SSD is loaded **once per analyte here** (via the
cached `.load_or_fit()` path) so that
[`compute_mspaf_per_sample()`](https://vorpalvorpal.github.io/hydroSense/reference/compute_mspaf_per_sample.md)
can evaluate every sample's PAF in a single vectorised
[`ssdtools::ssd_hp()`](https://bcgov.github.io/ssdtools/reference/ssd_hp.html)
call per analyte, rather than refitting / re-resolving the SSD inside a
per-(sample × analyte) loop.

## Details

Eligibility criteria:

- `ssd_available == TRUE` in the metadata

- `analyte` not in `.MSPAF_EXCLUDED_ANALYTES`

- [`ssd_hc50()`](https://vorpalvorpal.github.io/hydroSense/reference/ssd_hc50.md)
  returns a non-NA value (model fits successfully)
