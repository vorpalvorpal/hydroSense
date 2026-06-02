# Synthetic leachate-impacted water-quality data (demo)

Returns a small, fictional long-format water-quality dataset used in the
examples and vignettes. It is **not** real monitoring data —
concentrations are generated from a simple two-component mixing model so
the worked examples are reproducible and self-contained.

## Usage

``` r
leachate_demo()
```

## Value

A [tibble](https://tibble.tidyverse.org/reference/tibble.html) with one
row per sample x analyte and the columns:

- sample_id:

  Character. Unique sample identifier.

- site_id:

  Character. One of `"downstream"`, `"reference"`, `"leachate"`.

- datetime:

  Date. Sampling date.

- analyte:

  Character. Analyte name.

- value:

  Numeric. Concentration in `units.analyte`.

- detected:

  Logical. `FALSE` marks a below-detection-limit result.

- units.analyte:

  Character. Concentration units: major ions in `"mg/L"`;
  `Cu`/`Zn`/`NH3-N` in `"ug/L"` (the SSD scale); `pH` unitless;
  `temperature` in `"degC"`; `hardness` in `"mg/L CaCO3"`.

- valence.analyte:

  Numeric. Ionic charge (for the meq conversion in
  [`to_meq()`](https://vorpalvorpal.github.io/leachatetools/reference/to_meq.md));
  `NA` for non-ionic analytes.

- atomic_mass.analyte:

  Numeric. Molar/atomic mass in g/mol (for the meq conversion); `NA` for
  non-ionic analytes.

## Details

The data are bundled in the qs2 *qdata* format under `inst/extdata` (qs2
is already a dependency) and read on demand by this accessor, rather
than being lazy-loaded R data.

Three sites are included:

- `downstream`:

  A leachate-impacted site, ~15% leachate mixing fraction relative to
  the reference.

- `reference`:

  Clean background / upstream chemistry.

- `leachate`:

  The leachate end-member.

Each site has six bi-monthly samples through 2024. The analyte panel
carries the major ions used by
[`add_lmf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_lmf.md)
/
[`to_meq()`](https://vorpalvorpal.github.io/leachatetools/reference/to_meq.md)
(`Na`, `K`, `Ca`, `Mg`, `Cl`, sulfate, `F`, `NO3-N`, `NO2-N`,
`CO3-CaCO3`, `HCO3-CaCO3`), the toxicants assessed by
[`ssd_paf()`](https://vorpalvorpal.github.io/leachatetools/reference/ssd_paf.md)
/
[`add_amspaf()`](https://vorpalvorpal.github.io/leachatetools/reference/add_amspaf.md)
(`Cu`, `Zn`, `NH3-N`), and the co-analytes needed for the ammonia and
bioavailability normalisations (`pH`, `temperature`, `DOC`, `hardness`).

## Examples

``` r
demo <- leachate_demo()
table(demo$site_id, demo$analyte == "Cu")
#>             
#>              FALSE TRUE
#>   downstream   102    6
#>   leachate     102    6
#>   reference    102    6
head(subset(demo, site_id == "downstream"), 4)
#> # A tibble: 4 × 9
#>   sample_id site_id    datetime   analyte value detected units.analyte
#>   <chr>     <chr>      <date>     <chr>   <dbl> <lgl>    <chr>        
#> 1 DS-01     downstream 2024-01-15 Na      149.  TRUE     mg/L         
#> 2 DS-01     downstream 2024-01-15 K        30.0 TRUE     mg/L         
#> 3 DS-01     downstream 2024-01-15 Ca       75.5 TRUE     mg/L         
#> 4 DS-01     downstream 2024-01-15 Mg       25.9 TRUE     mg/L         
#> # ℹ 2 more variables: valence.analyte <dbl>, atomic_mass.analyte <dbl>
```
