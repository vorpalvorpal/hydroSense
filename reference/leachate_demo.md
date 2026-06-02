# Synthetic leachate-impacted water-quality data (demo)

A small, fictional long-format water-quality dataset used in examples
and the package vignettes. It is **not** real monitoring data —
concentrations are generated from a simple two-component mixing model so
the worked examples are reproducible and self-contained.

## Usage

``` r
leachate_demo
```

## Format

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
  [`to_meq()`](https://www.kedumba.com.au/leachatetools/reference/to_meq.md));
  `NA` for non-ionic analytes.

- atomic_mass.analyte:

  Numeric. Molar/atomic mass in g/mol (for the meq conversion); `NA` for
  non-ionic analytes.

## Details

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
[`add_lmf()`](https://www.kedumba.com.au/leachatetools/reference/add_lmf.md)
/
[`to_meq()`](https://www.kedumba.com.au/leachatetools/reference/to_meq.md)
(`Na`, `K`, `Ca`, `Mg`, `Cl`, `SO4`-with-charge, `F`, `NO3-N`, `NO2-N`,
`CO3-CaCO3`, `HCO3-CaCO3`), the toxicants assessed by
[`ssd_paf()`](https://www.kedumba.com.au/leachatetools/reference/ssd_paf.md)
/
[`add_amspaf()`](https://www.kedumba.com.au/leachatetools/reference/add_amspaf.md)
(`Cu`, `Zn`, `NH3-N`), and the co-analytes needed for the ammonia and
bioavailability normalisations (`pH`, `temperature`, `DOC`, `hardness`).

## Examples

``` r
# Leachate-mixing fraction per downstream sample:
# \donttest{
add_lmf(
  df             = subset(leachate_demo, site_id == "downstream"),
  leachate_data  = subset(leachate_demo, site_id == "leachate"),
  reference_data = subset(leachate_demo, site_id == "reference")
)
#> # A tibble: 114 × 33
#>    sample_id site_id    datetime   analyte     value detected units.analyte
#>    <chr>     <chr>      <date>     <chr>       <dbl> <lgl>    <chr>        
#>  1 DS-01     downstream 2024-01-15 Na        149.    TRUE     mg/L         
#>  2 DS-01     downstream 2024-01-15 K          30.0   TRUE     mg/L         
#>  3 DS-01     downstream 2024-01-15 Ca         75.5   TRUE     mg/L         
#>  4 DS-01     downstream 2024-01-15 Mg         25.9   TRUE     mg/L         
#>  5 DS-01     downstream 2024-01-15 Cl        295.    TRUE     mg/L         
#>  6 DS-01     downstream 2024-01-15 SO4²⁻      54.2   TRUE     mg/L         
#>  7 DS-01     downstream 2024-01-15 F           0.859 TRUE     mg/L         
#>  8 DS-01     downstream 2024-01-15 NO3-N       1.57  TRUE     mg/L         
#>  9 DS-01     downstream 2024-01-15 NO2-N       0.155 TRUE     mg/L         
#> 10 DS-01     downstream 2024-01-15 CO3-CaCO3  16.3   TRUE     mg/L         
#> # ℹ 104 more rows
#> # ℹ 26 more variables: valence.analyte <dbl>, atomic_mass.analyte <dbl>,
#> #   lmf_naive <dbl>, lmf_reason <chr>, n_ions_used <int>,
#> #   n_ions_downweighted <int>, sigma_lmf <dbl>, chi2_per_df <dbl>,
#> #   uuid.analyte <chr>, uuid <chr>, value.guideline_1 <dbl>,
#> #   level_name.guideline_1 <chr>, guideline.guideline_1 <chr>,
#> #   comments.guideline_1 <glue>, value.guideline_2 <dbl>, …
# }
```
